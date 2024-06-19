#include <omp.h>

#include <random>

#include "header/ball.h"
#include "header/cluster.h"
#include "header/k_MSR.h"
#include "header/point.h"
#include "header/welzl.h"
#include "header/yildirim.h"

using namespace std;

// Weist jedem Punkt im Vektor 'points' das nächstgelegene Zentrum im Vektor
// 'centers' zu
vector<Cluster> assignPointsToCluster(const vector<Point> &points,
                                      const vector<Point> &centers, int k) {
  int n = points.size();
  vector<Cluster> clusters(k);

  // Erstelle Cluster basierend auf den Zentren
  for (int i = 0; i < n; i++) {
    int closestCenter = -1;
    double minDist = numeric_limits<double>::max();

    // Finde das nächstgelegene Zentrum für jeden Punkt
    for (int j = 0; j < k; j++) {
      double dist = points[i].distanceTo(centers[j]);
      if (dist < minDist) {
        minDist = dist;
        closestCenter = j;
      }
    }
    // Füge den aktuellen Punkt zu seinem nächsten Cluster hinzu
    clusters[closestCenter].addPoint(points[i]);
  }
  return clusters;
}

// Überprüft, ob zwei Cluster sich überlappen oder berühren
bool clustersOverlap(const Cluster &c1, const Cluster &c2) {
  // Berechne die MEBs
  Ball b1 = findMinEnclosingBall(c1.getPoints());
  Ball b2 = findMinEnclosingBall(c2.getPoints());

  // Berechne die euklidische Distanz zwischen den Zentren der beiden Bälle
  double distance = Point::distance(b1.getCenter(), b2.getCenter());

  // Berechne die Summe der Radien der beiden Bälle
  double radiusSum = b1.getRadius() + b2.getRadius();

  // Überprüfe, ob die Distanz zwischen den Zentren kleiner oder gleich der
  // Summe der Radien ist
  return distance <= radiusSum;
}

// Merged überlappende oder berührende Cluster
vector<Cluster> mergeCluster(vector<Cluster> &clusters) {
  bool changed;

  // Wiederhole den Merge-Vorgang, bis keine Cluster mehr gemerged werden
  do {
    changed = false;
    vector<Cluster> mergedClusters;
    vector<bool> merged(clusters.size(), false);

    for (int i = 0; i < clusters.size(); i++) {
      if (merged[i]) {
        continue;  // Überspringe bereits gemergte Cluster
      }
      Cluster currentCluster = clusters[i];
      merged[i] = true;

      for (int j = i + 1; j < clusters.size(); j++) {
        if (merged[j]) {
          continue;  // Überspringe bereits gemergte Cluster
        }
        if (clustersOverlap(currentCluster, clusters[j])) {
          currentCluster.merge(clusters[j]);
          merged[j] = true;
          changed = true;  // Es gab eine Änderung
        }
      }
      mergedClusters.push_back(currentCluster);  // Füge das gemergte Cluster zu
                                                 // den gemergten Clustern hinzu
    }

    clusters = mergedClusters;  // Aktualisiere die Cluster-Liste

  } while (changed);

  return clusters;
}

// Berechnet den Schwerpunkt (Centroid) des Clusters
Point computeCentroid(const vector<Point> &points) {
  int dimension = points[0].getCoordinates().size();
  vector<double> centroidCoords(dimension, 0.0);

  // Summe der Koordinaten aller Punkte im Cluster
  for (const Point &p : points) {
    for (int i = 0; i < dimension; i++) {
      centroidCoords[i] += p.getCoordinates()[i];
    }
  }

  // Mittelwert der Koordinaten berechnen
  for (int i = 0; i < dimension; i++) {
    centroidCoords[i] /= points.size();
  }

  return Point(centroidCoords);
}

vector<Cluster> gonzales(vector<Point> &points, int k) {
  srand(1234);
  int n = points.size();
  vector<Point> centers;
  centers.push_back(points[rand() % n]);

  // Finde die restlichen k-1 Zentren
  for (int i = 1; i < k; i++) {
    int nextCenter = -1;
    double maxDist = -1.0;

    // Finde den Punkt, der am weitesten von seinem nächsten Zentrum entfernt
    // ist
    for (int j = 0; j < n; j++) {
      double dist = numeric_limits<double>::max();
      for (Point center : centers) {
        dist = min(dist, points[j].distanceTo(center));
      }
      if (dist > maxDist) {
        maxDist = dist;
        nextCenter = j;
      }
    }
    centers.push_back(points[nextCenter]);
  }

  // Weise die Punkte den Zentren zu und erstelle Cluster
  vector<Cluster> clusters = assignPointsToCluster(points, centers, k);

  // Merge überlappende oder berührende Cluster
  return mergeCluster(clusters);
}

vector<Cluster> kMeansPlusPlus(vector<Point> &points, int k) {
  int n = points.size();
  vector<Point> centers;
  mt19937 gen(1234);
  uniform_int_distribution<> dis(0, n - 1);

  // Wähle das erste Zentrum zufällig aus
  centers.push_back(points[dis(gen)]);

  // Wähle die restlichen Zentren basierend auf der Distanzverteilung
  for (int i = 1; i < k; i++) {
    vector<double> dist(n, numeric_limits<double>::max());

    for (int j = 0; j < n; j++) {
      for (const Point &center : centers) {
        dist[j] = min(dist[j], points[j].distanceTo(center));
      }
    }

    // Berechne die Wahrscheinlichkeitsverteilung für die Auswahl des nächsten
    // Zentrums
    vector<double> distSquared(n);
    double sumDist = 0.0;
    for (int j = 0; j < n; j++) {
      distSquared[j] = dist[j] * dist[j];
      sumDist += distSquared[j];
    }

    uniform_real_distribution<> disReal(0, sumDist);
    double r = disReal(gen);
    double cumulativeDist = 0.0;

    for (int j = 0; j < n; j++) {
      cumulativeDist += distSquared[j];
      if (cumulativeDist >= r) {
        centers.push_back(points[j]);
        break;
      }
    }
  }

  vector<Cluster> clusters;
  bool changed = true;
  while (changed) {
    changed = false;

    // Weise die Punkte den Zentren zu und erstelle Cluster
    clusters = assignPointsToCluster(points, centers, k);

    // Aktualisiere die Zentren basierend auf den Clustern
    for (int i = 0; i < k; i++) {
      Point newCenter = computeCentroid(clusters[i].getPoints());
      if (newCenter != centers[i]) {
        centers[i] = newCenter;
        changed = true;
      }
    }
  }

  // Merge überlappende oder berührende Cluster
  return mergeCluster(clusters);
}

vector<Cluster> heuristik(vector<Point> &points, int k) {
  int n = points.size();
  vector<Cluster> bestCluster;
  bestCluster.push_back(
      Cluster(points));  // Initialisiere mit allen Punkten in einem Cluster
  vector<vector<double>> distances(n, vector<double>(n, 0));

// Berechnung der Abstände zwischen allen Punkten
#pragma omp parallel for collapse(2)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      distances[i][j] = Point::distance(points[i], points[j]);
    }
  }

#pragma omp parallel
  {
    vector<Cluster> localBestCluster =
        bestCluster;  // Lokale Variable für die besten Cluster in jedem Thread
    double localBestCost =
        cost(localBestCluster);  // Kosten der lokalen besten Cluster

#pragma omp for
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        vector<Point> centers;
        Point largestCenter = points[i];
        centers.push_back(largestCenter);
        double radius = distances[i][j];

        // Finde k Zentren
        while (centers.size() != k) {
          int nextCenter = -1;
          double maxDist = -1.0;

          for (int h = 0; h < n; h++) {
            double dist = numeric_limits<double>::max();
            for (const Point &center : centers) {
              if (center == largestCenter) {
                dist = min(dist, points[h].distanceTo(center) - radius);
              } else {
                dist = min(dist, points[h].distanceTo(center));
              }
            }
            if (dist > maxDist) {
              maxDist = dist;
              nextCenter = h;
            }
          }
          centers.push_back(points[nextCenter]);
        }

        // Weise die Punkte den nächstgelegenen Zentren zu
        vector<Cluster> cluster = assignPointsToCluster(points, centers, k);
        double clusterCost =
            cost(cluster);  // Berechne die Kosten des aktuellen Clusters

        // Aktualisiere lokale beste Cluster, falls das aktuelle Cluster besser
        // ist
        if (clusterCost < localBestCost) {
          localBestCluster = cluster;
          localBestCost = clusterCost;
        }
      }
    }
#pragma omp critical
    {
      if (localBestCost < cost(bestCluster)) {
        bestCluster = localBestCluster;
      }
    }
  }

  // Merged überlappende oder berührende Cluster
  return mergeCluster(bestCluster);
}
