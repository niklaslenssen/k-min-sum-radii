#include <algorithm>
#include <omp.h>
#include <random>

#include "header/Ball.h"
#include "header/Cluster.h"
#include "header/Point.h"
#include "header/k_MSR.h"
#include "header/welzl.h"
#include "header/yildirim.h"

using namespace std;

vector<Cluster> assignPointsToCluster(vector<Point> &points, vector<Point> &centers, int k) {
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

    clusters[closestCenter].addPoint(points[i]);
  }
  return clusters;
}

bool clustersOverlap(Cluster &c1, Cluster &c2) {
  Ball b1 = findMinEnclosingBall(c1.points);
  Ball b2 = findMinEnclosingBall(c2.points);

  double distance = Point::distance(b1.center, b2.center);
  double radiusSum = b1.radius + b2.radius;
  return distance <= radiusSum;
}

vector<Cluster> mergeCluster(vector<Cluster> &clusters) {
  bool changed;

  do {
    changed = false;
    vector<Cluster> mergedClusters;
    vector<bool> merged(clusters.size(), false);

    for (int i = 0; i < clusters.size(); i++) {
      if (merged[i]) {
        continue;
      }
      Cluster currentCluster = clusters[i];
      merged[i] = true;

      for (int j = i + 1; j < clusters.size(); j++) {
        if (merged[j]) {
          continue;
        }
        if (clustersOverlap(currentCluster, clusters[j])) {
          currentCluster.merge(clusters[j]);
          merged[j] = true;
          changed = true;
        }
      }
      mergedClusters.push_back(currentCluster);
    }

    clusters = mergedClusters;

  } while (changed);

  return clusters;
}

vector<Cluster> gonzales(vector<Point> &points, int k) {
  int n = points.size();
  vector<Point> centers;
  centers.push_back(points[0]);

  for (int i = 1; i < k; i++) {
    int nextCenter = -1;
    double maxDist = -1.0;

    // Finde den Punkt, der am weitesten von seinem nächsten Zentrum entfernt ist
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
  random_device rd;
  mt19937 gen(rd());
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

    // Berechne die Wahrscheinlichkeitsverteilung für die Auswahl des nächsten Zentrums
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

  // Erstelle Cluster basierend auf den Zentren
  vector<Cluster> clusters = assignPointsToCluster(points, centers, k);
  return mergeCluster(clusters);
}

vector<Cluster> heuristik(vector<Point> &points, int k) {
  int n = points.size();
  vector<Cluster> bestCluster;
  bestCluster.push_back(Cluster(points));
  vector<vector<double>> distances(n, vector<double>(n, 0));

#pragma omp parallel for collapse(2)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      distances[i][j] = Point::distance(points[i], points[j]);
    }
  }

#pragma omp parallel
  {
    vector<Cluster> localBestCluster = bestCluster;
    double localBestCost = cost(localBestCluster);

#pragma omp for
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        vector<Point> centers;
        Point largestCenter = points[i];
        centers.push_back(largestCenter);
        double radius = distances[i][j];

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

        vector<Cluster> cluster = assignPointsToCluster(points, centers, k);
        double clusterCost = cost(cluster);

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

  return mergeCluster(bestCluster);
}