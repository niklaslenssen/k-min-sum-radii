#include "header/Cluster.h"
#include "header/Point.h"

#include <random>

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

vector<Cluster> gonzales(vector<Point> &points, int k) {
  int n = points.size();

  vector<Point> centers;
  centers.push_back(points[0]); // Initialisiere das erste Zentrum als ersten Punkt

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

  // Erstelle Cluster basierend auf den Zentren
  vector<Cluster> clusters = assignPointsToCluster(points, centers, k);
  return clusters;
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
  return clusters;
}