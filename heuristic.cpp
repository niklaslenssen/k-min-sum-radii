#include "header/Cluster.h"
#include "header/Point.h"

using namespace std;

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
  vector<Cluster> clusters(k);
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
