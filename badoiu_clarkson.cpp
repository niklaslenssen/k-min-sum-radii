#include "header/badoiu_clarkson.h"

using namespace std;

Ball findMEB(const vector<Point> &points, double epsilon) {
  // Starte mit dem ersten Punkt als Mittelpunkt und einem sehr kleinen Radius
  Point center = points[0];
  double radius = 0.0;

  while (true) {
    // Finde den am weitesten entfernten Punkt vom Zentrum
    Point farthest_point = Point({0, 0});
    double max_distance = 0.0;
    for (const Point &p : points) {
      double dist = Point::distance(center, p);
      if (dist > radius) {
        radius = dist;
        farthest_point = p;
      }
    }

    // Prüfe, ob der gefundene Punkt innerhalb des aktuellen Balls liegt
    double dist_to_farthest = Point::distance(center, farthest_point);
    if (dist_to_farthest <= radius) {
      break; // Der aktuelle Ball umschließt alle Punkte
    }

    // Aktualisiere den Mittelpunkt und Radius des Balls
    for (int i = 0; i < center.coordinates.size(); i++) {
      center.coordinates[i] += (farthest_point.coordinates[i] - center.coordinates[i]) * (epsilon / 2.0);
    }
    radius = (1 + epsilon) * radius;
  }

  return Ball(center, radius);
}