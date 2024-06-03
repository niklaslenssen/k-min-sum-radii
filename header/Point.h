// Point.h

#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <vector>
#include <stdexcept>

class Point {
 public:
  std::vector<double> coordinates;

  Point() {}
  Point(std::vector<double> coords) : coordinates(coords) {}

  // Euklidischer Abstand zwischen 'this' und einem anderen Punkt.
  double distanceTo(const Point& other) const {
    if (coordinates.size() != other.coordinates.size()) {
      throw std::invalid_argument("Die Punkte müssen dieselbe Dimension haben!");
    }
    double sum = 0.0;
    for (int i = 0; i < coordinates.size(); i++) {
      sum += (coordinates[i] - other.coordinates[i]) * (coordinates[i] - other.coordinates[i]);
    }
    return sqrt(sum);
  }

  // Euklidischer Abstand zwischen zwei Punkten.
  static double distance(const Point& p1, const Point& p2) {
    return p1.distanceTo(p2);
  }

  bool operator<(const Point& other) const {
    return coordinates < other.coordinates;
  }

  Point operator+(const Point& other) const {
    if (coordinates.size() != other.coordinates.size()) {
      throw std::invalid_argument("Die Punkte müssen dieselbe Dimension haben!");
    }
    std::vector<double> result_coords(coordinates.size());
    for (int i = 0; i < coordinates.size(); i++) {
      result_coords[i] = coordinates[i] + other.coordinates[i];
    }
    return Point(result_coords);
  }

  Point operator*(double scalar) const {
    std::vector<double> result_coords(coordinates.size());
    for (int i = 0; i < coordinates.size(); i++) {
      result_coords[i] = coordinates[i] * scalar;
    }
    return Point(result_coords);
  }
};

#endif