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

  double distanceTo(const Point& other) const {
    if (coordinates.size() != other.coordinates.size()) {
      throw std::invalid_argument("Die Punkte müssen dieselbe Dimension haben!");
    }
    double sum = 0.0;
    for (int i = 0; i < coordinates.size(); ++i) {
      sum += (coordinates[i] - other.coordinates[i]) *
             (coordinates[i] - other.coordinates[i]);
    }
    return sqrt(sum);
  }

  bool operator<(const Point& other) const {
    return coordinates < other.coordinates;
  }
};

#endif