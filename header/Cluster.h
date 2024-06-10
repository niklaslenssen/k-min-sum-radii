// Cluster.h

#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

#include "Point.h"

class Cluster {
 public:
  std::vector<Point> points;

  Cluster() {}
  Cluster(const std::vector<Point>& points) : points(points) {}

  void addPoint(const Point& p) {
    points.push_back(p);
  }

  void merge(const Cluster& other) {
    points.insert(points.end(), other.points.begin(), other.points.end());
  }
};

#endif