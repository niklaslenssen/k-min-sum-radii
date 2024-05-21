// Cluster.h

#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

#include "Point.h"

class Cluster {
 public:
  std::vector<Point> points;

  void addPoint(const Point& p) {
    points.push_back(p);
  }
};

#endif