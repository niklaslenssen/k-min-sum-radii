// Ball.h

#ifndef BALL_H
#define BALL_H

#include "Point.h"

class Ball {
 public:
  Point center;
  double radius;

  Ball() : center({0, 0}), radius(0.0) {}
  Ball(const Point& center, double radius) : center(center), radius(radius) {}

  bool contains(const Point& p) const { return center.distanceTo(p) <= radius; }
};

#endif