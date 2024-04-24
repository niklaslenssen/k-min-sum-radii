// welzl.h

#include "Ball.h"
#include "Point.h"

#ifndef WELZL_H
#define WELZL_H

Ball welzlHelper(std::vector<Point>& points, std::vector<Point> support, int n);
Ball makeBallThreePoints(const Point& A, const Point& B, const Point& C);
Ball makeBallTwoPoints(const Point& A, const Point& B);
Ball findMinEnclosingBall(std::vector<Point>& points);

#endif