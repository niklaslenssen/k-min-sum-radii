#include <algorithm>
#include <vector>

#include "header/Ball.h"
#include "header/Point.h"
#include "header/welzl.h"

using namespace std;

Ball makeBallTwoPoints(const Point &A, const Point &B) {
  double x = (A.coordinates[0] + B.coordinates[0]) / 2;
  double y = (A.coordinates[1] + B.coordinates[1]) / 2;
  vector<double> centerCoords = {x, y};
  Point center(centerCoords);
  double radius = A.distanceTo(B) / 2;
  return Ball(center, radius);
}

Ball makeBallThreePoints(const Point &A, const Point &B, const Point &C) {
  double ax = A.coordinates[0], ay = A.coordinates[1];
  double bx = B.coordinates[0], by = B.coordinates[1];
  double cx = C.coordinates[0], cy = C.coordinates[1];

  double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  if (std::abs(d) < 1e-10)
    return Ball(Point({0, 0}), 0);

  double x = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
  double y = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
  Point center({x, y});
  double radius = center.distanceTo(A);
  return Ball(center, radius);
}

Ball welzlHelper(vector<Point> &points, vector<Point> support, int n) {
  if (n == 0 || support.size() == 3) {
    switch (support.size()) {
    case 0:
      return Ball();
    case 1:
      return Ball(support[0], 0);
    case 2:
      return makeBallTwoPoints(support[0], support[1]);
    case 3:
      return makeBallThreePoints(support[0], support[1], support[2]);
    }
  }

  int idx = rand() % n;
  Point p = points[idx];
  std::swap(points[idx], points[n - 1]);

  Ball d = welzlHelper(points, support, n - 1);

  if (!d.contains(p)) {
    support.push_back(p);
    d = welzlHelper(points, support, n - 1);
    support.pop_back();
  }

  std::swap(points[idx], points[n - 1]);
  return d;
}

Ball findMinEnclosingBall(vector<Point> &points) {
  srand(static_cast<unsigned int>(time(nullptr)));
  vector<Point> support;
  std::random_shuffle(points.begin(), points.end());
  return welzlHelper(points, support, points.size());
}