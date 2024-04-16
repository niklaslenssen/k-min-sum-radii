// k_MSR.h

#include "Ball.h"
#include "Cluster.h"
#include "Point.h"

#ifndef K_MSR_H
#define K_MSR_H

bool containsPoint(const Point& p, const std::vector<Ball>& balls);

std::vector<Ball> selection(const std::vector<Point>& points, int k,
                            const std::vector<int>& u,
                            const std::vector<double>& radii, double epsilon);

std::vector<std::vector<double>> getRadii(double rmax, int k, double epsilon);

double logBase(double x, double b);

std::vector<std::vector<int>> getU(int k, double epsilon);

bool containsAllPoints(const std::vector<Point>& points,
                       const std::vector<Ball>& balls);

double cost(const std::vector<Cluster>& cluster);

std::vector<Cluster> clustering(const std::vector<Point>& points, int k,
                                double epsilon, double rmax);

#endif