// Heuristic.h

#ifndef HEURISTIC_H
#define HEURISTIC_H

#include "Point.h"
#include "Cluster.h"
#include "Ball.h"

std::vector<Cluster> gonzales(std::vector<Point> &points, int k);
std::vector<Cluster> kMeansPlusPlus(std::vector<Point> &points, int k);
Ball heuristik(std::vector<Point> &points, int k);

#endif