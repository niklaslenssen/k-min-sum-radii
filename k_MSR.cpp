#include "header/k_MSR.h"

#include <omp.h>

#include <iostream>
#include <random>
#include <set>

#include "header/MEB.h"

using namespace std;

double logBase(double x, double b) { return log(x) / log(b); }

bool containsAllPoints(const vector<Point> &points, const vector<Ball> &balls) {
  int count = 0;

  for (Point p : points) {
    for (Ball b : balls) {
      if (b.contains(p)) {
        count++;
        break;
      }
    }
  }

  return (count == points.size());
}

bool containsPoint(const Point &p, const vector<Ball> &balls) {
  for (Ball b : balls) {
    if (b.contains(p)) {
      return true;
    }
  }
  return false;
}

double cost(const vector<Cluster> &cluster) {
  double result = 0;

  for (Cluster c : cluster) {
    result += findMinEnclosingBall(c.points).radius;
  }
  return result;
}

vector<vector<double>> getRadii(double rmax, int k, double epsilon) {
  vector<vector<double>> result;
  vector<int> indices(k - 1, 0);
  vector<double> set;

  int limit = ceil(logBase((k / epsilon), (1 + epsilon)));

  for (int i = 0; i <= limit; i++) {
    set.push_back(pow((1 + epsilon), i) * (epsilon / k) * rmax);
  }

  while (true) {
    vector<double> current;

    current.push_back(rmax);

    for (int idx : indices) {
      current.push_back(set[idx]);
    }
    result.push_back(current);

    int next = k - 2;
    while (next >= 0 && ++indices[next] == set.size()) {
      indices[next] = 0;
      next--;
    }
    if (next < 0) break;
  }

  return result;
}

vector<vector<int>> getU(int k, double epsilon) {
  int length = 8;
  vector<vector<int>> result;
  vector<int> indices(length, 0);
  vector<int> set;

  for (int i = 0; i < k; i++) {
    set.push_back(i);
  }

  while (true) {
    vector<int> current;

    for (int idx : indices) {
      current.push_back(set[idx]);
    }
    result.push_back(current);

    int next = length - 1;
    while (next >= 0 && ++indices[next] == set.size()) {
      indices[next] = 0;
      next--;
    }
    if (next < 0) break;
  }

  return result;
}

vector<vector<int>> getRandomU(int k, double epsilon) {
  vector<vector<int>> result;
  int length = (32 * k * (1 + epsilon)) / (pow(epsilon, 3));
  int numVectors = 100;
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distrib(0, k - 1);

  for (int i = 0; i < numVectors; i++) {
    vector<int> currentVector;
    for (int j = 0; j < length; j++) {
      currentVector.push_back(distrib(gen));
    }
    result.push_back(currentVector);
  }

  return result;
}

vector<Ball> selection(const vector<Point> &points, int k, const vector<int> &u,
                       const vector<double> &radii, double epsilon) {
  vector<Ball> balls(k);

  vector<vector<Point>> Si(k);

  double lambda = 1 + epsilon + 2 * sqrt(epsilon);

  set<Point> X;

  vector<Ball> R;

  for (int j = 0; j < Si.size(); j++) {
    if (Si[j].size() == 1) {
      Point center = Si[j][0];
      double radius = (epsilon / (1 + epsilon)) * radii[j];
      R.push_back(Ball(center, radius));
    }
  }

  for (int i = 0; i < u.size(); i++) {
    
    if (X.size() >= points.size()) {
      return balls;
    }

    vector<Ball> Ri = R;

    for (Point p : points) {
      if ((X.count(p) == 0) && !containsPoint(p, Ri)) {
        Si[u[i]].push_back(p);
        break;
      }
    }

    if (Si[u[i]].size() >= 2) {
      Ball b = findMinEnclosingBall(Si[u[i]]);
      b.radius = b.radius * lambda;
      balls[u[i]] = b;

      for (Point p : points) {
        for (Ball b : balls) {
          if (b.center.coordinates.size() > 0) {
            if (b.contains(p)) {
              X.insert(p);
              break;
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < Si.size(); i++) {
    if (Si[i].size() == 1) {
      balls[i] = Ball(Si[i][0], 0);
    }
  }

  return balls;
}

vector<Cluster> clustering(const vector<Point> &points, int k, double epsilon,
                           double rmax) {
  vector<Cluster> cluster(k);
  vector<Cluster> bestCluster(k);
  vector<Ball> balls(k);

  vector<vector<double>> radii = getRadii(rmax, k, epsilon);
  vector<vector<int>> u = getRandomU(k, epsilon);

  cout << u.size() << " " << u[0].size() << endl;
  cout << radii.size() << " " << radii[0].size() << endl;

  bestCluster[0].points = points;

  double bestCost = 10000;

#pragma omp parallel for collapse(2) schedule(dynamic) \
    shared(bestCluster, bestCost)
  for (int i = 0; i < radii.size(); i++) {
    for (int j = 0; j < u.size(); j++) {
      vector<double> r = radii[i];
      vector<int> ui = u[j];

      vector<Ball> localBalls = selection(points, k, ui, r, epsilon);

      if (containsAllPoints(points, localBalls)) {
        vector<Cluster> localCluster(k);

        for (Point p : points) {
          for (int c = 0; c < k; c++) {
            if (localBalls[c].contains(p)) {
              localCluster[c].points.push_back(p);
              break;
            }
          }
        }

        double localCost = cost(localCluster);

#pragma omp critical
        {
          if (localCost < bestCost) {
            bestCost = localCost;
            bestCluster = localCluster;
          }
        }
      }
    }
  }

  return bestCluster;
}