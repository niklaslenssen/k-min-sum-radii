#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "header/Cluster.h"
#include "header/MEB.h"
#include "header/k_MSR.h"

using namespace std;

void saveBallsInCSV(const vector<Ball>& balls) {
  ofstream ballsfile("Data/balls.csv");

  if (!ballsfile) {
    cerr << "Fehler beim Öffnen der Datei!" << endl;
    exit;
  }

  for (Ball b : balls) {
    ballsfile << b.center.coordinates[0] << "," << b.center.coordinates[1]
              << "," << b.radius << "\n";
  }

  ballsfile.close();
  return;
}

void saveClusterInCSV(const vector<Cluster>& cluster) {
  ofstream clusterfile("Data/cluster.csv");

  if (!clusterfile) {
    cerr << "Fehler beim Öffnen der Datei!" << endl;
    exit;
  }

  for (int i = 0; i < cluster.size(); i++) {
    Cluster c = cluster[i];
    for (int j = 0; j < c.points.size(); j++) {
      Point p = c.points[j];
      clusterfile << p.coordinates[0] << "," << p.coordinates[1] << "," << i
                  << endl;
    }
  }

  clusterfile.close();
  return;
}

vector<Point> readPointsFromCSV(double& rmax) {
  vector<Point> points;
  ifstream file("Data/points.csv");
  string line;

  if (!file.is_open()) {
    cerr << "Datei konnte nicht geöffnet werden!" << endl;
    exit;
  }

  getline(file, line);

  rmax = stod(line);

  while (getline(file, line)) {
    vector<double> coords;
    stringstream stream(line);
    string value;

    while (getline(stream, value, ',')) {
      coords.push_back(stod(value));
    }

    coords.pop_back();

    points.push_back(Point(coords));
  }

  file.close();

  return points;
}

int main(int argc, char const* argv[]) {
  int k = 3;
  double epsilon = 0.5;
  double rmax;

  vector<Point> points = readPointsFromCSV(rmax);

  auto start = std::chrono::steady_clock::now();

  vector<Cluster> cluster = clustering(points, k, epsilon, rmax);

  auto end = std::chrono::steady_clock::now();

  vector<Ball> balls(k);

  for (int i = 0; i < cluster.size(); i++) {
    balls[i] = findMinEnclosingBall(cluster[i].points);
  }

  double radii = 0;
  for (Ball b : balls) {
    radii += b.radius;
  }

  cout << "Summe der Radien: " << radii << endl;

  auto diff = end - start;

  cout << "Dauer des Durchlaufs: " << chrono::duration<double>(diff).count()
       << " Sekunden" << endl;

  saveClusterInCSV(cluster);
  saveBallsInCSV(balls);

  return 0;
}