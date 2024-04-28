#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "header/Cluster.h"
#include "header/badoiu_clarkson.h"
#include "header/hochbaumShmyos.h"
#include "header/k_MSR.h"
#include "header/welzl.h"

using namespace std;

// Schreibt die übergebenen Bälle in eine CSV Datei.
void saveBallsInCSV(const vector<Ball> &balls) {
  ofstream ballsfile("Data/balls.csv");

  if (!ballsfile) {
    cerr << "Fehler beim Öffnen der Datei!" << endl;
    exit;
  }

  for (Ball b : balls) {
    ballsfile << b.center.coordinates[0] << "," << b.center.coordinates[1] << "," << b.radius << "\n";
  }

  ballsfile.close();
  return;
}

// Schreibt die übergebenen Cluster in eine CSV Datei.
void saveClusterInCSV(const vector<Cluster> &cluster) {
  ofstream clusterfile("Data/cluster.csv");

  if (!clusterfile) {
    cerr << "Fehler beim Öffnen der Datei!" << endl;
    exit;
  }

  for (int i = 0; i < cluster.size(); i++) {
    Cluster c = cluster[i];
    for (int j = 0; j < c.points.size(); j++) {
      Point p = c.points[j];
      clusterfile << p.coordinates[0] << "," << p.coordinates[1] << "," << i << endl;
    }
  }

  clusterfile.close();
  return;
}

// Liest eine Menge von Punkten und rmax von einer CSV Datei ein.
vector<Point> readPointsFromCSV(int k) {
  vector<Point> points;
  ifstream file("Data/points.csv");
  string line;

  if (!file.is_open()) {
    cerr << "Datei konnte nicht geöffnet werden!" << endl;
    exit;
  }

  while (getline(file, line)) {
    vector<double> coords;
    stringstream stream(line);
    string value;

    getline(stream, value, ',');

    while (getline(stream, value, ',')) {
      coords.push_back(stod(value));
    }

    points.push_back(Point(coords));
  }

  file.close();

  return points;
}

int main(int argc, char const *argv[]) {
  if (argc != 3) {
    cerr << "Bitte übergebe einen Wert für 'k' und 'epsilon'!" << endl;
    return -1;
  }

  int k = stod(argv[1]);
  double epsilon = stod(argv[2]);

  vector<Point> points = readPointsFromCSV(k);

  double rmax = hochbaumShmoysKCenter(points, k);

  auto start = std::chrono::steady_clock::now();

  vector<Cluster> cluster = clustering(points, k, epsilon, rmax);

  auto end = std::chrono::steady_clock::now();

  vector<Ball> balls;

  for (int i = 0; i < cluster.size(); i++) {
    Ball b = findMinEnclosingBall(cluster[i].points);
    if (b.radius != 0) {
      balls.push_back(b);
    }
  }

  double radii = 0;
  for (Ball b : balls) {
    radii += b.radius;
  }


  auto diff = end - start;

  cout << "Dauer des Durchlaufs: " << chrono::duration<double>(diff).count() << " Sekunden" << endl;
  cout << "Schmidt:                   " << radii << endl;

  saveClusterInCSV(cluster);
  saveBallsInCSV(balls);

  return 0;
}