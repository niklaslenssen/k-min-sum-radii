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
#include "header/heuristic.h"
#include "header/hochbaumShmyos.h"
#include "header/k_MSR.h"
#include "header/welzl.h"

using namespace std;

// Schreibt die übergebenen Bälle in eine CSV Datei.
void saveBallsInCSV(const vector<Ball> &balls, const string &filePath) {
  ofstream ballsfile(filePath);
  cout << filePath << endl;

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
void saveClusterInCSV(const vector<Cluster> &cluster, const string &filePath) {
  ofstream clusterfile(filePath);

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
vector<Point> readPointsFromCSV(const string &filePath) {
  vector<Point> points;
  ifstream file(filePath);
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
  if (argc != 7) {
    cerr << "Bitte übergebe einen Wert für 'k', 'epsilon' und die Anzahl an 'u'!" << endl;
    return -1;
  }

  int k = stod(argv[1]);
  double epsilon = stod(argv[2]);
  int numVectors = stod(argv[3]);
  string pointFilePath = argv[4];
  string ballFilePath = argv[5];
  string clusterFilePath = argv[6];

  vector<Point> points = readPointsFromCSV(pointFilePath);

  double rmax = hochbaumShmoysKCenter(points, k);

  auto start = std::chrono::steady_clock::now();

  vector<Cluster> cluster = clustering(points, k, epsilon, rmax, numVectors);

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

  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);

  return 0;
}