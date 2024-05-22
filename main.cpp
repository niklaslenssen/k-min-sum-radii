#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "header/Cluster.h"
#include "header/Heuristic.h"
#include "header/badoiu_clarkson.h"
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
    exit(1);
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
    exit(1);
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
    exit(1);
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

// Erzeugt die MEBs der übergeben Cluster.
vector<Ball> getBallsFromCluster(vector<Cluster> &cluster) {
  vector<Ball> balls;
  for (int i = 0; i < cluster.size(); i++) {
    Ball b = findMinEnclosingBall(cluster[i].points);
    if (b.radius != 0) {
      balls.push_back(b);
    }
  }
  return balls;
}

// Berechnet die Summe der Radien einer Liste von Bällen.
double sumOfRadii(vector<Ball> &balls) {
  double radii = 0;
  for (Ball b : balls) {
    radii += b.radius;
  }
  return radii;
}

void analyseSchmidt(vector<Point> &points, int k, double epsilon, int numVectors, string clusterFilePath, string ballFilePath) {
  double rmax = hochbaumShmoysKCenter(points, k);
  auto start = std::chrono::steady_clock::now();
  vector<Cluster> cluster = clustering(points, k, epsilon, rmax, numVectors);
  auto end = std::chrono::steady_clock::now();
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  auto diff = end - start;
  cout << "Dauer des Durchlaufs: " << chrono::duration<double>(diff).count() << " Sekunden" << endl;
  cout << "Schmidt:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

void analyseGonzales(vector<Point> &points, int k, string clusterFilePath, string ballFilePath) {
  vector<Cluster> cluster = gonzales(points, k);
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  cout << "Gonzales:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

int main(int argc, char const *argv[]) {
  if (string(argv[1]) == "s") {
    int k = stod(argv[2]);
    double epsilon = stod(argv[3]);
    int numVectors = stod(argv[4]);
    string pointFilePath = argv[5];
    string ballFilePath = argv[6];
    string clusterFilePath = argv[7];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseSchmidt(points, k, epsilon, numVectors, clusterFilePath, ballFilePath);
  } else if (string(argv[1]) == "g") {
    int k = stod(argv[2]);
    string pointFilePath = argv[3];
    string ballFilePath = argv[4];
    string clusterFilePath = argv[5];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseGonzales(points, k, clusterFilePath, ballFilePath);
  }
  return 0;
}