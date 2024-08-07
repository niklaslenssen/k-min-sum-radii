#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <vector>

#include "kMSR/header/cluster.h"
#include "kMSR/header/heuristic.h"
#include "kMSR/header/k_MSR.h"
#include "kMSR/header/welzl.h"
#include "kMSR/header/yildirim.h"

using namespace std;

// Schreibt die übergebenen Bälle in eine CSV Datei.
void saveBallsInCSV(const vector<Ball> &balls, const string &filePath) {
  ofstream ballsfile(filePath);

  if (!ballsfile) {
    cerr << "Fehler beim Öffnen der Datei!" << endl;
    exit(1);
  }

  for (Ball b : balls) {
    ballsfile << b.getCenter().getCoordinates()[0] << ","
              << b.getCenter().getCoordinates()[1] << "," << b.getRadius()
              << "\n";
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
    for (int j = 0; j < c.getPoints().size(); j++) {
      Point p = c.getPoints()[j];
      clusterfile << p.getCoordinates()[0] << "," << p.getCoordinates()[1]
                  << "," << i << endl;
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
    Ball b = findMinEnclosingBall(cluster[i].getPoints());
    if (b.getRadius() != 0) {
      balls.push_back(b);
    }
  }
  return balls;
}

// Berechnet die Summe der Radien einer Liste von Bällen.
double sumOfRadii(vector<Ball> &balls) {
  double radii = 0;
  for (Ball b : balls) {
    radii += b.getRadius();
  }
  return radii;
}

void analyseSchmidt(vector<Point> &points, int k, double epsilon,
                    int numVectors, string &clusterFilePath,
                    string &ballFilePath) {
  auto start = std::chrono::steady_clock::now();
  vector<Cluster> cluster = clustering(points, k, epsilon, numVectors, 10);
  auto end = std::chrono::steady_clock::now();
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  auto diff = end - start;
  cout << "Dauer des Durchlaufs: " << chrono::duration<double>(diff).count()
       << " Sekunden" << endl;
  cout << "Schmidt:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

void analyseGonzales(vector<Point> &points, int k, string &clusterFilePath,
                     string &ballFilePath) {
  vector<Cluster> cluster = gonzales(points, k);
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  cout << "Gonzales:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

void analyseKMeansPlusPlus(vector<Point> &points, int k,
                           string &clusterFilePath, string &ballFilePath) {
  vector<Cluster> cluster = kMeansPlusPlus(points, k);
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  cout << "KMeansPlusPlus:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

void analyseHeuristik(vector<Point> &points, int k, string &clusterFilePath,
                      string &ballFilePath) {
  vector<Cluster> cluster = heuristik(points, k);
  vector<Ball> balls = getBallsFromCluster(cluster);
  double radii = sumOfRadii(balls);
  cout << "Heuristik:                   " << radii << endl;
  saveClusterInCSV(cluster, clusterFilePath);
  saveBallsInCSV(balls, ballFilePath);
}

int main(int argc, char const *argv[]) {
  if (string(argv[1]) == "Schmidt") {
    int k = stod(argv[2]);
    double epsilon = stod(argv[3]);
    int numVectors = stod(argv[4]);
    string pointFilePath = argv[5];
    string ballFilePath = argv[6];
    string clusterFilePath = argv[7];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseSchmidt(points, k, epsilon, numVectors, clusterFilePath,
                   ballFilePath);
  } else if (string(argv[1]) == "Gonzales") {
    int k = stod(argv[2]);
    string pointFilePath = argv[3];
    string ballFilePath = argv[4];
    string clusterFilePath = argv[5];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseGonzales(points, k, clusterFilePath, ballFilePath);
  } else if (string(argv[1]) == "KMeansPlusPlus") {
    int k = stod(argv[2]);
    string pointFilePath = argv[3];
    string ballFilePath = argv[4];
    string clusterFilePath = argv[5];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseKMeansPlusPlus(points, k, clusterFilePath, ballFilePath);
  } else if (string(argv[1]) == "Heuristik") {
    int k = stod(argv[2]);
    string pointFilePath = argv[3];
    string ballFilePath = argv[4];
    string clusterFilePath = argv[5];

    vector<Point> points = readPointsFromCSV(pointFilePath);
    analyseHeuristik(points, k, clusterFilePath, ballFilePath);
  }
  return 0;
}