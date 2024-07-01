#include "header/heuristic.h"
#include "header/k_MSR.h"
#include "header/point.h"

using namespace std;

struct PointData {
  double* coordinates;
  int dimension;
};

struct ClusterData {
  PointData* points;
  int num_points;
};

vector<Point> array_to_vector(double* array, int num_points, int dimension) {
  vector<Point> points;
  for (int i = 0; i < num_points; i++) {
    vector<double> coordinates;
    for (int j = 0; j < dimension; j++) {
      coordinates.push_back(array[i * dimension + j]);
    }
    points.push_back(Point(coordinates));
  }
  return points;
}

ClusterData* cluster_to_array(vector<Cluster> clusters, int* num_clusters) {
  *num_clusters = clusters.size();

  ClusterData* cluster_data = new ClusterData[*num_clusters];
  for (int i = 0; i < *num_clusters; i++) {
    const vector<Point>& cluster_points = clusters[i].getPoints();
    cluster_data[i].num_points = cluster_points.size();
    cluster_data[i].points = new PointData[cluster_points.size()];
    for (int j = 0; j < cluster_points.size(); j++) {
      const vector<double>& coords = cluster_points[j].getCoordinates();
      cluster_data[i].points[j].dimension = coords.size();
      cluster_data[i].points[j].coordinates = new double[coords.size()];
      for (int k = 0; k < coords.size(); k++) {
        cluster_data[i].points[j].coordinates[k] = coords[k];
      }
    }
  }
  return cluster_data;
}

extern "C" {

ClusterData* schmidt_wrapper(double* point_array, int num_points,
                                int dimension, int k, double epsilon,
                                int numUVectors, int numRadiiVectors, int* num_clusters) {
  vector<Point> points = array_to_vector(point_array, num_points, dimension);

  vector<Cluster> cluster = clustering(points, k, epsilon, numUVectors, numRadiiVectors);

  return cluster_to_array(cluster, num_clusters);
}

ClusterData* heuristic_wrapper(double* point_array, int num_points,
                               int dimension, int k, int* num_clusters) {
  vector<Point> points = array_to_vector(point_array, num_points, dimension);

  vector<Cluster> cluster = heuristik(points, k);

  return cluster_to_array(cluster, num_clusters);
}

ClusterData* gonzales_wrapper(double* point_array, int num_points,
                              int dimension, int k, int* num_clusters) {
  vector<Point> points = array_to_vector(point_array, num_points, dimension);

  vector<Cluster> cluster = gonzales(points, k);

  return cluster_to_array(cluster, num_clusters);
}

ClusterData* kmeans_wrapper(double* point_array, int num_points, int dimension,
                            int k, int* num_clusters) {
  vector<Point> points = array_to_vector(point_array, num_points, dimension);

  vector<Cluster> cluster = kMeansPlusPlus(points, k);

  return cluster_to_array(cluster, num_clusters);
}

void free_cluster_data(ClusterData* cluster_data, int num_clusters) {
  for (int i = 0; i < num_clusters; i++) {
    for (int j = 0; j < cluster_data[i].num_points; j++) {
      delete[] cluster_data[i].points[j].coordinates;
    }
    delete[] cluster_data[i].points;
  }
  delete[] cluster_data;
}
}
