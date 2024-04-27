"""
This file is used to generate random clusters by generating points around random centers.
Always creates new files and does not override old ones.
Run with -h for help.
"""

import argparse
import math
import matplotlib.pyplot as plt
import random
import os

def generate_centers(config):
    centers = []
    for _ in range(config["number_centers"]):
        x = random.uniform(0, config["max_x"])
        y = random.uniform(0, config["max_y"])
        centers.append((x, y))

    return centers


def generate_points_around_center(center, config):
    # generate config for this cluster
    number_points = random.randint(
        config["min_points_per_cluster"], config["max_points_per_cluster"]
    )
    cluster_radius = random.uniform(
        config["min_cluster_radius"], config["max_cluster_radius"]
    )

    # generate points for this cluster
    points = []
    count = 0
    while count < number_points:
        radius = random.uniform(0, cluster_radius)
        angle = random.uniform(0, 2 * math.pi)
        x = radius * math.cos(angle) + center[0]
        y = radius * math.sin(angle) + center[1]
        # (x,y) out of bounds? -> try again
        if x < 0 or x > config["max_x"] or y < 0 or y > config["max_y"]:
            continue
        points.append((x, y))
        count += 1

    return points


def generate_clusters(config):
    centers = generate_centers(config)

    points = []
    for center in centers:
        points = points + generate_points_around_center(center, config)

    return centers, points


def write_clusters_to_file(points, config):
    # save config in file name
    file_name = (
        "RandomCenter"
        + f"-{config['max_x']}"
        + f"-{config['max_y']}"
        + f"-{config['number_centers']}"
        + f"-{config['min_points_per_cluster']}"
        + f"-{config['max_points_per_cluster']}"
        + f"-{config['min_cluster_radius']}"
        + f"-{config['max_cluster_radius']}"
        + "-nr-"
    )

    nr = 0
    generated_file = False
    while not generated_file:
        # make sure to create a new file
        if os.path.isfile(file_name + str(nr) + ".csv"):
            nr += 1
            continue

        # write points to file
        file = open("Data/points.csv", "w")
        for point in points:
            file.write(f"2,{point[0]},{point[1]}\n")
        file.close()
        generated_file = True


def handle_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-p",
        "--only-print",
        action="store_true",
        help="print data and do not create files",
        default=False,
    )
    parser.add_argument(
        "-nf",
        "--number-files",
        type=int,
        help="number of files that will be generated with the given config",
        default=1,
    )
    parser.add_argument(
        "-x",
        "--max-x",
        type=int,
        help="generated points will have an x-value from 0 to max-x",
        default=1,
    )
    parser.add_argument(
        "-y",
        "--max-y",
        type=int,
        help="generated points will have a y-value from 0 to max-y",
        default=1,
    )
    parser.add_argument(
        "-n",
        "--number-centers",
        type=int,
        help="number of centers/clusters that will be generated",
        default=3,
    )
    parser.add_argument(
        "-minp",
        "--min-points-per-cluster",
        type=int,
        help="minimum number of points per cluster (center excluded)",
        default=30,
    )
    parser.add_argument(
        "-maxp",
        "--max-points-per-cluster",
        type=int,
        help="maximum number of points per cluster (center excluded)",
        default=50,
    )
    parser.add_argument(
        "-minr",
        "--min-cluster-radius",
        type=float,
        help="minimum radius of generated clusters",
        default=0.05,
    )
    parser.add_argument(
        "-maxr",
        "--max-cluster-radius",
        type=float,
        help="minimum radius of generated clusters",
        default=0.1,
    )
    args = parser.parse_args()

    # sanity checks
    assert args.max_x > 0
    assert args.max_y > 0
    assert args.number_centers > 0
    assert args.min_points_per_cluster > 0
    assert args.min_points_per_cluster <= args.max_points_per_cluster
    assert args.min_cluster_radius > 0
    assert args.min_cluster_radius <= args.max_cluster_radius

    return vars(args)


def print_clusters(centers, points, config):
    plt.xlim(0, config["max_x"])
    plt.ylim(0, config["max_y"])
    plt.plot([p[0] for p in centers], [p[1] for p in centers], "r.")
    plt.plot([p[0] for p in points], [p[1] for p in points], "b.")
    plt.savefig('points.png')


def main():
    config = handle_arguments()
    print(f"Generating cluster data with the following config:\n{config}")
    for _ in range(config["number_files"]):
        centers, points = generate_clusters(config)
        if config["only_print"]:
            print_clusters(centers, points, config)
        else:
            write_clusters_to_file(centers + points, config)


if __name__ == "__main__":
    main()
