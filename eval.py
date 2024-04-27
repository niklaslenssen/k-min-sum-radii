import SHK.Algorithms as Algorithms
import SHK.RandomCenters as RandomCenters
import subprocess
import plot

def main():
    k = 3
    epsilon = 0.5

    config = RandomCenters.handle_arguments()
    print(f"Generating cluster data with the following config:\n{config}")
    for _ in range(config["number_files"]):
        centers, points = RandomCenters.generate_clusters(config)
        if config["only_print"]:
            RandomCenters.print_clusters(centers, points, config)
        else:
            RandomCenters.write_clusters_to_file(centers + points, config)



    filename = "Data/points.csv"

    # subprocess.run(['g++', '-fopenmp', '-o', 'main', 'main.cpp', 'k_MSR.cpp', 'badoiu_clarkson.cpp', 'welzl.cpp', 'hochbaumShmyos.cpp' ])
    subprocess.run(['./main', f'{k}', f'{epsilon}'])


    plot.plot_cluster()



    return








if __name__ == "__main__":
    main()