import SHK.Algorithms as Algorithms
import SHK.RandomCenters as RandomCenters
import SHK.Points as Points
import subprocess
import plot
import matplotlib.pyplot as plt

def main():
    config = RandomCenters.handle_arguments()
    k = config['number_centers']
    epsilon = 0.5

    
    for _ in range(config["number_files"]):
        centers, points = RandomCenters.generate_clusters(config)
        if config["only_print"]:
            RandomCenters.print_clusters(centers, points, config)
        else:
            RandomCenters.write_clusters_to_file(centers + points, config)



    filename = "Data/points.csv"

    imageOutputSimpleGonzalez = f"Data/Plots/PlainGonzalez.png"
    imageOutputWithoutMerge = f"Data/Plots/WithoutMerge.png"
    imageOutputWithMerge = f"Data/Plots/WithMerge.png"
    imageOutputWithMergeIterative = f"Data/Plots/WithMergeIterative.png"

    subprocess.run(['g++', '-fopenmp', '-o', 'main', 'main.cpp', 'k_MSR.cpp', 'badoiu_clarkson.cpp', 'welzl.cpp', 'hochbaumShmyos.cpp' ])
    subprocess.run(['./main', f'{k}', f'{epsilon}'])
    plot.plot_cluster()

    points = Points.read_points(filename)

    _, clustered_points, rad_sum = Algorithms.gonzalez_k_msr(points, k)
    radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    Points.save_plot(imageOutputSimpleGonzalez, clustered_points, showRadii=True, radii=radii)
    print(f"Gonzales:                  {round(rad_sum, 6)}")

    _, clustered_points, rad_sum = Algorithms.iterative_gonzalez_k_msr(points, k)
    radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    Points.save_plot(imageOutputWithoutMerge, clustered_points, showRadii=True, radii=radii)
    print(f"Gonzales iterativ:         {round(rad_sum, 6)}")

    _, clustered_points, rad_sum = Algorithms.merging_k_msr(points, k, use_iterative_gonzalez=False)
    radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    Points.save_plot(imageOutputWithMerge, clustered_points, showRadii=True, radii=radii)
    print(f"Merging Gonzales:          {round(rad_sum, 6)}")

    _, clustered_points, rad_sum = Algorithms.merging_k_msr(points, k)
    radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    Points.save_plot(imageOutputWithMergeIterative, clustered_points, showRadii=True, radii=radii)
    print(f"Merging Gonzales iterativ: {round(rad_sum, 6)}")

    return


if __name__ == "__main__":
    main()