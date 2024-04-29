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

    



    filename = "Data/points.csv"

    imageOutputSimpleGonzalez = f"Data/Plots/PlainGonzalez.png"
    imageOutputWithoutMerge = f"Data/Plots/WithoutMerge.png"
    imageOutputWithMerge = f"Data/Plots/WithMerge.png"
    imageOutputWithMergeIterative = f"Data/Plots/WithMergeIterative.png"

    #subprocess.run(['g++', '-fopenmp', '-o', 'main', 'main.cpp', 'k_MSR.cpp', 'badoiu_clarkson.cpp', 'welzl.cpp', 'hochbaumShmyos.cpp' ])
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

    # _, clustered_points, rad_sum = Algorithms.merging_k_msr(points, k, use_iterative_gonzalez=False)
    # radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    # Points.save_plot(imageOutputWithMerge, clustered_points, showRadii=True, radii=radii)
    # print(f"Merging Gonzales:          {round(rad_sum, 6)}")

    # _, clustered_points, rad_sum = Algorithms.merging_k_msr(points, k)
    # radii = Algorithms.get_radii_of_clustering(clustered_points, k)
    # Points.save_plot(imageOutputWithMergeIterative, clustered_points, showRadii=True, radii=radii)
    # print(f"Merging Gonzales iterativ: {round(rad_sum, 6)}")

    fig, axs = plt.subplots(2, 2)

    # Plot 1 zu Subplot (0, 0) hinzufügen
    axs[0, 0].imshow(plt.imread(imageOutputSimpleGonzalez))
    axs[0, 0].set_title("Gonzales")
    axs[0, 0].axis('off')

    # Plot 2 zu Subplot (0, 1) hinzufügen
    axs[0, 1].imshow(plt.imread(imageOutputSimpleGonzalez))
    axs[0, 1].set_title("Gonzales")
    axs[0, 1].axis('off')

    # Plot 3 zu Subplot (1, 0) hinzufügen
    axs[1, 0].imshow(plt.imread(imageOutputWithoutMerge))
    axs[1, 0].set_title("Gonzales iterativ")
    axs[1, 0].axis('off')

    # Plot 4 zu Subplot (1, 1) hinzufügen
    axs[1, 1].imshow(plt.imread(imageOutputWithoutMerge))
    axs[1, 1].set_title("Gonzales")
    axs[1, 1].axis('off')

    # Bild speichern
    plt.savefig("combined_plots.png")

    return


if __name__ == "__main__":
    main()