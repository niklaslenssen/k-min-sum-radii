import generator
import plot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re
import csv
import ctypes
import time
import miniball

lib = ctypes.CDLL('./clustering.so')


class PointData(ctypes.Structure):
    _fields_ = [('coordinates', ctypes.POINTER(ctypes.c_double)),
                ('dimension', ctypes.c_int)]


class ClusterData(ctypes.Structure):
    _fields_ = [('points', ctypes.POINTER(PointData)),
                ('num_points', ctypes.c_int)]


class BallData(ctypes.Structure):
    _fields_ = [('center', ctypes.POINTER(ctypes.c_double)),
                ('radius', ctypes.c_double)]


lib.clustering_wrapper.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
lib.clustering_wrapper.restype = ctypes.POINTER(ClusterData)

lib.heuristic_wrapper.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
lib.heuristic_wrapper.restype = ctypes.POINTER(ClusterData)

lib.gonzales_wrapper.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
lib.gonzales_wrapper.restype = ctypes.POINTER(ClusterData)

lib.kmeans_wrapper.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int)
]
lib.kmeans_wrapper.restype = ctypes.POINTER(ClusterData)

lib.free_cluster_data.argtypes = [ctypes.POINTER(ClusterData), ctypes.c_int]
lib.free_cluster_data.restype = None


def calculate_miniball(points):
    mb = miniball.Miniball(points)
    radius = np.sqrt(mb.squared_radius())
    return mb.center(), radius


def save_clusters_to_csv(clusters, output_file):
    data_to_save = []
    for cluster_idx, cluster in enumerate(clusters):
        for point in cluster:
            data_to_save.append(list(point) + [cluster_idx])

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['x', 'y', 'cluster'])
        writer.writerows(data_to_save)


def extract_clusters(cluster_ptr, num_clusters_value):
    clusters = cluster_ptr[:num_clusters_value]
    cluster_points_list = []

    for cluster_idx in range(num_clusters_value):
        cluster = clusters[cluster_idx]
        cluster_points = []
        for point_idx in range(cluster.num_points):
            point = cluster.points[point_idx]
            coords = [point.coordinates[dim] for dim in range(point.dimension)]
            cluster_points.append(coords)
        cluster_points_list.append(np.array(cluster_points))

    return np.array(cluster_points_list, dtype=object)


def save_mebs_to_csv(clusters, output_file):
    meb_data_to_save = []
    sum_of_radii = 0
    for _, cluster in enumerate(clusters):
        if len(cluster) > 0:
            center, radius = calculate_miniball(np.array(cluster))
            sum_of_radii += radius
            meb_data_to_save.append([center[0], center[1], radius])

    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['center_x', 'center_y', 'radius'])
        writer.writerows(meb_data_to_save)

    return sum_of_radii


def read_points_from_csv(input_file):
    points_array = np.loadtxt(input_file, delimiter=',', usecols=[1, 2])
    points_array = np.ascontiguousarray(points_array)

    c_array = points_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    numPoints = len(points_array)

    return c_array, numPoints


def call_clustering_function(clustering_func, c_array, numPoints, dimensions, k, num_clusters):
    cluster_ptr = clustering_func(
        c_array,
        numPoints,
        dimensions,
        k,
        num_clusters
    )
    return cluster_ptr


def cluster(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory, algorithm, c_function):
    # Anzahl der Cluster
    k = config['k']

    # Verzeichnisse erstellen, falls sie nicht existieren
    os.makedirs(os.path.join(ball_directory, algorithm), exist_ok=True)
    os.makedirs(os.path.join(cluster_directory, algorithm), exist_ok=True)
    os.makedirs(os.path.join(plot_directory, algorithm), exist_ok=True)
    os.makedirs(os.path.join(
        f'Data/{config['dimensions']}/Results', algorithm), exist_ok=True)

    # Regex-Muster zum Extrahieren der Nummer aus dem Dateinamen
    pattern = r"points_(\d+)\.csv"
    # Regex-Muster zur Extraktion der Radii-Information
    radii_pattern = re.compile(rf"{algorithm}:\s+([\d\.]+)")

    results = []

    for point_file in point_files:
        # Nummer aus dem Dateinamen extrahieren
        match = re.search(pattern, point_file)
        if match:
            number = match.group(1)

        point_path = os.path.join(point_directory, point_file)
        ball_path = os.path.join(
            ball_directory, algorithm, f'balls_{number}.csv')
        cluster_path = os.path.join(
            cluster_directory, algorithm, f'cluster_{number}.csv')
        plot_path = os.path.join(
            plot_directory, algorithm, f'plot_{number}.png')

        c_array, numPoints = read_points_from_csv(point_path)

        num_clusters = ctypes.c_int()

        start_time = time.time()

        cluster_ptr = call_clustering_function(
            c_function,
            c_array,
            numPoints,
            config['dimensions'],
            config['k'],
            ctypes.byref(num_clusters)
        )

        end_time = time.time()

        duration = end_time - start_time

        num_clusters_value = num_clusters.value

        cluster = extract_clusters(cluster_ptr, num_clusters_value)

        lib.free_cluster_data(cluster_ptr, num_clusters)

        save_clusters_to_csv(cluster, cluster_path)

        radii = save_mebs_to_csv(cluster, ball_path)

        # Ergebnisse speichern
        results.append((point_file, duration,  radii))

        # Cluster plotten und speichern
        plot.plot_cluster(cluster_path, ball_path, plot_path)

    # Ergebnisse sortieren nach Dateiname
    results.sort(key=lambda x: (x[0]))

    # Ergebnisse in eine CSV-Datei schreiben
    with open(f'Data/{config['dimensions']}/Results/{algorithm}/results.csv', 'w') as f:
        f.write('Datei,Dauer (Sekunden),Radii\n')
        for point_file, duration, radii in results:
            f.write(f"{point_file},{duration},{radii}\n")


def schmidt(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory):
    # Werte für epsilon und u definieren
    epsilon_values = [0.5]
    u_values = [1, 10, 100, 1000, 2000, 3000]
    # Anzahl der Cluster
    k = config['k']
    dimension = config['dimensions']

    # Verzeichnisse erstellen, falls sie nicht existieren
    os.makedirs(os.path.join(ball_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(cluster_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(plot_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(
        f'Data/{dimension}/Results', "Schmidt"), exist_ok=True)

    # Regex-Muster zum Extrahieren der Nummer aus dem Dateinamen
    pattern = r"points_(\d+)\.csv"

    results = []
    count = 1
    for point_file in point_files:
        # Nummer aus dem Dateinamen extrahieren
        match = re.search(pattern, point_file)
        if match:
            number = match.group(1)

        point_path = os.path.join(point_directory, point_file)

        for epsilon in epsilon_values:
            for u in u_values:
                # Definieren der Pfade für die Ball-, Cluster- und Plot-Dateien
                ball_path = os.path.join(ball_directory, 'Schmidt', f'balls_{
                                         number}_u{u}_epsilon{epsilon}.csv')
                cluster_path = os.path.join(cluster_directory, 'Schmidt', f'cluster_{
                                            number}_u{u}_epsilon{epsilon}.csv')
                plot_path = os.path.join(plot_directory, 'Schmidt', f'plot_{
                                         number}_u{u}_epsilon{epsilon}.png')

                c_array, numPoints = read_points_from_csv(point_path)

                num_clusters = ctypes.c_int()

                start_time = time.time()

                # Rufe die C-Funktion auf
                cluster_ptr = lib.clustering_wrapper(
                    c_array,
                    numPoints,
                    dimension,
                    k,
                    epsilon,
                    u,
                    ctypes.byref(num_clusters)
                )

                end_time = time.time()

                duration = end_time - start_time

                num_clusters_value = num_clusters.value

                cluster = extract_clusters(cluster_ptr, num_clusters_value)

                lib.free_cluster_data(cluster_ptr, num_clusters)

                save_clusters_to_csv(cluster, cluster_path)

                radii = save_mebs_to_csv(cluster, ball_path)

                # Ergebnisse speichern
                results.append((point_file, u, epsilon, duration, radii))

                # Cluster plotten und speichern
                plot.plot_cluster(cluster_path, ball_path, plot_path)

        print(count)
        count += 1

    # Ergebnisse nach Dateiname, 'u' und 'epsilon' sortieren
    results.sort(key=lambda x: (x[0], x[1], -x[2]))

    # Ergebnisse in eine CSV-Datei schreiben
    with open(f'Data/{dimension}/Results/Schmidt/results.csv', 'w') as f:
        f.write('Datei,u,epsilon,Dauer (Sekunden),Radii\n')
        for point_file, u, epsilon, duration, radii in results:
            f.write(f"{point_file},{u},{epsilon},{duration},{radii}\n")


def analyze_results_schmidt(config):
    # Ergebnisse in einen DataFrame umwandeln
    df = pd.read_csv(
        f'Data/{config['dimensions']}/Results/Schmidt/results.csv')

    # Boxplot der Radien nach 'u' und 'epsilon'
    plt.figure(figsize=(10, 6))
    df.boxplot(column='Radii', by=['u', 'epsilon'])
    plt.title('Verteilung der Radien nach u und epsilon')
    plt.suptitle('')
    plt.xlabel('u, epsilon')
    plt.ylabel('Summe der Radien')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'Data/{config['dimensions']
                        }/Results/Schmidt/radii_boxplot_all.png')
    plt.close()

    # Boxplots für konstantes u und variierendes epsilon
    for u_val in df['u'].unique():
        df_u = df[df['u'] == u_val]
        plt.figure(figsize=(10, 6))
        df_u.boxplot(column='Radii', by='epsilon')
        plt.title(f'Verteilung der Radien für u={
                  u_val} und variierendes epsilon')
        plt.suptitle('')
        plt.xlabel('epsilon')
        plt.ylabel('Summe der Radien')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(
            f'Data/{config['dimensions']}/Results/Schmidt/radii_boxplot_u{u_val}.png')
        plt.close()

    # Boxplots für konstantes epsilon und variierendes u
    for epsilon_val in df['epsilon'].unique():
        df_epsilon = df[df['epsilon'] == epsilon_val]
        plt.figure(figsize=(10, 6))
        df_epsilon.boxplot(column='Radii', by='u')
        plt.title(f'Verteilung der Radien für epsilon={
                  epsilon_val} und variierendes u')
        plt.suptitle('')
        plt.xlabel('u')
        plt.ylabel('Summe der Radien')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(
            f'Data/{config['dimensions']}/Results/Schmidt/radii_boxplot_epsilon{epsilon_val}.png')
        plt.close()

    # Boxplot der Dauer der Durchläufe nach 'u'
    plt.figure(figsize=(10, 6))
    df.boxplot(column='Dauer (Sekunden)', by='u')
    plt.title('Verteilung der Dauer nach u')
    plt.suptitle('')
    plt.xlabel('u')
    plt.ylabel('Dauer (Sekunden)')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f'Data/{config["dimensions"]
                        }/Results/Schmidt/duration_boxplot.png')
    plt.close()

    # Berechnung der Verbesserung des Radius
    min_u = min(df['u'])
    max_epsilon = max(df['epsilon'])
    df_min = df[(df['u'] == min_u) & (df['epsilon'] == max_epsilon)]

    improvements = []
    for u in df['u'].unique():
        for epsilon in df['epsilon'].unique():
            if u != min_u or epsilon != max_epsilon:
                df_u_epsilon = df[(df['u'] == u) & (df['epsilon'] == epsilon)]
                merged = pd.merge(df_min, df_u_epsilon, on='Datei', suffixes=(
                    '_min', f'_u{u}_epsilon{epsilon}'))
                merged['improvement'] = (
                    merged['Radii_min'] - merged[f'Radii_u{u}_epsilon{epsilon}']) / merged['Radii_min'] * 100
                improvements.append((u, epsilon, merged['improvement'].mean()))

    improvement_df = pd.DataFrame(
        improvements, columns=['u', 'epsilon', 'Improvement (%)'])

    # Verbesserungen anzeigen
    print(f"Verbesserung des Radius im Vergleich zu u={
          min_u} und epsilon={max_epsilon}:")
    print(improvement_df)

    # Verbesserungen in eine CSV-Datei schreiben
    improvement_df.to_csv(
        f'Data/{config['dimensions']}/Results/Schmidt/radius_improvement.csv', index=False)

    # Beste Kombination von u und epsilon basierend auf dem kleinsten Radius
    best_combination = df.groupby(['u', 'epsilon'])['Radii'].mean().idxmin()
    best_u, best_epsilon = best_combination
    best_mean_radius = df[(df['u'] == best_u) & (
        df['epsilon'] == best_epsilon)]['Radii'].mean()
    print(f"Die beste Kombination ist u={best_u} und epsilon={
          best_epsilon} mit dem kleinsten durchschnittlichen Radius={best_mean_radius:.4f}.")

    # Markdown-Datei erstellen und Ergebnisse formatieren
    with open(f"Data/{config['dimensions']}/Results/Schmidt/radius_improvement.md", 'w') as md_file:
        md_file.write("# Verbesserung des Radius\n")
        md_file.write(f"Vergleich der Verbesserungen des Radius im Vergleich zu u={
                      min_u} und epsilon={max_epsilon}:\n\n")
        md_file.write("| u | epsilon | Improvement (%) |\n")
        md_file.write("|---|---------|-----------------|\n")
        for row in improvements:
            md_file.write(f"| {row[0]} | {row[1]} | {row[2]:.2f} |\n")
        md_file.write("\n")
        md_file.write(f"Die beste Kombination ist u={best_u} und epsilon={
                      best_epsilon} mit dem kleinsten durchschnittlichen Radius={best_mean_radius:.4f}.\n")


def compare_algorithms(config):
    # Lade die Ergebnisse
    schmidt_results = pd.read_csv(
        f'Data/{config["dimensions"]}/Results/Schmidt/results.csv')
    gonzales_results = pd.read_csv(
        f'Data/{config["dimensions"]}/Results/Gonzales/results.csv')
    kmeans_results = pd.read_csv(
        f'Data/{config["dimensions"]}/Results/KMeansPlusPlus/results.csv')
    heuristik_results = pd.read_csv(
        f'Data/{config["dimensions"]}/Results/Heuristik/results.csv')

    all_comparison_results = []

    for u_val in schmidt_results['u'].unique():
        # Finden der besten Kombination von epsilon für das aktuelle u
        best_combination = schmidt_results[schmidt_results['u'] == u_val].groupby(
            ['u', 'epsilon'])['Radii'].mean().idxmin()
        best_u, best_epsilon = best_combination
        best_schmidt_results = schmidt_results[(schmidt_results['u'] == best_u) & (
            schmidt_results['epsilon'] == best_epsilon)]

        # Paarweise Vergleiche der Radien
        paired_results = pd.merge(best_schmidt_results[['Datei', 'Radii']], gonzales_results[[
                                  'Datei', 'Radii']], on='Datei', suffixes=('_Schmidt', '_Gonzales'))
        paired_results = pd.merge(
            paired_results, kmeans_results[['Datei', 'Radii']], on='Datei')
        paired_results.rename(columns={'Radii': 'Radii_KMeans'}, inplace=True)
        paired_results = pd.merge(
            paired_results, heuristik_results[['Datei', 'Radii']], on='Datei')
        paired_results.rename(
            columns={'Radii': 'Radii_Heuristik'}, inplace=True)

        # Zählen, wie oft Schmidt besser als alle anderen ist
        paired_results['Schmidt_better_than_all'] = paired_results.apply(
            lambda row: row['Radii_Schmidt'] < row[[
                'Radii_Gonzales', 'Radii_KMeans', 'Radii_Heuristik']].min(),
            axis=1
        )

        # Zählen, wie oft Schmidt schlechter als alle anderen ist
        paired_results['Schmidt_worse_than_all'] = paired_results.apply(
            lambda row: row['Radii_Schmidt'] > row[[
                'Radii_Gonzales', 'Radii_KMeans', 'Radii_Heuristik']].max(),
            axis=1
        )

        schmidt_better_count = paired_results['Schmidt_better_than_all'].sum()
        schmidt_worse_count = paired_results['Schmidt_worse_than_all'].sum()

        # Paarweise Vergleiche
        schmidt_vs_gonzales_better = (
            paired_results['Radii_Schmidt'] < paired_results['Radii_Gonzales']).sum()
        schmidt_vs_gonzales_worse = (
            paired_results['Radii_Schmidt'] > paired_results['Radii_Gonzales']).sum()

        schmidt_vs_kmeans_better = (
            paired_results['Radii_Schmidt'] < paired_results['Radii_KMeans']).sum()
        schmidt_vs_kmeans_worse = (
            paired_results['Radii_Schmidt'] > paired_results['Radii_KMeans']).sum()

        schmidt_vs_heuristik_better = (
            paired_results['Radii_Schmidt'] < paired_results['Radii_Heuristik']).sum()
        schmidt_vs_heuristik_worse = (
            paired_results['Radii_Schmidt'] > paired_results['Radii_Heuristik']).sum()

        total_count = paired_results.shape[0]

        schmidt_better_percentage = (schmidt_better_count / total_count) * 100
        schmidt_worse_percentage = (schmidt_worse_count / total_count) * 100
        schmidt_vs_gonzales_better_percentage = (
            schmidt_vs_gonzales_better / total_count) * 100
        schmidt_vs_gonzales_worse_percentage = (
            schmidt_vs_gonzales_worse / total_count) * 100

        schmidt_vs_kmeans_better_percentage = (
            schmidt_vs_kmeans_better / total_count) * 100
        schmidt_vs_kmeans_worse_percentage = (
            schmidt_vs_kmeans_worse / total_count) * 100

        schmidt_vs_heuristik_better_percentage = (
            schmidt_vs_heuristik_better / total_count) * 100
        schmidt_vs_heuristik_worse_percentage = (
            schmidt_vs_heuristik_worse / total_count) * 100

        # Ergebnisse sammeln
        all_comparison_results.append({
            'u': u_val,
            'Schmidt vs Alle Besser (%)': schmidt_better_percentage,
            'Schmidt vs Alle Schlechter (%)': schmidt_worse_percentage,
            'Schmidt vs Gonzales Besser (%)': schmidt_vs_gonzales_better_percentage,
            'Schmidt vs Gonzales Schlechter (%)': schmidt_vs_gonzales_worse_percentage,
            'Schmidt vs KMeans++ Besser (%)': schmidt_vs_kmeans_better_percentage,
            'Schmidt vs KMeans++ Schlechter (%)': schmidt_vs_kmeans_worse_percentage,
            'Schmidt vs Heuristik Besser (%)': schmidt_vs_heuristik_better_percentage,
            'Schmidt vs Heuristik Schlechter (%)': schmidt_vs_heuristik_worse_percentage
        })

    # Ergebnisse in ein DataFrame packen für den Vergleich
    comparison_df = pd.DataFrame(all_comparison_results)

    # Ergebnisse als Tabelle speichern
    comparison_df.to_csv(
        f'Data/{config["dimensions"]}/Results/comparison_all_us.csv', index=False)

    # Markdown-Datei erstellen
    with open(f'Data/{config["dimensions"]}/Results/comparison_all_us.md', 'w') as md_file:
        md_file.write(
            "# Paarweiser Vergleich der Clustering-Algorithmen für alle u-Werte\n\n")
        md_file.write("| u  | Schmidt vs Alle Besser (%) | Schmidt vs Alle Schlechter (%) | Schmidt vs Gonzales Besser (%) | Schmidt vs Gonzales Schlechter (%) | Schmidt vs KMeans++ Besser (%) | Schmidt vs KMeans++ Schlechter (%) | Schmidt vs Heuristik Besser (%) | Schmidt vs Heuristik Schlechter (%) |\n")
        md_file.write("|----|----------------------------|--------------------------------|--------------------------------|------------------------------------|--------------------------------|------------------------------------|---------------------------------|-------------------------------------|\n")
        for _, row in comparison_df.iterrows():
            md_file.write(f"| {int(row['u'])} | {row['Schmidt vs Alle Besser (%)']:.2f} | {row['Schmidt vs Alle Schlechter (%)']:.2f} | {row['Schmidt vs Gonzales Besser (%)']:.2f} | {row['Schmidt vs Gonzales Schlechter (%)']:.2f} | {
                          row['Schmidt vs KMeans++ Besser (%)']:.2f} | {row['Schmidt vs KMeans++ Schlechter (%)']:.2f} | {row['Schmidt vs Heuristik Besser (%)']:.2f} | {row['Schmidt vs Heuristik Schlechter (%)']:.2f} |\n")


def main():
    # Argumente aus der Konfiguration holen
    config = generator.handle_arguments()

    # Verzeichnisse definieren
    directory = f'Data/{config["dimensions"]}'
    point_directory = os.path.join(directory, 'Points')
    ball_directory = os.path.join(directory, 'Balls')
    cluster_directory = os.path.join(directory, 'Cluster')
    plot_directory = os.path.join(directory, 'Plots')

    # Überprüfen, ob die Punkte-Dateien existieren, andernfalls generieren
    if not os.path.exists(os.path.join(point_directory, f'points_{config["number_files"] - 1}.csv')):
        generator.generate_data(config)

    # Liste der Punkte-Dateien im Verzeichnis
    point_files = [f for f in os.listdir(
        point_directory) if f.endswith('.csv')]

    # Ausführung der Algorithmen
    schmidt(point_files, config, ball_directory,
            cluster_directory, plot_directory, point_directory)
    cluster(point_files, config, ball_directory, cluster_directory,
            plot_directory, point_directory, 'Heuristik', lib.heuristic_wrapper)
    cluster(point_files, config, ball_directory, cluster_directory,
            plot_directory, point_directory, 'Gonzales', lib.gonzales_wrapper)
    cluster(point_files, config, ball_directory, cluster_directory,
            plot_directory, point_directory, 'KMeansPlusPlus', lib.kmeans_wrapper)

    # Analyse und Vergleich der Ergebnisse von Schmidt
    analyze_results_schmidt(config)

    # Vergleich der Ergebnisse der verschiedenen Algorithmen
    compare_algorithms(config)


if __name__ == "__main__":
    config = generator.handle_arguments()
    main()
