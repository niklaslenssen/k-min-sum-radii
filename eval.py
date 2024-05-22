import generator
import subprocess
import plot
import matplotlib.pyplot as plt
import pandas as pd
import os
import re


def main():
    # Argumente aus der Konfiguration holen
    config = generator.handle_arguments()

    # Kompilieren des C++ Codes, falls erforderlich
    if config['compile']:
        result = subprocess.run(['g++', '-fopenmp', '-o', 'main', 'main.cpp', 'k_MSR.cpp', 'badoiu_clarkson.cpp', 'welzl.cpp', 'hochbaumShmyos.cpp', 'heuristic.cpp'])
        if result.returncode != 0:
            print('Fehler bei der Kompilierung.')
            return
        print('Erfolgreich kompiliert.')

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
    point_files = [f for f in os.listdir(point_directory) if f.endswith('.csv')]

    # Ausführung der Schmidt und Gonzales Algorithmen
    schmidt(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory)
    gonzales(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory)

    # Vergleich der Ergebnisse der beiden Algorithmen
    compare_schmidt_gonzales(config)


def gonzales(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory):
    # Anzahl der Cluster
    k = config['k']

    # Verzeichnisse erstellen, falls sie nicht existieren
    os.makedirs(os.path.join(ball_directory, 'Gonzales'), exist_ok=True)
    os.makedirs(os.path.join(cluster_directory, 'Gonzales'), exist_ok=True)
    os.makedirs(os.path.join(plot_directory, 'Gonzales'), exist_ok=True)
    os.makedirs(os.path.join(f'Data/{config['dimensions']}/Results', "Gonzales"), exist_ok=True)

    # Regex-Muster zum Extrahieren der Nummer aus dem Dateinamen
    pattern = r"points_(\d+)\.csv"

    results = []

    for point_file in point_files:
        # Nummer aus dem Dateinamen extrahieren
        match = re.search(pattern, point_file)
        if match:
            number = match.group(1)

        point_path = os.path.join(point_directory, point_file)
        ball_path = os.path.join(ball_directory, 'Gonzales', f'balls_{number}.csv')
        cluster_path = os.path.join(cluster_directory, 'Gonzales', f'cluster_{number}.csv')
        plot_path = os.path.join(plot_directory, 'Gonzales', f'plot_{number}.png')

        # Ausführen des C++ Programms mit den entsprechenden Parametern
        result = subprocess.run(['./main', 'g', str(k), point_path, ball_path, cluster_path], capture_output=True, text=True)

        # Ergebnis aus der Programmausgabe parsen
        radii = parse_output(result.stdout)

        # Ergebnisse speichern
        results.append((point_file, radii))

        # Cluster plotten und speichern
        plot.plot_cluster(cluster_path, ball_path, plot_path)
    
    # Ergebnisse sortieren nach Dateiname
    results.sort(key=lambda x: (x[0]))

    # Ergebnisse in eine CSV-Datei schreiben
    with open(f'Data/{config['dimensions']}/Results/Gonzales/results.csv', 'w') as f:
        f.write('Datei,Radii\n')
        for point_file, radii in results:
            f.write(f"{point_file},{radii}\n")


def schmidt(point_files, config, ball_directory, cluster_directory, plot_directory, point_directory):
    # Werte für epsilon und u definieren
    epsilon_values = [0.5]
    u_values = [1] + list(range(100, 201, 100))
    # Anzahl der Cluster
    k = config['k']

    # Verzeichnisse erstellen, falls sie nicht existieren
    os.makedirs(os.path.join(ball_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(cluster_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(plot_directory, "Schmidt"), exist_ok=True)
    os.makedirs(os.path.join(f'Data/{config['dimensions']}/Results', "Schmidt"), exist_ok=True)

    # Regex-Muster zum Extrahieren der Nummer aus dem Dateinamen
    pattern = r"points_(\d+)\.csv"

    results = []
    count = 0
    for point_file in point_files:
        print("Schmidt: " + str(count))
        # Nummer aus dem Dateinamen extrahieren
        match = re.search(pattern, point_file)
        if match:
            number = match.group(1)

        point_path = os.path.join(point_directory, point_file)

        for epsilon in epsilon_values:
            for u in u_values:
                # Definieren der Pfade für die Ball-, Cluster- und Plot-Dateien
                ball_path = os.path.join(ball_directory, 'Schmidt', f'balls_{number}_u{u}_epsilon{epsilon}.csv')
                cluster_path = os.path.join(cluster_directory, 'Schmidt', f'cluster_{number}_u{u}_epsilon{epsilon}.csv')
                plot_path = os.path.join(plot_directory, 'Schmidt', f'plot_{number}_u{u}_epsilon{epsilon}.png')

                # Ausführen des C++ Programms mit den entsprechenden Parametern
                result = subprocess.run(['./main', 's', str(k), str(epsilon), str(u), point_path, ball_path, cluster_path], capture_output=True, text=True)

                # Ergebnisse aus der Programmausgabe parsen
                duration, radii = parse_output_schmidt(result.stdout)

                # Ergebnisse speichern
                results.append((point_file, u, epsilon, duration, radii))

                # Cluster plotten und speichern
                plot.plot_cluster(cluster_path, ball_path, plot_path)
        count += 1

    # Ergebnisse nach Dateiname, 'u' und 'epsilon' sortieren
    results.sort(key=lambda x: (x[0], x[1], -x[2]))

    # Ergebnisse in eine CSV-Datei schreiben
    with open(f'Data/{config['dimensions']}/Results/Schmidt/results.csv', 'w') as f:
        f.write('Datei,u,epsilon,Dauer (Sekunden),Radii\n')
        for point_file, u, epsilon, duration, radii in results:
            f.write(f"{point_file},{u},{epsilon},{duration},{radii}\n")

    # Analysieren und Vergleichen der Ergebnisse
    analyze_results_schmidt(results, config)


def parse_output_schmidt(output):
    duration = None
    radii = None

    # Regex-Muster zur Extraktion der relevanten Informationen
    duration_pattern = re.compile(r"Dauer des Durchlaufs: ([\d\.]+) Sekunden")
    radii_pattern = re.compile(r"Schmidt:\s+([\d\.]+)")

    # Suche nach den Mustern in der Ausgabe
    for line in output.split('\n'):
        duration_match = duration_pattern.search(line)
        radii_match = radii_pattern.search(line)

        if duration_match:
            duration = float(duration_match.group(1))
        if radii_match:
            radii = float(radii_match.group(1))

    return duration, radii


def parse_output(output):
    radii = None

    # Regex-Muster zur Extraktion der Radii-Information
    radii_pattern = re.compile(r"Gonzales:\s+([\d\.]+)")

    # Suche nach dem Muster in der gesamten Ausgabe
    radii_match = radii_pattern.search(output)

    if radii_match:
        radii = float(radii_match.group(1))
    
    return radii


def analyze_results_schmidt(results, config):
    # Ergebnisse in einen DataFrame umwandeln
    df = pd.DataFrame(results, columns=['Datei', 'u', 'epsilon', 'Dauer (Sekunden)', 'Radii'])

    # Boxplot der Radien nach 'u' und 'epsilon'
    plt.figure(figsize=(10, 6))
    df.boxplot(column='Radii', by=['u', 'epsilon'])
    plt.title('Verteilung der Radien nach u und epsilon')
    plt.xlabel('u, epsilon')
    plt.ylabel('Radius')
    plt.xticks(rotation=90)
    plt.savefig(f'Data/{config['dimensions']}/Results/Schmidt/radii_boxplot_all.png')
    plt.close()

    # Boxplots für konstantes u und variierendes epsilon
    for u_val in df['u'].unique():
        df_u = df[df['u'] == u_val]
        plt.figure(figsize=(10, 6))
        df_u.boxplot(column='Radii', by='epsilon')
        plt.title(f'Verteilung der Radien für u={u_val} und variierendes epsilon')
        plt.xlabel('epsilon')
        plt.ylabel('Radius')
        plt.xticks(rotation=90)
        plt.savefig(f'Data/{config['dimensions']}/Results/Schmidt/radii_boxplot_u{u_val}.png')
        plt.close()

    # Boxplots für konstantes epsilon und variierendes u
    for epsilon_val in df['epsilon'].unique():
        df_epsilon = df[df['epsilon'] == epsilon_val]
        plt.figure(figsize=(10, 6))
        df_epsilon.boxplot(column='Radii', by='u')
        plt.title(f'Verteilung der Radien für epsilon={epsilon_val} und variierendes u')
        plt.xlabel('u')
        plt.ylabel('Radius')
        plt.xticks(rotation=90)
        plt.savefig(f'Data/{config['dimensions']}/Results/Schmidt/radii_boxplot_epsilon{epsilon_val}.png')
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
                merged = pd.merge(df_min, df_u_epsilon, on='Datei', suffixes=('_min', f'_u{u}_epsilon{epsilon}'))
                merged['improvement'] = (merged['Radii_min'] - merged[f'Radii_u{u}_epsilon{epsilon}']) / merged['Radii_min'] * 100
                improvements.append((u, epsilon, merged['improvement'].mean()))

    improvement_df = pd.DataFrame(improvements, columns=['u', 'epsilon', 'Improvement (%)'])

    # Verbesserungen anzeigen
    print(f"Verbesserung des Radius im Vergleich zu u={min_u} und epsilon={max_epsilon}:")
    print(improvement_df)

    # Verbesserungen in eine CSV-Datei schreiben
    improvement_df.to_csv(f'Data/{config['dimensions']}/Results/Schmidt/radius_improvement.csv', index=False)

    # Beste Kombination von u und epsilon basierend auf dem kleinsten Radius
    best_combination = df.groupby(['u', 'epsilon'])['Radii'].mean().idxmin()
    best_u, best_epsilon = best_combination
    best_mean_radius = df[(df['u'] == best_u) & (df['epsilon'] == best_epsilon)]['Radii'].mean()
    print(f"Die beste Kombination ist u={best_u} und epsilon={best_epsilon} mit dem kleinsten durchschnittlichen Radius={best_mean_radius:.4f}.")

    # Finde das u, bei dem die Verbesserung der Radien < 1% ist
    stable_u = find_last_significant_improvement_u(df)
    print(f"Das stabile u, bei dem die Verbesserung der Radien < 1% ist: {stable_u}")

    # Markdown-Datei erstellen und Ergebnisse formatieren
    with open(f"Data/{config['dimensions']}/Results/Schmidt/radius_improvement.md", 'w') as md_file:
        md_file.write("# Verbesserung des Radius\n")
        md_file.write(f"Vergleich der Verbesserungen des Radius im Vergleich zu u={min_u} und epsilon={max_epsilon}:\n\n")
        md_file.write("| u | epsilon | Improvement (%) |\n")
        md_file.write("|---|---------|-----------------|\n")
        for row in improvements:
            md_file.write(f"| {row[0]} | {row[1]} | {row[2]:.2f} |\n")
        md_file.write("\n")
        md_file.write(f"Die beste Kombination ist u={best_u} und epsilon={best_epsilon} mit dem kleinsten durchschnittlichen Radius={best_mean_radius:.4f}.\n")
        md_file.write(f"Das stabile u, bei dem die Verbesserung der Radien < 1% ist: {stable_u}")


def find_last_significant_improvement_u(df):
    previous_mean_radius = None
    last_significant_improvement_u = None

    # Schleife über die sortierten Werte von 'u'
    for u in sorted(df['u'].unique()):
        current_mean_radius = df[df['u'] == u]['Radii'].mean()
        if previous_mean_radius is not None:
            improvement = (previous_mean_radius - current_mean_radius) / previous_mean_radius * 100
            if improvement > 1:
                last_significant_improvement_u = u
        previous_mean_radius = current_mean_radius

    return last_significant_improvement_u

def compare_schmidt_gonzales(config):
    # Lade die Ergebnisse
    schmidt_results = pd.read_csv(f'Data/{config["dimensions"]}/Results/Schmidt/results.csv')
    gonzales_results = pd.read_csv(f'Data/{config["dimensions"]}/Results/Gonzales/results.csv')

    # Finden der besten Kombination von u und epsilon für Schmidt
    best_combination = schmidt_results.groupby(['u', 'epsilon'])['Radii'].mean().idxmin()
    best_u, best_epsilon = best_combination
    best_schmidt_results = schmidt_results[(schmidt_results['u'] == best_u) & (schmidt_results['epsilon'] == best_epsilon)]

    # Berechne den durchschnittlichen Radius für jede Methode
    schmidt_avg_radius = best_schmidt_results['Radii'].mean()
    gonzales_avg_radius = gonzales_results['Radii'].mean()

    print(f"Durchschnittlicher Radius für Schmidt (beste Kombination u={best_u}, epsilon={best_epsilon}): {schmidt_avg_radius}")
    print(f"Durchschnittlicher Radius für Gonzales: {gonzales_avg_radius}")

    # Paarweise Vergleiche der Radien
    paired_results = pd.merge(
        best_schmidt_results[['Datei', 'Radii']],
        gonzales_results[['Datei', 'Radii']],
        on='Datei',
        suffixes=('_Schmidt', '_Gonzales')
    )

    paired_results['Better_Algorithm'] = paired_results.apply(
        lambda row: 'Schmidt' if row['Radii_Schmidt'] < row['Radii_Gonzales'] else 'Gonzales',
        axis=1
    )

    # Berechnung der Prozentwerte für den besseren Algorithmus
    schmidt_better_count = paired_results[paired_results['Better_Algorithm'] == 'Schmidt'].shape[0]
    gonzales_better_count = paired_results[paired_results['Better_Algorithm'] == 'Gonzales'].shape[0]
    total_count = paired_results.shape[0]

    schmidt_better_percentage = (schmidt_better_count / total_count) * 100
    gonzales_better_percentage = (gonzales_better_count / total_count) * 100

    print(f"Schmidt liefert bessere Ergebnisse in {schmidt_better_percentage:.2f}% der Fälle.")
    print(f"Gonzales liefert bessere Ergebnisse in {gonzales_better_percentage:.2f}% der Fälle.")

    # Ergebnisse in ein DataFrame packen für den Vergleich
    comparison_df = pd.DataFrame({
        'Methode': ['Schmidt', 'Gonzales'],
        'Durchschnittlicher Radius': [schmidt_avg_radius, gonzales_avg_radius],
        'Besser in % der Fälle': [schmidt_better_percentage, gonzales_better_percentage]
    })

    # Ergebnisse als Tabelle speichern
    comparison_df.to_csv(f'Data/{config["dimensions"]}/Results/comparison.csv', index=False)



if __name__ == "__main__":
    main()