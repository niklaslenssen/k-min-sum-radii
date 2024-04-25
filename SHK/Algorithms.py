# -*- coding: utf-8 -*-
"""
Algorithem-Collection
"""

# Hier sollen später die Algorithmen Implementiert werden
import numpy as np
import copy


# # # =========== Gonzalez =========== # # #
# Einfacher Gonzalez mit Zufallsindex oder festgewähltem Startindex
def gonzalez(points_in, k, firstIndex = -1):
    # Deep-Copy, damit die Ursprungspunkte nicht verändert werden
    points = copy.deepcopy(points_in)
    
    
    n = len(points)
    centers = []
    center_distances = np.zeros(n)
    
    # Erstes Zentrum Zufällig wählen
    if (firstIndex == -1):
        first_center_index = np.random.randint(n)
    # Oder den Vorgegebenen Index benutzen
    else:
        first_center_index = firstIndex
    centers.append(points[first_center_index])
    
    # Neues Centrum bei sich als Zentrum vermerken
    points[first_center_index].isCenter = True
    # Zentum sich selbst zuweisen
    points[first_center_index].cluster = 0
    
    
    # Distanz von jedem Punkt zu dem ersten Zentrum berechnen
    for i, point in enumerate(points):
        # Distanz berechnen
        center_distances[i] = point.dist_to(points[first_center_index])
        
        # Clusterzugehörigkeit setzen
        points[i].cluster = 0
    
    
    # Restliche Zentren berechnen
    for clu in range(1, k):
        # Sucht den Punkt mit dem maximalen Abstand zu allen Zentren
        max_index = np.argmax(center_distances)
        # Diesen neuen Punkt zu der Liste von Zentren hinzufügen
        centers.append(points[max_index])
        
        # Neues Centrum bei sich als Zentrum vermerken
        points[max_index].isCenter = True
        # Zentum sich selbst zuweisen
        points[max_index].cluster = clu
        
        
        # Updaten der minimalen Distanz zu einem Cluster
        for i, point in enumerate(points):
            # Muss Cluster aktuallisiert werden?
            if (center_distances[i] > point.dist_to(points[max_index])):
                # Distanz aktuallisieren
                center_distances[i] = point.dist_to(points[max_index])
                # Cluster Aktuallisieren
                points[i].cluster = clu

    
    return centers, points



# Gonzalez-Variante bei der die Zentren nochmal optimiert werden
def optimized_gonzalez(points, k, start_index = -1):
    # Gonzalez normal machen    
    _, clustered_points = gonzalez(points, k, firstIndex=start_index)
    
    # Zentren optimieren
    optimized_centers, clustered_points = optimize_centers(clustered_points, k)
    
    # Optimierte Punkte und Zentren zurückgeben
    return optimized_centers, clustered_points



























# # # =========== simple k-msr =========== # # #
# Berechent erst optimierten Gonzalez und bildet dann die Summe
def gonzalez_k_msr(points, k, first_index=0):
      
    # Modifizierter Gonzalez machen
    centers, clustered_points = optimized_gonzalez(points, k, start_index=first_index)
    
    
    # Maximale Radien zu jedem Zentrum berechnen
    max_radii = get_radii_of_clustering(clustered_points, k)

    # Summe der Radien bilden
    radii_sum = 0
    for rad in max_radii:
        radii_sum += rad

    return centers, clustered_points, radii_sum


# # # Iterative-Gonzalez for k-msr # # #
# Berechnet auch den optimierten k-msr aber testet für jeden Startindex
def iterative_gonzalez_k_msr(points, k):
    
    # Variablen anlegen
    n = len(points)
    max_sums = np.full(n, -1, dtype=np.float64)
    
    # Jeden Startindex ausprobieren
    for start_index in range(0, n):
    
        # k-msr Gonzalez auf diesem Startindex Ausprobieren
        _, _, radii_sum = gonzalez_k_msr(points, k, first_index=start_index)
        
        # Potentielle Radius-Summe abspeichern
        max_sums[start_index] = radii_sum
        
        
    # Bester Index Finden
    best_start_index = np.argmin(max_sums)
    
    # Diesen Index nochmal zum richtigen Clustering benutzen
    centers, clustered_points, radii_sum_best = gonzalez_k_msr(points, k, 
                                                                  first_index=best_start_index)

    return centers, clustered_points, radii_sum_best
    











# # # =========== Mergender k-msr Algorithmus =========== # # #
# Auf Basis von gonzalez
def merging_k_msr(points, k, use_iterative_gonzalez=True):
    
    # Welchen Basisalgorithmus
    if(use_iterative_gonzalez):
        centers, clustered_points, _ = iterative_gonzalez_k_msr(points, k)
    else:
        centers, clustered_points, _ = gonzalez_k_msr(points, k)
    
    # Die Radien holen um sie später miteinander zu vergleichen
    all_radii = get_radii_of_clustering(clustered_points, k)
    
    # schauen ob man Cluster mergen kann
    merge_clusters(centers, clustered_points, all_radii)
    
    # neue Centrenliste holen
    centers = get_centers_of_points(clustered_points)
    
    # neue Radien berechnen
    all_radii = get_radii_of_clustering(clustered_points, len(centers))
    # -1 rausschmeißen (invaliede Clusternummern)
    for i in range(0, len(all_radii)):
        if (all_radii[i] < 0):
            all_radii[i] = 0
    
    # Aufsummieren     
    radii_sum = sum(all_radii)
    
    return centers, clustered_points, radii_sum


# Merged alle Cluster welche sich mit ihren Radien berühren
def merge_clusters(centers, clustered_points, all_radii):
    
    # Anzahl an Cluster bestimmen
    k = len(all_radii)
    
    newK = k
    
    
    # k mal durchgehen, weil potenziell das selbe Cluster k mal verschmolzen werden kann
    for i in range(0, k):
        # Alle Zentren durchgehen
        for first_index in range(0, k):
            # mit allen anderen Zentren vergleichen
            for second_index in range(0, k):
                # Nicht mit sich selbst verschmelzen
                if (first_index == second_index):
                   continue 
                # Cluster mit negativem Radius existieren nicht mehr (schon mal verschmolzen)
                if (all_radii[first_index] < 0 or all_radii[second_index] < 0):
                   continue 
                
                # Ist der Abstand kleiner als die Summe der Radien
                sum_of_radii = all_radii[first_index] + all_radii[second_index]
                
                # verschmelzen
                if (centers[first_index].dist_to(centers[second_index]) < sum_of_radii):
                    # Die zwei Cluster werden gemerged zu dem Cluster mit ID first_index
                    merge_cluster_by_index(first_index, second_index, clustered_points)
                    
                    # markieren, das dieses Cluster nicht mehr existiert
                    all_radii[second_index] = -1
                    
                    # Anzahl der Cluster vermerken
                    newK -= 1
                    
                    # Zentrum in neuem Cluster optimieren
                    opt_center = optimize_single_center(clustered_points, first_index)
                    
                    # optimiertes Zentrum in Zentrenliste schreiben
                    centers[first_index] = opt_center
                    
                    # Neue Clustergröße vom verschmelztem Cluster abspeichern
                    all_radii[first_index] = get_radius_of_cluster(first_index, clustered_points)
    
    return newK


# Fügt alle Punkte des zweiten Clusters dem ersten zu
def merge_cluster_by_index(first_index, second_index, clustered_points):    
    # Durch jeden Punkt durchgehen
    for i in range(0, len(clustered_points)):
        # Wenn er zu dem zweiten Cluster gehört, wird er dem ersten hinzugefügt
        if (clustered_points[i].cluster == second_index):
            # Clusterzugeörigkeit ändern
            clustered_points[i].cluster = first_index
            # Zentrum des zweiten clusters löschen
            if (clustered_points[i].isCenter):
                print(f"Zentrum deaktiviert. Das Zentrum {second_index} wurde zum Cluster {first_index} hinzugefügt")
                clustered_points[i].isCenter = False
        
    return











# # # =========== Utility Methoden =========== # # #
# Erstellt liste aller Zentren in den Punkten
def get_centers_of_points(clustered_points):
    # Leere Liste erstellen
    new_centers = []
    
    # Durch alle Punkte durchgehen
    for i in range(0, len(clustered_points)):
        # Ist dieser Punkt ein Zentrum?
        if (clustered_points[i].isCenter):
            # Zur Liste hinzufügen
            new_centers.append(clustered_points[i])
    
    return new_centers






# # # Alle Radien der Cluster bestimmen # # #
# Berechne den Radius eines bestimmten Clusters
def get_radius_of_cluster(clusterID, clustered_points):
    n = len(clustered_points)
    
    # Zentrum des CLuster finden
    center = -1
    for i in range(0, n):
        if (clustered_points[i].cluster == clusterID):    
            # Zentrum abspeichern
            if (clustered_points[i].isCenter):
                center = clustered_points[i]
    
    # Fehlerfall kein Zentrum gefunden
    if (center == -1):
        print(f"Das Cluster {clusterID} besitz kein Zentrum")
        return -1
    
    # Abstand zu allen Clustermitgliedern messen
    max_radius = -1
    for i in range(0, n):
        if (clustered_points[i].cluster == clusterID):    
            # Nach neuem Radius prüfen
            if (center.dist_to(clustered_points[i]) > max_radius):
                # Radius aktuallisieren
                max_radius = center.dist_to(clustered_points[i])
    
    return max_radius


# Berechnet alle Radien
def get_radii_of_clustering(clustered_points, k):
    # Liste für die Radien
    radii = []
    
    # Laufvariablen
    clu = 0
    highes_cluster = k
    
    # Durch jedes Cluster gehen und radius suchen
    while (clu < highes_cluster):
        # Radius eines Clusters mit der ID clu berechnen
        radius_of_cluster = get_radius_of_cluster(clu, clustered_points)
        # Diesen Radius der Liste hinzufügen
        radii.append(radius_of_cluster)
        
        # Gelöschtes Cluster ist nicht mehr vorhanden belegt aber den Index
        # dadurch ist ggf highes_cluster größer als k        
        if(radius_of_cluster < 0):
            highes_cluster += 1
        
        clu += 1
    
    return radii







# # # Zentrum eines bestimmten Clusters Optimieren # # #
# Gibt das optimierte Zentrum oder -1 bei leerem Cluster zurück
def optimize_single_center(clustered_points, cluster_ID):
    
    # Alle Punkte aus diesem Cluster als Liste holen
    points_of_cluster = []
    for i in range(len(clustered_points)):
        # Falls der Punkt zum Cluster gehört, hinzufügen
        if(clustered_points[i].cluster == cluster_ID):           
            # Falls er das Zentrum ist, das weglöschen
            if (clustered_points[i].isCenter):
                clustered_points[i].isCenter = False
            
            points_of_cluster.append(clustered_points[i])
    
    # Für den Fall das das Cluster leer ist, dann wurde es schon verschmolzen
    # Es soll übersprungen werden
    if(len(points_of_cluster) == 0):
        return -1
    
    # Erstes Zentrum beliebig wählen
    best_center_index = 0    
    
    # Für alle Punkte Testen ob sie ein Gutes Zentrum sind
    max_dist_of_point = np.full(len(points_of_cluster), -1, dtype=np.float64)
    
    for i in range(0, len(points_of_cluster)):
        # Array wo ich die Distanzen zu jedem Punkt im Cluster speichere
        radii = np.full(len(points_of_cluster), -1, dtype=np.float64)
        
        # Radien Berechnen
        for index, member in enumerate(points_of_cluster):
            radii[index] = points_of_cluster[i].dist_to(member)
    
        # Maximum suchen, da dies die Clustergröße wäre
        max_dist_of_point[i] = max(radii)
 
    # Bestes Zentrum rausfischen
    best_center_index = np.argmin(max_dist_of_point)
    
    # Zentrum wählen
    points_of_cluster[best_center_index].isCenter = True
        
    # Besseres Zentrum zurückgeben
    return points_of_cluster[best_center_index]
    


# Optimiert alle Cluster durch
def optimize_centers(clustered_points, k):
    
    clu = 0
    highest_cluster = k
    
    # Neue Zentren abspeichern
    optimized_centers = []   

    # Durch jedes Cluster gehen und optimieren
    while (clu < highest_cluster):
        # Optimierte Zentrum berechen
        new_center = optimize_single_center(clustered_points, clu)
        
        # Für den Fall das das Cluster leer ist, dann wurde es schon verschmolzen
        # Es soll übersprungen werden
        if(new_center == -1):
            highest_cluster += 1
            clu += 1
            if(highest_cluster > 100):
                break
            continue
        
        # Zentrum einhängen in Liste
        optimized_centers.append(new_center)
        
        # Laufindizes
        clu += 1
        
    return optimized_centers, clustered_points














# # # =========== Douglas Ansatz =========== # # #
# # # =========== Basis Complete Linkage =========== # # #

# naiver Ansatz - überhaupt nicht effizient
def link_clusters(clustered_points, k):
    # erzeuge Cluster für jeden Punkt und merge Cluster bis nur noch k Cluster übrig
    # objective function: r(A u B) - r(A) - r(B)
    clusters = []
    cluster_centers = []
    # n-tes Clusterzentrum entspricht n-tem Cluster
    for i in range(len(clustered_points)):
        clusters.append([clustered_points[i]])
        cluster_centers.append(clustered_points[i])

    for i in range(len(clustered_points)-k-1):
        best_cost_reduction = np.inf
        for A in clusters:
            for B in clusters:
                if A == B:
                    continue
                # objective function
                new_center, r_AuB = get_radius_of_cluster_link_version(A + B)
                center_A, r_A = get_radius_of_cluster_link_version(A)
                center_B, r_B = get_radius_of_cluster_link_version(B)
                cost_reduction = r_AuB - r_A - r_B

                # find best pair of Clusters to merge
                if cost_reduction < best_cost_reduction:
                    best_cost_reduction = cost_reduction
                    to_merge = (A, B)
                    old_center_indices = (clusters.index(A), clusters.index(B))
                    best_center = new_center

        # nur wenn es sich verbessern würde
        if best_cost_reduction >= 0:
            # passe Clustermenge an
            clusters.remove(to_merge[0])
            clusters.remove(to_merge[1])
            clusters.append(to_merge[0] + to_merge[1])

            # passe Zentrenmenge an
            old_centers = (cluster_centers[old_center_indices[0]], cluster_centers[old_center_indices[1]])
            cluster_centers.remove(old_centers[0])
            cluster_centers.remove(old_centers[1])
            cluster_centers.append(best_center)
            
            # Debug print
            print(f"Link Clusters together in iteration {i}")

    return cluster_centers


# von GPT-4 optimierter Code
def optimize_link_clusters(clustered_points, k):
    # Precompute the initial radii and centers for all clusters
    cluster_info = [get_radius_of_cluster_link_version([pt]) for pt in clustered_points]
    clusters = [[i] for i in range(len(clustered_points))]  # Use indices instead of points directly

    while len(clusters) > k:
        best_cost_reduction = np.inf
        to_merge_indices = None
        best_new_info = None
        
        # Debug print
        print(f"Link Clusters in while-Schleife. Len Clusters = {len(clusters)}")

        # Iterate through pairs of clusters
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                cluster_i = clusters[i]
                cluster_j = clusters[j]

                # Calculate merged cluster info
                merged_cluster = [clustered_points[idx] for idx in cluster_i + cluster_j]
                new_center, r_AuB = get_radius_of_cluster_link_version(merged_cluster)
                _, r_A = cluster_info[cluster_i[0]]
                _, r_B = cluster_info[cluster_j[0]]
                cost_reduction = r_AuB - r_A - r_B

                if cost_reduction < best_cost_reduction:
                    best_cost_reduction = cost_reduction
                    to_merge_indices = (i, j)
                    best_new_info = (new_center, r_AuB)
               


        if best_cost_reduction < 0:
            # Merge clusters
            clusters.append(clusters[to_merge_indices[0]] + clusters[to_merge_indices[1]])
            cluster_info.append(best_new_info)  # Add new cluster info

            # Remove old clusters and their info
            for idx in sorted(to_merge_indices, reverse=True):
                del clusters[idx]
                del cluster_info[idx]
            
    # Reconstruct final cluster centers from indices
    final_clusters = [cluster_info[clusters[i][0]] for i in range(len(clusters))]
    return [info[0] for info in final_clusters]


def get_radius_of_cluster_link_version(cluster):
    # suche minimalen Radius und bestes Zentrum für eine Liste von Punkten

    if len(cluster) == 1:
        return cluster[0], 0

    best_radius = {}
    for p1 in cluster:
        best_for_p1 = np.inf
        for p2 in cluster:
            if p1 == p2:
                continue
            radius = p1.dist_to(p2)
            if radius < best_for_p1:
                best_for_p1 = radius

        best_radius[p1] = best_for_p1

    best_center = min(best_radius, key=best_radius.get)
    radius = best_radius[best_center]

    return best_center, radius