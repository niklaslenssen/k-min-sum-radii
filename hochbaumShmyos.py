import numpy as np
import pandas as pd

def hochbaum_shmoys_k_center_from_csv(df, k):
    # Extrahiere die Koordinatenpunkte aus dem DataFrame
    points = df[['x', 'y']].values
    
    # Wähle den ersten Punkt zufällig aus den Eingabepunkten
    centers = [points[np.random.randint(len(points))]]
    
    # Berechne den maximalen Abstand für die initiale Zuordnung
    distances = np.array([np.linalg.norm(point - centers[0]) for point in points])
    max_distance = np.max(distances)

    # Iteriere, um die restlichen k - 1 Zentren zu finden
    while len(centers) < k:
        # Wähle das nächste Zentrum basierend auf dem größten Abstand
        next_center_index = np.argmax(distances)
        next_center = points[next_center_index]
        centers.append(next_center)
        
        # Aktualisiere die Abstände für jeden Punkt zu seinem nächsten Zentrum
        for i, point in enumerate(points):
            distances[i] = min(distances[i], np.linalg.norm(point - next_center))
        
        # Berechne den neuen maximalen Abstand
        max_distance = np.max(distances)
    
    return max_distance