import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.datasets import make_blobs

from hochbaumShmyos import hochbaum_shmoys_k_center_from_csv

# Anzahl der Punkte
n_samples = 100

# Anzahl der Cluster
k = 3

# Dimension der Punkte
n_features = 2

# Standardabweichung der Punkte zu ihren Centren
cluster_std_global = 0.3

# Bereich der Centren
center_box = (0, 50)

# Seed
random_state = 20

# Generiere Daten mit einer globalen Standardabweichung.
X, y = make_blobs(n_samples=n_samples, centers=k, n_features=n_features,
                  cluster_std=cluster_std_global, center_box=center_box, random_state=random_state)

# Individuelles Rauschen für jeden Cluster definieren.
cluster_stds = [5, 5, 5]

# Füge unterschiedliches Rauschen zu jedem Cluster hinzu.
for i in range(k):
    cluster_points = X[y == i]
    adjusted_noise = np.random.normal(0, cluster_stds[i], cluster_points.shape)
    X[y == i] = cluster_points + adjusted_noise


# Speichern der Punkte in einer CSV Datei.
df = pd.DataFrame(X, columns=['x', 'y'])
df['Center'] = y

rmax = hochbaum_shmoys_k_center_from_csv(df, k)

with open("Data/points.csv", "w") as file:
    file.write(f"{rmax}\n")
df.to_csv("Data/points.csv", mode='a', index=False, header=False)
