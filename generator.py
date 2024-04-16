import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.datasets import make_blobs

from hochbaumShmyos import hochbaum_shmoys_k_center_from_csv

k = 3

X, y = make_blobs(n_samples=100, centers=k, n_features=2, cluster_std=0.8, random_state=1) # type: ignore

df = pd.DataFrame(X, columns=['x', 'y'])
df['Center'] = y

rmax = hochbaum_shmoys_k_center_from_csv(df, k)

with open("Data/points.csv", "w") as file:
    file.write(f"{rmax}\n")
df.to_csv("Data/points.csv", mode='a', index=False, header=False)
