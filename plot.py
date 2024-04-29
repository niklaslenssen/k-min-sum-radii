import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def plot_cluster():
    points = pd.read_csv('Data/cluster.csv', header=None, names=['x', 'y', 'Center'])
    balls = pd.read_csv('Data/balls.csv', header=None, names=['x', 'y', 'radius'])

    plt.scatter(points['x'], points['y'], c=points['Center'])
    ax = plt.gca()

    for _, ball in balls.iterrows():
        kreis = Circle((ball['x'], ball['y']), ball['radius'], fill=False, edgecolor='black')
        ax.add_patch(kreis)
        ax.plot(ball['x'], ball['y'], '+', color='black')


    plt.axis('equal')
    plt.savefig('Data/Plots/schmidt.png')