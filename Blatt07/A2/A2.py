import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from mpl_toolkits import mplot3d

def plot_traje(fname, winkel1, winkel2):
    """Plotte die Datei, die in fname gegeben ist in einer 3D-Darstellung mit
    dem Azimuth(?)-Winkel winkel."""
    path = "A2/build/"+fname
    Y = np.genfromtxt(path)
    print("Plot "+fname)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.view_init(winkel1, winkel2)
    p = ax.scatter3D(Y[1, :], Y[2, :], Y[3, :], c=Y[0, :], cmap='inferno')
    cbar = fig.colorbar(p)
    cbar.set_label(r'Zeitschritte $t_n$', rotation=270, labelpad=10)
    destination = "A2/build/Plots/"+fname
    fig.savefig(destination.split('.')[0]+'.pdf')
    fig.clf()


def plot_energie(fname):
    print("Plot "+fname)
    df = pd.read_csv(fname, decimal='.',
                     delimiter=' ')
    fig = plt.figure()
    plt.plot(df.Zeit, df.Energie, 'r-')
    plt.xlabel("Zeit")
    plt.ylabel("Energie")
    fig.savefig("A2/build/Plots/Energie.pdf")
    fig.clf()


plot_traje("aufg1_a1.txt", 20, 20)
plot_traje("aufg1_a2.txt", 20, 20)
plot_traje("aufg1_a3.txt", 20, 20)
plot_traje("aufg1_a4.txt", 20, 20)
plot_energie("A2/build/aufg1_b.txt")
