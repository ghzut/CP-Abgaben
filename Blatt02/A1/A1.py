import numpy as np
import matplotlib.pyplot as plt

x = ["bild", "A10", "A20", "A50"]

for xn in x:
    mat = np.genfromtxt("A1/build/"+xn+".txt", unpack=True)

    plt.clf()
    plt.gray()
    plt.imsave("A1/build/"+xn+".pdf", mat)
