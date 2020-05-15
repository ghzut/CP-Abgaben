import numpy as np
import matplotlib.pyplot as plt

a,b,c,d,e,f,g,h = np.genfromtxt("A2/build/spektrum.txt", unpack=True)


plt.plot(1,a,'rx')
