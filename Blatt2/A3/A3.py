import numpy as np
import matplotlib.pyplot as plt

N, t1, t2, t3 = np.genfromtxt("A3/build/times.txt", unpack=True, skip_header = 2)

plt.plot(N, t1, "b-", label = "Lösung Mx=b")
plt.plot(N, t2, "r-", label = "Partielle LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Vollständige LU-Zerlegung")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A3/build/timers.pdf')
