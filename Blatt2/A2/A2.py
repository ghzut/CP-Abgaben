import numpy as np
import matplotlib.pyplot as plt

N, t1, t2, t3 = np.genfromtxt("A2/build/times.txt", unpack=True)

plt.plot(N, t1, "b-", label = "Initialisierung")
plt.plot(N, t2, "r-", label = "LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Lösung Mx=b")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A2/build/timers.pdf')

N, t1, t2, t3 = np.genfromtxt("A2/build/times_lin.txt", unpack=True)

plt.cla()
plt.clf()
plt.plot(N, t1, "b-", label = "Initialisierung")
plt.plot(N, t2, "r-", label = "LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Lösung Mx=b")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A2/build/timers_lin.pdf')


