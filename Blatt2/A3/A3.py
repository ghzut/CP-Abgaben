import numpy as np
import matplotlib.pyplot as plt

N, t1, t2, t3 = np.genfromtxt("A3/build/times.txt", unpack=True, skip_header = 2)

plt.plot(N, t1, "b-", label = "Lösung Mx=b")
plt.plot(N, t2, "r-", label = "Partielle LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Vollständige LU-Zerlegung")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A3/build/timers.pdf')
plt.clf()


N, dev1, dev2 = np.genfromtxt("A3/build/deviations.txt", unpack=True, skip_header = 2)

plt.plot(N, dev1, "r-", label = "Abweichung der Partielle LU-Zerlegung")
plt.plot(N, dev2, "g-", label = "Abweichung der Vollständige LU-Zerlegung")
plt.xlabel(r'$N$')
plt.ylabel(r'Abweichung in Prozent')
plt.yscale("log")
plt.legend(loc='best')
plt.savefig('A3/build/devs.pdf')
