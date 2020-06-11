import numpy as np
import matplotlib.pyplot as plt

h, i, tol = np.genfromtxt("A1/build/A1_0.txt", unpack=True)

i_10 = i[0:9]
tol_10 = tol[0:9]
plt.plot(i_10, tol_10, 'r.')
plt.plot(i_10, tol_10, 'r-', label=r'$h=1$')
tol_10 = tol[10:19]
plt.plot(i_10, tol_10, 'k.')
plt.plot(i_10, tol_10, 'k-', label=r'$h=10^{-1}$')
tol_10 = tol[20:29]
plt.plot(i_10, tol_10, 'b.')
plt.plot(i_10, tol_10, 'b-', label=r'$h=10^{-2}$')
tol_10 = tol[30:39]
plt.plot(i_10, tol_10, 'g.')
plt.plot(i_10, tol_10, 'g-', label=r'$h=10^{-3}$')
tol_10 = tol[40:49]
plt.plot(i_10, tol_10, 'y.')
plt.plot(i_10, tol_10, 'y-', label=r'$h=10^{-4}$')
plt.xlabel(r'$i$')
plt.ylabel(r'$|r_i-r_0|$')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('A1/build/ab_plot.pdf')

h, i, tol = np.genfromtxt("A1/build/A1_1.txt", unpack=True)

i_10 = i[0:9]
tol_10 = tol[0:9]
plt.clf()
plt.plot(i_10, tol_10, 'r.')
plt.plot(i_10, tol_10, 'r-', label=r'$h=1$')
tol_10 = tol[10:19]
plt.plot(i_10, tol_10, 'k.')
plt.plot(i_10, tol_10, 'k-', label=r'$h=10^{-1}$')
tol_10 = tol[20:29]
plt.plot(i_10, tol_10, 'b.')
plt.plot(i_10, tol_10, 'b-', label=r'$h=10^{-2}$')
tol_10 = tol[30:39]
plt.plot(i_10, tol_10, 'g.')
plt.plot(i_10, tol_10, 'g-', label=r'$h=10^{-3}$')
tol_10 = tol[40:49]
plt.plot(i_10, tol_10, 'y.')
plt.plot(i_10, tol_10, 'y-', label=r'$h=10^{-4}$')
tol_10 = tol[50:59]
plt.plot(i_10, tol_10, 'm.')
plt.plot(i_10, tol_10, 'm-', label=r'$h=10^{-5}$')
plt.xlabel(r'$i$')
plt.ylabel(r'$|r_i-r_0|$')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('A1/build/ab_plot2.pdf')

k,E = np.genfromtxt("A1/build/A1_c_0.txt", unpack=True)
E-=E[0]
plt.clf()
plt.plot(k, E, 'r.')
plt.plot(k, E, 'r-', label='Abweichung der Gesamtenergie, $h=10^{-3}$')
k,E = np.genfromtxt("A1/build/A1_c_1.txt", unpack=True)
E-=E[0]
plt.plot(k, E, 'b.')
plt.plot(k, E, 'b-', label='Abweichung der Gesamtenergie, $h=10^{-4}$')
plt.xlabel(r'$i$')
plt.ylabel(r'$\Delta E$')
plt.legend(loc='best')
plt.savefig('A1/build/E.pdf')
