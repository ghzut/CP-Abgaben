import matplotlib.pyplot as plt
import numpy as np

# Aufgabenteil b)

t_1f, e_1f = np.genfromtxt('A2/build/15kbt-b-fest.txt', unpack=True)
t_1z, e_1z = np.genfromtxt('A2/build/15kbt-b-zufall.txt', unpack=True)

t_2f, e_2f = np.genfromtxt('A2/build/22kbt-b-fest.txt', unpack=True)
t_2z, e_2z = np.genfromtxt('A2/build/22kbt-b-zufall.txt', unpack=True)

t_3f, e_3f = np.genfromtxt('A2/build/3kbt-b-fest.txt', unpack=True)
t_3z, e_3z = np.genfromtxt('A2/build/3kbt-b-zufall.txt', unpack=True)

plt.plot(t_1z, e_1z, 'g-', label=r'$1.5 k_{\mathrm{B}}T$, Zufall')
plt.plot(t_1f, e_1f, 'g--', label=r'$1.5 k_{\mathrm{B}}T$, fest')

plt.plot(t_2z, e_2z, 'r-', label=r'$2.27 k_{\mathrm{B}}T$, Zufall')
plt.plot(t_2f, e_2f, 'r--', label=r'$2.27 k_{\mathrm{B}}T$, fest')

plt.plot(t_3z, e_3z, 'b-', label=r'$3 k_{\mathrm{B}}T$, Zufall')
plt.plot(t_3f, e_3f, 'b--', label=r'$3 k_{\mathrm{B}}T$, fest')


plt.xlabel(r'Simulationszeit $t$')
plt.ylabel(r'mittlere Energie pro Spin $e(t)$')
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('A2/build/plot_b.pdf')
plt.close()
