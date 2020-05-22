import numpy as np
import matplotlib.pyplot as plt

n, int, int2, err = np.genfromtxt("A1/build/A1b.txt", unpack=True)

#b)
x = np.linspace(0, np.max(n))
plt.plot(n, int, "b-", label= r'Integral 1 $e^{-x}/x$ Simpsonregel')
plt.plot(n+1., int2, "g-", label= r'Integral 2 $e^{-x}/x$ z. Fehlerabsch.')
plt.plot(x, 1.77245385+0*x, "r--")
plt.xlabel(r'$b$')
plt.ylabel(r'$I_b$')
plt.legend(loc="best")
plt.savefig("A1/build/1b.pdf")

plt.clf()
plt.plot(n, err, "b-", label= r'Fehler $e^{-x}/x$ Simpsonregel')
plt.xlabel(r'$b$')
plt.ylabel(r'$|I_{2b}-I_b|$')
plt.legend(loc="best")
plt.yscale("log")
plt.grid()
plt.savefig("A1/build/1b_err.pdf")

#c)

n, int, int, err = np.genfromtxt("A1/build/A1c.txt", unpack=True)

plt.clf()
x = np.linspace(0.1, np.max(n))
plt.plot(n, int, "b-", label= r'Integral 1 $\sin{x}/x$ Mittelpunktsregel')
plt.plot(n+2.5, int2, "g-", label= r'Integral 2 $\sin{x}/x$ z. Fehlerabsch.')
plt.plot(x, np.pi+0*x, "r--")
plt.xlabel(r'$b$')
plt.ylabel(r'$I_b$')
plt.legend(loc="best")
plt.xscale("log")
plt.savefig("A1/build/1c.pdf")

plt.clf()
plt.plot(n, err, "b-", label= r'Fehler $e^{-x}/x$ Simpsonregel')
plt.xlabel(r'$b$')
plt.ylabel(r'$|I_{2b}-I_b|$')
plt.legend(loc="best")
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.savefig("A1/build/1c_err.pdf")
