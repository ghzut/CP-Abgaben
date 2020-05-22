import numpy as np
import matplotlib.pyplot as plt

n, int, err = np.genfromtxt("A1/build/A1b.txt", unpack=True)

#b)
x = np.linspace(0, np.max(n))
plt.plot(n, int, "b-", label= r'Integral $e^{-x}/x$ Simpsonregel')
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
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.savefig("A1/build/1b_err.pdf")

#c)

n, int, err = np.genfromtxt("A1/build/A1c.txt", unpack=True)

plt.clf()
x = np.linspace(0.1, np.max(n))
plt.plot(n, int, "b-", label= r'Integral $\sin{x}/x$ Mittelpunktsregel')
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
