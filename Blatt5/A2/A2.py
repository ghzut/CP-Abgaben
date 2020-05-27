import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

x, phi = np.genfromtxt("A2/build/ausserhalb_a.txt", unpack=True)

def f(x, a):
    return a/x



xx = np.linspace(1.1, 8, 100000)
param, cov = curve_fit(f, x, phi)
plt.plot(xx, f(xx, *param), color='b', linestyle='-', label=r"$\frac{1}{x}$ Näherung")


# plt.plot(xx, f(xx), "b-")
plt.plot(x, phi, "rx", label=r"$\tilde{\Phi}(x)$")

plt.xlabel(r"$x$")
plt.ylabel(r"$\Phi$")

plt.legend()
plt.tight_layout()
plt.savefig("A2/build/ausserhalb_a.pdf") 

plt.clf()

x, phi = np.genfromtxt("A2/build/innerhalb_a.txt", unpack=True)

plt.plot(x, phi, "rx")

plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{\Phi}(x)$")

plt.tight_layout()
plt.savefig("A2/build/innerhalb_a.pdf") 


plt.clf()

x, phi = np.genfromtxt("A2/build/ausserhalb_b.txt", unpack=True)

def g(x, a):
    return a/x**2


plt.plot(x, phi, "rx", label=r"$\tilde{\Phi}(x)$")

plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{\Phi}(x)$")

plt.tight_layout()
plt.savefig("A2/build/ausserhalb_b.pdf") 

plt.clf()

xx = np.linspace(1.1, 8, 100000)
param, cov = curve_fit(g, x, phi)
plt.plot(xx, g(xx, *param), color='b', linestyle='-', label=r"$\frac{1}{x^2}$ Näherung")

plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{\Phi}(x)$")

plt.legend()
plt.tight_layout()
plt.savefig("A2/build/ausserhalb_b_fit.pdf") 

plt.clf()

x, phi = np.genfromtxt("A2/build/innerhalb_b.txt", unpack=True)

plt.plot(x, phi, "rx")

plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{\Phi}(x)$")

plt.tight_layout()
plt.savefig("A2/build/innerhalb_b.pdf") 