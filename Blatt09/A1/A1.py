import numpy as np
import matplotlib.pyplot as plt

h, m = np.genfromtxt("A1/build/mc.txt", unpack=True)

def f(x):
    return np.tanh(x)

x = np.linspace(-5,5, 10000)
plt.plot(x, f(x), "r-", label=r"$\tanh{(\beta H)}$")

plt.plot(h, m, "b-", label=r"Monte-Carlo")

plt.xlabel(r"$H$")
plt.ylabel(r"$m$")

plt.legend()
plt.tight_layout()
plt.savefig("A1/build/vergleich.pdf")

plt.clf()

diff = np.zeros(10000)

for i in range(0,10000-1):
    diff[i] = (f(x[i]) - m[i])/f(x[i])

plt.plot(x, diff, "b-")

plt.xlabel(r"$H$")
plt.ylabel(r"$\frac{m_{\mathrm{analytisch}} - m_{\mathrm{MC}}}{m_{\mathrm{analytisch}}}$")

plt.tight_layout()
plt.savefig("A1/build/rel_differenz.pdf")
