import numpy as np
import matplotlib.pyplot as plt


n, a, b, c = np.genfromtxt("A2/build/intervallhalbierung.txt", unpack=True)

plt.plot(n[50:], a[50:], "rx", label=r"$a$")
plt.plot(n[50:], b[50:], "bx", label=r"$b$")
plt.plot(n[50:], c[50:], "gx", label=r"$c$")

plt.xlabel(r"Anzahl der Interationen")
plt.ylabel(r"$x$")

plt.tight_layout()
plt.legend()
plt.savefig("A2/build/intervallhalbierung2.pdf")

plt.clf()

plt.plot(n[0:50], a[0:50], "rx", label=r"$a$")
plt.plot(n[0:50], b[0:50], "bx", label=r"$b$")
plt.plot(n[0:50], c[0:50], "gx", label=r"$c$")


plt.xlabel(r"Anzahl der Interationen")
plt.ylabel(r"$x$")

plt.tight_layout()
plt.legend()
plt.savefig("A2/build/intervallhalbierung.pdf")

plt.clf()


n, x_0 = np.genfromtxt("A2/build/newton.txt", unpack=True)

plt.plot(n,x_0, "rx")

plt.xlabel(r"Anzahl der Interationen")
plt.ylabel(r"$x_0$")

plt.tight_layout()
plt.savefig("A2/build/newton.pdf")

plt.clf()