import numpy as np
import matplotlib.pyplot as plt

def rosenbrock(x, y):
    return (1-x)**2+100*(y-x**2)**2
    
def f(x, y):
    return 1/(1+np.exp(-10*(x*y-3)**2)/(x**2 + y**2))


n, x1, x2, g1, g2, eps = np.genfromtxt("A1/build/gradientenverfahren.txt", unpack = True)

plt.plot(n, eps, "b.")

plt.yscale("log")
plt.xlabel(r"Iterationen $n$")
plt.ylabel(r"$\varepsilon$")

plt.tight_layout()
plt.savefig("A1/build/gradientenverfahren.pdf")

plt.clf()

X, Y = np.meshgrid(x1, x2)
Z = rosenbrock(X, Y)

plt.contour(X, Y, Z, cmap="winter")
plt.colorbar()
plt.plot(x1, x2, 'r.', markersize=1, label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")

plt.tight_layout()
plt.savefig("A1/build/gradientenverfahren_contour.pdf")
plt.clf()

n, x1, x2, g1, g2, eps = np.genfromtxt("A1/build/konjugiert_rosenbrock.txt", unpack = True)

X, Y = np.meshgrid(x1[:-4], x2[:-4])
Z = rosenbrock(X, Y)

plt.contour(X, Y, Z, cmap="winter")
plt.colorbar()
plt.plot(x1[:-4], x2[:-4], 'r.', markersize=1, label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")

plt.tight_layout()
plt.savefig("A1/build/konjugiert_rosenbrock_contour.pdf")
plt.clf()

plt.plot(n[:-4], eps[:-4], "b.")

# plt.yscale("log")
plt.xlabel(r"Iterationen $n$")
plt.ylabel(r"$\varepsilon$")

plt.tight_layout()
plt.savefig("A1/build/konjugiert_rosenbrock.pdf")

plt.clf()

n, x1, x2, g1, g2, eps = np.genfromtxt("A1/build/konjugiert_b1.txt", unpack = True)

X, Y = np.meshgrid(x1, x2)
Z = f(X, Y)

plt.contour(X, Y, Z, cmap="winter")
plt.colorbar()
plt.plot(x1, x2, 'r.', markersize=1, label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")

plt.tight_layout()
plt.savefig("A1/build/b1.pdf")
plt.clf()

n, x1, x2, g1, g2, eps = np.genfromtxt("A1/build/konjugiert_b2.txt", unpack = True)

X, Y = np.meshgrid(x1, x2)
Z = f(X, Y)

plt.contour(X, Y, Z, cmap="winter")
plt.colorbar()
plt.plot(x1, x2, 'r.', markersize=1, label="Schritte")
plt.legend(loc='best')
plt.xlabel(r"$x_1$")
plt.ylabel(r"$x_2$")

plt.tight_layout()
plt.savefig("A1/build/b2.pdf")
plt.clf()

