import numpy as np
import matplotlib.pyplot as plt

# Plotten der zu integrierenden Funktionen
def f1(x):
    return np.exp(-x)/x

def f2(x):
    return x * np.sin(1/x)

x = np.linspace(1,100,1000000)

plt.plot(x, f1(x))
plt.xlabel("x")
plt.ylabel(r"$f_1(x)$")

plt.tight_layout()
plt.savefig("A3/build/f1.pdf")

plt.clf()

x = np.linspace(1e-9,1,1000000)

plt.plot(x, f2(x))
plt.xlabel("x")
plt.ylabel(r"$f_2(x)$")

plt.tight_layout()
plt.savefig("A3/build/f2.pdf")

plt.clf()


# Vergleich der Integrationsroutinen für I1
n_mit, err_mit = np.genfromtxt("A3/build/mittelpunktsregel_f1.txt", unpack = True)
n_tra, err_tra = np.genfromtxt("A3/build/trapezregel_f1.txt", unpack = True)
n_sim, err_sim = np.genfromtxt("A3/build/simpsonregel_f1.txt", unpack = True)

plt.plot(n_mit, err_mit, marker="x", label="Mittelpunktsregel")
plt.plot(n_tra, err_tra, marker="x", label="Trapezregel")
plt.plot(n_sim, err_sim, marker="x", label="Simpsonregel")

plt.xscale('log', basex=2)
plt.yscale("log")

plt.xlabel("n")
plt.ylabel(r"$\epsilon_{\mathrm{rel.}}$")

plt.legend()
plt.tight_layout()
plt.savefig("A3/build/I1.pdf")

plt.clf()


#Vergleich der Integrationsroutinen für I2
n_mit, err_mit = np.genfromtxt("A3/build/mittelpunktsregel_f2.txt", unpack = True)
n_tra, err_tra = np.genfromtxt("A3/build/trapezregel_f2.txt", unpack = True)
n_sim, err_sim = np.genfromtxt("A3/build/simpsonregel_f2.txt", unpack = True)

plt.plot(n_mit, err_mit, marker="x", label="Mittelpunktsregel")
plt.plot(n_tra, err_tra, marker="x", label="Trapezregel")
plt.plot(n_sim, err_sim, marker="x", label="Simpsonregel")

plt.xscale('log', basex=2)
plt.yscale("log")

plt.xlabel("n")
plt.ylabel(r"$\epsilon_{\mathrm{rel.}}$")

plt.legend()
plt.tight_layout()
plt.savefig("A3/build/I2.pdf")

