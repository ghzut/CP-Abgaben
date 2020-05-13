import numpy as np
import matplotlib.pyplot as plt

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
plt.savefig("A3/build/f1.pdf")

plt.clf()

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
plt.savefig("A3/build/f2.pdf")

