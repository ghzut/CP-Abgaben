import numpy as np
import matplotlib.pyplot as plt

k1, err1, r1 = np.genfromtxt("A2/build/A2_1.txt",unpack=True)
k2, err2, r2 = np.genfromtxt("A2/build/A2_2.txt",unpack=True)
k3, err3, r3 = np.genfromtxt("A2/build/A2_3.txt",unpack=True)


plt.plot(k1, err1, "r-", label="Exakte Hesse Init.")
plt.plot(k2, err2, "b-", label="Diag. Hesse Init.")
plt.plot(k3, err3, "g-", label="Diag. f(x) Init.")
plt.legend(loc="best")
plt.yscale("log")
plt.xlabel(r'$k$')
plt.ylabel(r'$|b_k|$')
plt.tight_layout()
plt.savefig("A2/build/A2.pdf")

plt.clf()
plt.plot(k1, r1, "r-", label="Exakte Hesse Init.")
plt.plot(k2, r2, "b-", label="Diag. Hesse Init.")
plt.plot(k3, r3, "g-", label="Diag. f(x) Init.")
plt.legend(loc="best")
plt.yscale("log")
plt.xlabel(r'$k$')
plt.ylabel(r'$r$')
plt.tight_layout()
plt.savefig("A2/build/A2r.pdf")



k1, err1, r1 = np.genfromtxt("A2/build/A2_1_l.txt",unpack=True)
k2, err2, r2 = np.genfromtxt("A2/build/A2_2_l.txt",unpack=True)
k3, err3, r3 = np.genfromtxt("A2/build/A2_3_l.txt",unpack=True)

plt.clf()
plt.plot(k1, err1, "r-", label="Exakte Hesse Init.")
plt.plot(k2, err2, "b-", label="Diag. Hesse Init.")
plt.plot(k3, err3, "g-", label="Diag. f(x) Init.")
plt.legend(loc="best")
plt.yscale("log")
plt.xlabel(r'$k$')
plt.ylabel(r'$|b_k|$')
plt.tight_layout()
plt.savefig("A2/build/A2_l.pdf")

plt.clf()
plt.plot(k1, r1, "r-", label="Exakte Hesse Init.")
plt.plot(k2, r2, "b-", label="Diag. Hesse Init.")
plt.plot(k3, r3, "g-", label="Diag. f(x) Init.")
plt.legend(loc="best")
plt.yscale("log")
plt.xlabel(r'$k$')
plt.ylabel(r'$r$')
plt.tight_layout()
plt.savefig("A2/build/A2_lr.pdf")
