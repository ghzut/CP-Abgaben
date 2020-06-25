import numpy as np
import matplotlib.pyplot as plt

t, vx, vy, ekin, epot, T = np.genfromtxt("A1/build/aequi1.txt", unpack=True)

plt.plot(t, T)
plt.savefig("A1/build/aequi1_T.pdf")

plt.clf()

plt.plot(t[:50], T[:50])
plt.savefig("A1/build/aequi1_TT.pdf")

plt.clf()

plt.plot(t, ekin)
plt.plot(t, epot)
plt.plot(t, ekin+epot)

plt.savefig("A1/build/aequi1_E.pdf")

plt.clf()

E = ekin+epot

plt.plot(t[:50], ekin[:50])
plt.plot(t[:50], epot[:50])
plt.plot(t[:50], E[:50])

plt.savefig("A1/build/aequi1_EE.pdf")

plt.clf()

plt.plot(t, vx)
plt.plot(t, vy)

plt.savefig("A1/build/aequi1_V.pdf")

plt.clf()

plt.plot(t[:50], vx[:50])
plt.plot(t[:50], vy[:50])

plt.savefig("A1/build/aequi1_VV.pdf")


plt.clf()

t, T = np.genfromtxt("A1/build/messung1.txt", unpack=True)

plt.plot(t, T)

plt.savefig("A1/build/messung1_T.pdf")

plt.clf()

T_mean = np.mean(T)
T_std = np.std(T)

print(T_mean, " +- ", T_std)

plt.clf()

t, T = np.genfromtxt("A1/build/messung001.txt", unpack=True)

plt.plot(t, T)

plt.savefig("A1/build/messung001_T.pdf")

plt.clf()

T_mean = np.mean(T)
T_std = np.std(T)

print(T_mean, " +- ", T_std)

plt.clf()

t, T = np.genfromtxt("A1/build/messung100.txt", unpack=True)

plt.plot(t, T)

plt.savefig("A1/build/messung100_T.pdf")

plt.clf()

T_mean = np.mean(T)
T_std = np.std(T)

print(T_mean, " +- ", T_std)

t, vx, vy, ekin, epot, T = np.genfromtxt("A1/build/aequi_isokinetisch001.txt", unpack=True)

plt.plot(t, ekin)
plt.plot(t, epot)
plt.plot(t, ekin+epot)

plt.savefig("A1/build/aequi1_E.pdf")

plt.clf()

E = ekin+epot

plt.plot(t[:50], ekin[:50])
plt.plot(t[:50], epot[:50])
plt.plot(t[:50], E[:50])

plt.savefig("A1/build/aequi1_EE.pdf")

plt.clf()
