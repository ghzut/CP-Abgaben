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



t, vx, vy, ekin, epot, T = np.genfromtxt("A1/build/aequi001.txt", unpack=True)

plt.plot(t, T)
plt.savefig("A1/build/aequi001_T.pdf")

plt.clf()

plt.plot(t[:200], T[:200])
plt.savefig("A1/build/aequi001_TT.pdf")

plt.clf()

plt.plot(t, ekin)
plt.plot(t, epot)
plt.plot(t, ekin+epot)

plt.savefig("A1/build/aequi001_E.pdf")

plt.clf()

E = ekin+epot

plt.plot(t[:200], ekin[:200])
plt.plot(t[:200], epot[:200])
plt.plot(t[:200], E[:200])

plt.savefig("A1/build/aequi001_EE.pdf")

plt.clf()

plt.plot(t, vx)
plt.plot(t, vy)

plt.savefig("A1/build/aequi001_V.pdf")

plt.clf()

plt.plot(t[:200], vx[:200])
plt.plot(t[:200], vy[:200])

plt.savefig("A1/build/aequi001_VV.pdf")



t, vx, vy, ekin, epot, T = np.genfromtxt("A1/build/aequi100.txt", unpack=True)

plt.plot(t, T)
plt.savefig("A1/build/aequi100_T.pdf")

plt.clf()

plt.plot(t[:200], T[:200])
plt.savefig("A1/build/aequi100_TT.pdf")

plt.clf()

plt.plot(t, ekin)
plt.plot(t, epot)
plt.plot(t, ekin+epot)

plt.savefig("A1/build/aequi100_E.pdf")

plt.clf()

E = ekin+epot

plt.plot(t[:200], ekin[:200])
plt.plot(t[:200], epot[:200])
plt.plot(t[:200], E[:200])

plt.savefig("A1/build/aequi100_EE.pdf")

plt.clf()

plt.plot(t, vx)
plt.plot(t, vy)

plt.savefig("A1/build/aequi100_V.pdf")

plt.clf()

plt.plot(t[:200], vx[:200])
plt.plot(t[:200], vy[:200])

plt.savefig("A1/build/aequi100_VV.pdf")