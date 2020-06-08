import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

N, t1, t2, t3 = np.genfromtxt("A2/build/times.txt", unpack=True)

plt.plot(N, t1, "b-", label = "Initialisierung")
plt.plot(N, t2, "r-", label = "LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Lösung Mx=b")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A2/build/timers.pdf')

N, t1, t2, t3 = np.genfromtxt("A2/build/times_lin.txt", unpack=True)

plt.cla()
plt.clf()
plt.plot(N, t1, "b-", label = "Initialisierung")
plt.plot(N, t2, "r-", label = "LU-Zerlegung")
plt.plot(N, t3, "g-", label = "Lösung Mx=b")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig('A2/build/timers_lin.pdf')

a,b = opt.curve_fit(lambda x,a,b,c,d: a*x**3+b*x**2+c*x+d,  N,  t1)

plt.cla()
plt.clf()
x = np.linspace(200,2000,10000)
plt.plot(N[200:], t1[200:], "kx", label = "Initialisierung Daten")
plt.plot(x, a[0]*x**3+a[1]*x**2+a[2]*x+a[3],'b-', label = "Initialisierung Fit")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig("A2/build/t1_200.pdf")

c,d = opt.curve_fit(lambda x,a,b,c,d: a*x**3+b*x**2+c*x+d,  N,  t2)

plt.cla()
plt.clf()
plt.plot(N[200:], t2[200:], "kx", label = "LU-Zerlegung Daten")
plt.plot(x, c[0]*x**3+c[1]*x**2+c[2]*x+c[3],'b-', label = "LU-Zerlegung Fit")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig("A2/build/t2_200.pdf")

e,f = opt.curve_fit(lambda x,a,b,c,d: a*x**3+b*x**2+c*x+d,  N,  t3)

plt.cla()
plt.clf()
plt.plot(N[200:], t3[200:], "kx", label = "Lösung Mx = b Daten")
plt.plot(x, e[0]*x**3+e[1]*x**2+e[2]*x+e[3],'b-', label = "Lösung Mx = b Fit")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r'$N$')
plt.ylabel(r'$t/s$')
plt.legend(loc='best')
plt.savefig("A2/build/t3_200.pdf")
