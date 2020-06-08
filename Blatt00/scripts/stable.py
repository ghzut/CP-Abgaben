import matplotlib.pyplot as plt
import numpy as np

x, df = np.genfromtxt("build/0_1a.txt", unpack=True)

plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil Wurzelproblem')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{|f_{is}|}$')
plt.xscale('log')
plt.savefig('build/Wurzel.pdf')

x, df = np.genfromtxt("build/0_1b.txt", unpack=True)

plt.clf()
plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil 1-cos-Problem')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{|f_{is}|}$')
plt.xscale('log')
plt.savefig('build/cos.pdf')

x, df = np.genfromtxt("build/0_1c.txt", unpack=True)

plt.clf()
plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil sin-Problem')
plt.legend(loc='best')
plt.xlabel(r'$\delta$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{|f_{is}|}$')
plt.xscale('log')
plt.savefig('build/sin.pdf')
