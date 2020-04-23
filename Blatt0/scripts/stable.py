import matplotlib.pyplot as plt
import numpy as np

x, df = np.genfromtxt("0_1a.txt", unpack=True)

plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil Wurzelproblem')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{f_{st}}$')
plt.xscale('log')
plt.savefig('Wurzel.pdf')

x, df = np.genfromtxt("0_1b.txt", unpack=True)

plt.clf()
plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil 1-cos-Problem')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{f_{st}}$')
plt.xscale('log')
plt.savefig('cos.pdf')

x, df = np.genfromtxt("0_1c.txt", unpack=True)

plt.clf()
plt.plot(x, df, 'b-', label='rel. Fehler stabil-instabil sin-Problem')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$\frac{|f_{is}-f_{st}|}{f_{st}}$')
plt.xscale('log')
plt.savefig('sin.pdf')
