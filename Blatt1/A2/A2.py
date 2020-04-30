import numpy as np
import matplotlib.pyplot as plt
from table2 import makeTable

m,n = np.genfromtxt("A2/data/py_input_m_n.txt", unpack=True)
x,y = np.genfromtxt("A2/data/py_input_x_y.txt", unpack=True)

makeTable([x,y], r'{$x$} & {$y$}', 'data_tab', ['S[table-format=2.1]' , 'S[table-format=2.1]'], ["%2.1f", "%2.1f"])

def f(x,m,n):
    return m*x+n

x_plot = np.linspace(x.min()-1.,x.max()+1.,1000)

plt.plot(x,y,'rx', label='Data')
plt.plot(x_plot, f(x_plot, m, n), 'b-', label='Regression')
plt.legend(loc='best')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('A2/build/Regression.pdf')
