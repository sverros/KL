import numpy as np

x = np.linspace(0.99999,1,1000)
A = np.zeros([71, 1000])

for l in range(0, 71):
    c = np.zeros(l+1)
    c[-1] = 1.0
    A[l,:] = np.polynomial.legendre.legval(x,c)

print A
