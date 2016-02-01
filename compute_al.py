import numpy as np
import math

l_range = 100

h = 1E-10
radius = 6371.0
a = np.float64(0.99999)
b = np.float64(1.0)
#h = 2E-04
#radius = 1.0
#a = np.float64(-1.0)
#b = np.float64(1.0)

num_points = int(round((b-a)/h))

pts = np.zeros(num_points+1, dtype = 'float64')
pts[0] = a

for i in range(1, num_points+1):
    pts[i] = a+i*h

A_l = np.zeros(l_range+1, dtype = 'float64')

for l in range(0, l_range+1):
    c = np.zeros(l+1)
    c[-1] = 1.0
    op = np.zeros(num_points+1)
    op[0] = np.exp(-3*radius*math.acos(pts[0])/8.5)*np.polynomial.legendre.legval(pts[0],c)
    
    for i in range(1, num_points+1):
        op[i] = np.exp(-3*radius*math.acos(pts[i])/8.5)*np.polynomial.legendre.legval(pts[i],c)
        A_l[l] += (op[i]+op[i-1])*h*math.pi

g = open('A_l_L%i_H%2.2e_A%2.2e.txt' % (l_range, h, a), 'w')

for l in range(0, l_range+1):
    g.write(str(A_l[l]) + '\n')

g.close()

