import numpy as np
import math
h = 1.0E-10
a = 0.999
b = 1.0
N = int(round((b-a)/h))
radius = 6371.0


pts = np.zeros([N+1,1])
pts[0] = a

for i in range(1, N+1):
    pts[i] = a+ i*h

integral1 = 0
integral2 = 0

l = 0
c = 1.0
op = np.zeros([N+1,1])
op[0] = np.exp(-3*radius*math.acos(pts[0])/8.5)*np.polynomial.legendre.legval(pts[0], c)

for i in range(0, N+1):
    op[i] = np.exp(-3*radius*math.acos(pts[i])/8.5)*np.polynomial.legendre.legval(pts[i], c)
    integral1 = integral1 + (op[i] + op[i-1])*h*math.pi

print integral1

l = 70
c = np.zeros(71)
c[-1] = 1.0
op = np.zeros([N+1,1])
op[0] = np.exp(-3*radius*math.acos(pts[0])/8.5)*np.polynomial.legendre.legval(pts[0], c)

for i in range(0, N+1):
    op[i] = np.exp(-3*radius*math.acos(pts[i])/8.5)*np.polynomial.legendre.legval(pts[i], c)
    integral2 = integral2 + (op[i] + op[i-1])*h*math.pi

print integral2
