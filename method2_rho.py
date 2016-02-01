import numpy as np
import math
from scipy.special import lpmn

l_range = 100
M = 2000
N = 1000
num_points = 10000
radius = 1.0
#radius = np.float64(6371.0)

# 1km = 0.000157

theta = np.reshape(np.linspace(0, math.pi, N), [1,N])
phi = np.reshape(np.linspace(0, 2*math.pi, M), [M,1])
#theta = np.reshape(np.linspace(0.5759, 0.6196, N), [1,N])
#phi =   np.reshape(np.linspace(2.0420, 2.094,  M), [M,1])

def compute_Ylm(l, theta, phi):
    Y = np.zeros([np.size(phi), np.size(theta), 2*l+1])
    L = np.zeros([l+1, np.size(theta)]) 

    for the in range(0, np.size(theta)):
        Pmn_z, Pmn_dz = lpmn(l, l, np.cos(theta[0, the]))
        L[:,the] = Pmn_z[:,-1]

    for m in range(0, l+1):
        L[m, :] = L[m, :]*((-1)**m)*sqrt_function(l,m)

    L = np.sqrt(1/(2*math.pi))*L

    for m in range(-l,l+1):
        if (m>0):
            Y[:,:, l+m] = np.sqrt(2.0)*np.cos(phi*m)*L[m, :]
        elif (m<0):
            Y[:,:, l+m] = np.sqrt(2.0)*np.sin(phi*m)*L[-m, :]
        else:
            Y[:,:, l+m] = np.ones([M,1])*L[0, :]

    #for i in range(0, 2*l+1):
    #    print Y[:,:, i]
    return Y
    
def sqrt_function(l,m):
    output = np.sqrt(l+0.5)
    for i in range(l-m+1, l+m+1):
        output = output/np.sqrt(i)

    return output

def factorial_function(l,m):
    output = np.float64(1.0)
    for i in range(l-m+1, l+m+1):
        output = output*np.float64(i)
    return np.float64(1.0)/output

X = np.zeros([np.size(phi), np.size(theta)])
A = np.zeros(l_range+1)

h = 2E-4
a = -1.0
#h = 1E-10
#a = 0.99999

g = open('A_l_L%i_H%2.2e_A%2.2e.txt' %(l_range, h, a), 'r')
for l in range(0, l_range+1):
    A[l] = np.float(g.readline())
g.close()

x = np.random.randn((2*l_range+1)*(l_range+1))
x = np.reshape(x, [2*l_range+1, l_range+1])

for l in range(0, l_range+1):
    Y = compute_Ylm(l, theta, phi)

    for m in range(-l, l+1):
        X = X + np.sqrt(A[l])*x[l+m,l]*Y[:,:, l+m]
    if np.mod(l, 20) == 0:
        print 'Finished ', l, ' of ', l_range+1
#for ind in range(0, np.size(phi)):
#    print X[ind,:]

g = open('T_N%i_rho_R1.txt'%(N), 'w')
s = str(np.size(theta)) + '\n' + str(np.size(phi)) + '\n'
g.write(s)
for i in range(0, np.size(theta)):
    for j in range(0, np.size(phi)):
        s = str(X[j,i]) + '\n'
        g.write(s)
g.close()
