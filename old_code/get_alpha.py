import numpy as np
import math
from scipy.special import sph_harm

def compute_alm(l_range, num_pts):
    alphas = get_alpha_vec(l_range, num_pts)
    alm = [None]*l_range

    #Compute a00
    alm[0] = np.random.normal(0.0, alphas[0])

    # Compute alm
    for l in range(1, l_range):
        vec = np.zeros(l+1, dtype = 'complex')
        num_rand_numbers = 2*l
        rand = np.random.normal(0.0, alphas[l]/2, num_rand_numbers)
        rand0 = np.random.normal(0.0, alphas[l], 1)
        vec[0] = rand0

        for m in range(1, l+1):
            ind = (m-1)*2
            vec[m] = complex(rand[ind], rand[ind+1])

        alm[l] = vec

    return alm
 
def get_alpha_vec(l_range, num_pts):
    # Get the vector of alpha for l = 0, l_range
    avec = np.zeros(l_range, dtype='float')
    for i in range(0, l_range):
        avec[i] = alpha(i,num_pts)

    return avec

def alpha(l, N):
    # Compute alpha for a particular l and number of quadrature points N
    radius = np.float64(6371.0)
    a = -1.0
    b = 1.0
    pts = np.zeros(N+1, dtype='float64')
    pts[0] = a
    int_approx = np.float64(0.0)
    
    for i in range(1, N+1):
        pts[i] = a + (b-a)*i/N
        int_approx += (func(pts[i-1], l) + func(pts[i],l))
    return int_approx*(2.0*math.pi*radius/N)

def func(x, l):
    # Integrand
    rho = np.exp(-3.0*math.acos(x)/8.5)
    c = np.zeros(l+1, dtype = 'float')
    c[-1] = 1.0
    L = np.polynomial.legendre.legval(x,c)

    return rho*L

#--------------------------------------------------------------------------

l_range = 80
quad_precision = 1000
radius = np.float64(6371.0)

alm = compute_alm(l_range, quad_precision)

theta = np.arange(0, 2.0*math.pi, 0.1)
phi = np.arange(0, math.pi, 0.1)

X = np.zeros([np.size(theta), np.size(phi)], dtype = 'float64')

g = open('X.txt', 'w')
s = str(np.size(theta)) + '\n' + str(np.size(phi)) + '\n'
g.write(s)

for t in range(0, np.size(theta)):
    for p in range(0, np.size(phi)):
        th = theta[t]
        ph = phi[p]

        X[t,p] += np.real(alm[0]*sph_harm(0,0,th,p))/radius

        for l in range(1, l_range):
            X[t,p] += alm[l][0] *sph_harm(0, l, th, ph)/radius
            
            for m in range(1, l+1):
                X[t,p] += 2*np.real(alm[l][m]*sph_harm(m,l,th,ph))/radius
        s = str(np.real(X[t,p])) + '\n'
        g.write(s)


g.close()

print np.mean(X), np.std(X)

