import numpy as np
import math
from scipy.special import sph_harm, lpmv

##################
# Set parameters
##################

l_range = 100             # maximum l
quad_precision = 10000    # number of intervals in quadrature
radius = np.float64(1.0)  # radius of sphere

h = 0.001                  # grid spacing
theta = np.arange(0, math.pi, h)
phi = np.arange(0,2*math.pi,  h)
write_X = True            # write output matrix to file

#theta = np.array([0.5, 0.6, 0.7])
#phi = np.array([0.75, 0.8, 0.85])
#theta = np.arange(0.5759, 0.6196, h)
#phi = np.arange(2.0420,2.094,  h)

######################
# Function definitions
######################

def get_approx(theta, phi, l_range, quad_precision, radius):
    ''' 
    Get the matrix X
    '''
    alpha_vec = get_alpha_vec(l_range, quad_precision, radius)
    output = np.zeros([np.size(theta), np.size(phi)], dtype = 'float64')
    
    # For testing with MATLAB
    #rand_nums = np.ones((l_range+1)*(l_range+2)*np.size(output))
    rand_nums = np.random.randn((l_range*(l_range+1)+l_range+1))


    coef1 = (0.5/np.sqrt(math.pi))*np.sqrt(alpha_vec[0])
    coef2 = np.zeros(l_range, dtype = 'float64')
    coef3 = np.zeros([l_range, l_range], dtype = 'float64')

    for t in range(0, np.size(theta)):
        if np.mod(t,20) == 0:
            print 'Finished with ', t, 'of ', np.size(theta)
        th = theta[t]

        for l in range(1, l_range):
            coef2[l] =  np.sqrt(alpha_vec[l])*L_fun(th, l, 0)
            for m in range(1, l+1):
                coef3[l,m] = np.sqrt(2*alpha_vec[l])*L_fun(th, l, m)

        for p in range(0, np.size(phi)):
            count = 0
            ph = phi[p]
            output[t,p] += coef1*rand_nums[count]
            count += 1

            for l in range(1, l_range):
                output[t,p] += coef2[l]*rand_nums[count]
                count += 1
                
                for m in range(1, l+1):
                    output[t,p] += (rand_nums[count]*np.cos(m*ph) + \
                                      rand_nums[count+1]*np.sin(m*ph))*coef3[l,m]
                    count += 2
    return output/radius

def L_fun(theta, l, m):
    '''
    Returns the normalized associated legendre with l,m, and x = cos(theta)
    '''
    if (l+m) >= 170:
        return np.float64(0.0)
    else:
        output = np.sqrt((2*l+1)/(4*math.pi))*\
            np.sqrt(np.float64(math.factorial(l-m))/np.float64(math.factorial(l+m)))
        l_out = lpmv(m,l,np.cos(theta))
        return output*l_out
    
def get_alpha_vec(l_range, num_pts, radius):
    '''
    Gets the vector of alpha_l for l = 0, ..., l_range
    '''
    avec = np.zeros(l_range, dtype='float')
    for i in range(0, l_range):
        avec[i] = alpha(i,num_pts, radius)
        if avec[i] < 0:
            avec[i] = 0.0
    return avec

def alpha(l, N, radius):
    '''
    Returns the value of alpha_l
    Compute alpha for a particular l and number of quadrature points N
    '''
    
    a = -1.0
    b = 1.0
    pts = np.zeros(N+1, dtype='float64')
    pts[0] = a
    int_approx = np.float64(0.0)
    
    for i in range(1, N+1):
        pts[i] = a + (b-a)*i/N
        int_approx += (func(pts[i-1], l) + func(pts[i],l))
    
    return int_approx*(2.0*math.pi*radius/N)
    #return (np.float64(l+1))**(-3)

def func(x, l):
    '''
    Returns the integrand for quadrature: rho(arccos(x))*L_p(x)
    '''
    rho = np.exp(-3.0*math.acos(x)/8.5)
    c = np.zeros(l+1, dtype = 'float')
    c[-1] = 1.0
    L = np.polynomial.legendre.legval(x,c)

    return rho*L

#####################
# Function calls
#####################
X = get_approx(theta, phi, l_range, quad_precision, radius)

if write_X == True:
    g = open('X_l%i_h%i.txt'%(l_range, h), 'w')
    s = str(np.size(theta)) + '\n' + str(np.size(phi)) + '\n'
    g.write(s)
    for i in range(0, np.size(theta)):
        for j in range(0, np.size(phi)):
            s = str(X[i,j]) + '\n'
            g.write(s)
    g.close()
