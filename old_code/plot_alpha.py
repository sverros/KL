import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
#from mayavi import mlab
import matlab.engine

g = open('X_2.txt', 'r')
X_dim = int(g.readline())
Y_dim = int(g.readline())

X = np.zeros([X_dim, Y_dim], dtype = 'float64')
for i in range(0, X_dim):
    for j in range(0, Y_dim):
        X[i,j] = float(g.readline())
g.close()
print np.mean(X), np.std(X)

theta = np.linspace(0, np.pi, num = X_dim)
phi = np.linspace(0, 2*np.pi, num = Y_dim)




#fig = plt.figure(figsize = (10,10))
#map = plt.imshow(np.transpose(X))
#ch=plt.colorbar(map, shrink=0.7)
#plt.savefig('figures/X.png')

#(n, m) = (X_dim, Y_dim)
#
## Meshing a unit sphere according to n, m
#theta = np.linspace(0, 2 * np.pi, num=n, endpoint=False)
#phi = np.linspace(np.pi * (-0.5 + 1./(m+1)), np.pi*0.5, num=m, endpoint=False)
#theta, phi = np.meshgrid(theta, phi)
#theta, phi = theta.ravel(), phi.ravel()
#theta = np.append(theta, [0.])
#phi = np.append(phi, [np.pi*0.5])
#mesh_x, mesh_y = ((np.pi*0.5 - phi)*np.cos(theta), (np.pi*0.5 - phi)*np.sin(theta))
#triangles = mtri.Triangulation(mesh_x, mesh_y).triangles
#x, y, z = np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta), np.sin(phi)
#
## Defining a custom color scalar field
#vals = np.zeros(X_dim*Y_dim+1)
#vals[:X_dim*Y_dim] = np.reshape(X, [X_dim*Y_dim])#np.sin(6*phi) * np.sin(3*theta)
#colors = np.mean(vals[triangles], axis=1)
#
## Plotting
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#cmap = plt.get_cmap('jet')
#triang = mtri.Triangulation(x, y, triangles)
#collec = ax.plot_trisurf(triang, z, cmap=cmap, shade=False, linewidth=0.)
#collec.set_array(colors)
#collec.autoscale()
#plt.show()


