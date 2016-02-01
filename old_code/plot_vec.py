import numpy as np
import matplotlib.pyplot as plt

g = open('saved_data/alpha_vec.txt', 'r')

alpha_vec = np.zeros(80)
comp_vec1 = np.zeros(80)
comp_vec2 = np.zeros(80)
for l in range(0, 80):
    alpha_vec[l] = float(g.readline())
    comp_vec1[l] = (l+1)**(-2)
    comp_vec2[l] = (l+1)**(-3)
g.close()

print 'I am here'

plt.hold(True)

plt.plot(comp_vec1[2:], label = 'alpha = 2')
plt.plot(comp_vec2[2:], label = 'alpha = 3')
plt.plot(alpha_vec[2:], label = 'rho')
plt.legend()
plt.savefig('figures/Convergence>2.png')
plt.clf()


plt.plot(comp_vec1[:3], label = 'alpha = 2')
plt.plot(comp_vec2[:3], label = 'alpha = 3')
plt.plot(alpha_vec[:3], label = 'rho')
plt.legend()
plt.savefig('figures/Convergence<3.png')


