# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

R = np.linspace(0, 0.5, 100)
h = 5
u = np.linspace(0,  2*np.pi, 100)

x = np.outer(R, np.cos(u))
y = np.outer(R, np.sin(u))

print len(y), len(x)

z = np.zeros((len(x), len(y)))

for i in range(len(x)):
	for j in range(len(y)):
		z[i,j]  = -2*(x[i,j]**2) - 2*(y[i,j]**2) + 2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x,y,z) # z in case of disk which is parallel to XY plane is constant and you can directly use h
plt.show()
