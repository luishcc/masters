import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
 
dt = 0.01
tempo = 1200
beta = 8./3
sig = 10
rho = 28
 
x = sp.zeros(tempo)
y = sp.zeros(tempo)
z = sp.zeros(tempo)

x2 = sp.zeros(tempo)
y2 = sp.zeros(tempo)
z2 = sp.zeros(tempo)

x[0] = 0.101
y[0] = 0.101
z[0] = 0.101

x2[0] = 0.1
y2[0] = 0.1
z2[0] = 0.1
 
for i in range(tempo-1):
    x[i+1] = sig*dt*(y[i]-x[i]) +x[i]
    y[i+1] = x[i]*dt*(rho - z[i]) + y[i]*(1-dt)
    z[i+1] = dt*(x[i]*y[i]-beta*z[i]) + z[i]

    x2[i+1] = sig * dt * (y2[i] - x2[i]) + x2[i]
    y2[i+1] = x2[i] * dt * (rho - z2[i]) + y2[i] * (1-dt)
    z2[i+1] = dt * (x2[i] * y2[i] - beta * z2[i]) + z2[i]

 
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot3D(x, y, z, 'b-', linewidth=0.2)
ax.scatter(x[0], y[0], z[0], 'bo')
ax.scatter(x[tempo-1], y[tempo-1], z[tempo-1], 'bo')

ax.plot3D(x2, y2, z2, 'k-', linewidth=0.2)
ax.scatter(x2[0], y2[0], z2[0], c='k')
ax.scatter(x2[tempo-1], y2[tempo-1], z2[tempo-1], c='k')

ax.grid('off')
plt.show()
