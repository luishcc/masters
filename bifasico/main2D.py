# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import semiLagrangean as sl
import libmalha as lm
import GMesh as gm
import os
import forces as ff
import random
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cwd = os.getcwd()

# Definição da malha

h = 1
L = 10
fine = 20
coarse = 10
d_fine = 0.2*h
d_coarse = 0.8*h


x_fine = sp.linspace(0, d_fine, fine, endpoint=False)
x_coarse = sp.linspace(d_fine, d_coarse, coarse, endpoint=False)
x_fine2 = sp.linspace(d_coarse, h, fine+1)
x = sp.append(x_fine, x_coarse)
x = sp.append(x, x_fine2)

nodes = len(x)
elem = nodes-1

viscosity = 0.0091
rho_fld = 1000.0
rho_part = 2500.0
grad_p = -120.
dt = 0.01
tempo = 100
g = 9.81
radius = 0.005
diameter = 2 * radius
volume = (4/3.) * sp.pi * radius**3



ien = sp.zeros((elem, 2), dtype='int32')
dx = sp.zeros(elem)
for i in range(0, elem):
    dx[i] = abs(x[i+1] - x[i])
    for j in [0, 1]:
        ien[i, j] = int(i + j + 1)

k = sp.array([[1,-1], [-1,1]])
m = sp.array([[2, 1], [1, 2]])
b = sp.array([1,1])

B = sp.zeros(nodes)
K = sp.zeros((nodes, nodes))
M = sp.zeros((nodes, nodes))


for i in range(elem):
    B[i:i+2] = B[i:i+2] + b * (-grad_p * dx[i] * 0.5)
    M[i:i+2, i:i+2] = M[i:i+2, i:i+2] + m * (dx[i])/6
    K[i:i+2, i:i+2] = K[i:i+2, i:i+2] + k / dx[i]

Mdt = M / dt
LHS = Mdt + viscosity * K

v0 = 0
vh = 0

for i in range(nodes):
    B[i] = B[i] - v0 * LHS[i, 0]
    B[i] = B[i] - vh * LHS[i, elem]


LHS = sp.delete(LHS, 0, axis=0)
LHS = sp.delete(LHS, 0, axis=1)
LHS = sp.delete(LHS, elem-1, axis=0)
LHS = sp.delete(LHS, elem-1, axis=1)
B = sp.delete(B, 0, axis=0)
B = sp.delete(B, elem-1, axis=0)

v_last = sp.zeros(nodes-2)
Mdt_s = sp.delete(Mdt, 0, axis=0)
Mdt_s = sp.delete(Mdt_s, 0, axis=1)
Mdt_s = sp.delete(Mdt_s, elem-1, axis=0)
Mdt_s = sp.delete(Mdt_s, elem-1, axis=1)


n_particule = 4
posx_particule = sp.zeros(n_particule)
posy_particule = sp.zeros(n_particule)


for i in range(n_particule):
    posx_particule[i] = random.randint(1, 10)*0.01
    posy_particule[i] = random.randint(1, 99)*0.01

vx_p = sp.zeros(n_particule)
vy_p = sp.zeros(n_particule)
vy_p_last = sp.zeros(n_particule)
vx_p_last = sp.zeros(n_particule)

constant = dt / ((rho_part + 0.5 * rho_fld) * volume)

f_g = sp.ones(n_particule) * rho_part * volume * g
f_dragx = sp.zeros(n_particule)
f_dragy = sp.zeros(n_particule)
f_massx = sp.zeros(n_particule)
f_massy = sp.zeros(n_particule)
f_lift = sp.zeros(n_particule)

vx_fld_last = sp.zeros(nodes-2)

vy_fld = sp.zeros(nodes)
vy_fld_last = sp.zeros(nodes)

fld_accx = sp.zeros(n_particule)
fld_accy = sp.zeros(n_particule)
vy_fld_last_big = sp.zeros(nodes)

for t in range(tempo):
    RHS = sp.dot(Mdt_s, vx_fld_last) + B
    vx_fld_new = sp.linalg.solve(LHS, RHS)


    vx_fld_last_big = sp.insert(vx_fld_last, 0, v0)
    vx_fld_last_big = sp.insert(vx_fld_last_big, elem, vh)

    velx_last_interp = sp.interpolate.interp1d(x, vx_fld_last_big)
    vely_last_interp = sp.interpolate.interp1d(x, vy_fld_last_big)

    vx_fld_last = sp.copy(vx_fld_new)

    vx_fld = sp.insert(vx_fld_new, 0, v0)
    vx_fld = sp.insert(vx_fld, elem, vh)

    velx_interp = sp.interpolate.interp1d(x, vx_fld)
    vely_interp = sp.interpolate.interp1d(x, vy_fld)


    for i in range(0, n_particule):

        fld_accx[i] = (velx_last_interp(posy_particule[i]) - velx_interp(posy_particule[i])) * dt
        fld_accy[i] = (vely_last_interp(posy_particule[i]) - vely_interp(posy_particule[i])) * dt

        relative_velx = velx_interp(posy_particule[i]) - vx_p[i]
        relative_vely = vely_interp(posy_particule[i]) - vy_p[i]

        f_dragx[i] = 3 * sp.pi * viscosity * 2 * radius * relative_velx
        f_dragy[i] = 3 * sp.pi * viscosity * 2 * radius * relative_vely

        f_massx[i] = 0.5 * rho_fld * volume * fld_accx[i]
        f_massy[i] = 0.5 * rho_fld * volume * fld_accy[i]

        f_lift[i] = 1.61 * sp.sqrt(viscosity * rho_part) * (diameter**2) * abs(relative_velx)

        vy_p[i] = vy_p_last[i] + (f_g[i] + f_lift[i] + f_dragy[i] + f_massy[i]) * constant
        vx_p[i] = vx_p_last[i - 1] + (f_dragx[i] + f_massx[i]) * constant

        print (f_g[i] + f_lift[i] + f_dragy[i] + f_massy[i])

        posx_particule[i] = vx_p[i] * dt + posx_particule[i - 1]
        posy_particule[i] = vy_p[i] * dt + posy_particule[i - 1]

        if posy_particule[i] < 0 or posy_particule[i] > h:
            vy_p[i] = -vy_p[i]
        if posx_particule[i] > L:
            posx_particule[i] = posx_particule[i] - L


    ax = plt.figure(t)
    for i in range(n_particule):
        ax.scatter(posx_particule[i], posy_particule[i], 'k.')
    plt.xlim(-1, 11)
    plt.ylim(-0.1, 1.1)
    plt.savefig('test'+str(t)+'.jpg', format='jpg')


def virtualMass():



    return f


def lift():

    return f

def drag():

    return f
