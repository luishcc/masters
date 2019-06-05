# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
from scipy.linalg import solve
import os
import sys
import csv
import matrix as ma
import matplotlib.pyplot as plt

cwd = os.getcwd()


# --------------------------------------------------
#   Problem Parameters

dt = 10
tempo = 1000


# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_din = 0.8e-1
viscosity_kin = viscosity_din / rho_fld
grad_p = -12

# boundary
v0 = 0
vh = 0



# --------------------------------------------------
#   Fluid Mesh and Matrices Generation

h = 1   # Duct height
L = 5*h     # Duct length
fine = 20
coarse = 500
d_fine = 0.15*h
d_coarse = 0.8*h

x_fine = sp.linspace(0, d_fine, fine, endpoint=False)
x_coarse = sp.linspace(d_fine, d_coarse, coarse, endpoint=False)
x_fine2 = sp.linspace(d_coarse, h, fine+1)
x = sp.append(x_fine, x_coarse)
x = sp.append(x, x_fine2)

x = sp.linspace(0, h, 100)

nodes = len(x)
elem = nodes-1

ien = sp.zeros((elem, 2), dtype='int32')
dx = sp.zeros(elem)
for i in range(0, elem):
    dx[i] = abs(x[i+1] - x[i])
    for j in [0, 1]:
        ien[i, j] = int(i + j)


k = sp.array([[1, -1], [-1, 1]])
m = sp.array([[2, 1], [1, 2]])
grad = sp.array([[-0.5, 0.5], [-0.5, 0.5]])
Q = sp.ones(nodes) * (-grad_p / rho_fld)
M = sp.zeros((nodes, nodes))
G = sp.zeros((nodes, nodes))

for e in range(elem):
    for i_local in range(2):
        i_global = ien[e, i_local]
        for j_local in range(2):
            j_global = ien[e, j_local]
            M[i_global, j_global] += m[i_local, j_local] * (dx[e]) / 6.
            G[i_global, j_global] += grad[i_local, j_local]

Mdt = M / dt
M_inv = sp.linalg.inv(M)


# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)
du_dy = sp.zeros(nodes)
viscosity_turb = sp.zeros(nodes)
y_plus = sp.zeros(nodes)
u_plus = sp.zeros(nodes)
u_plus_a = sp.zeros(nodes)
t_wall = sp.zeros(nodes)
u_t = sp.zeros(nodes)
lc = sp.zeros(nodes)



for i in range(nodes):
    if 0 <= x[i] <= 0.1 * h:
        lc[i] = 0.41 * x[i]
    elif h >= x[i] >= 0.9 * h:
        lc[i] = 0.41 * (h - x[i])
    else:
        lc[i] = 0.41 * 0.1 * h


def getViscosity(_u, _lc):
    du_dy = sp.dot(G, _u)
    for i in range(nodes):
        viscosity_turb[i] = _lc[i]**2 * abs(du_dy[i])

    return viscosity_turb, du_dy



with open('turbulent.csv', mode='w') as turbulent:
    writer = csv.writer(turbulent, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow([viscosity_turb, u_plus, u_plus_a, y_plus])

with open('flowResultT.csv', mode='w') as flowResultT:
    writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------



for t in range(tempo):
    print t, " / ", tempo
    # with open('flowResultT.csv', mode='a') as flowResultT:
    #     writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     writer.writerow(u_last)

    viscosity_turb, du_dy = getViscosity(u_last, lc)
    viscosity_turb = sp.dot(M_inv, viscosity_turb)
    K = sp.zeros((nodes, nodes))
    for e in range(elem):
        v1 = ien[e, 0]
        v2 = ien[e, 1]
        visc = (viscosity_turb[v1] + viscosity_turb[v2]) * 0.5
        for i_local in range(2):
            i_global = ien[e, i_local]
            for j_local in range(2):
                j_global = ien[e, j_local]
                K[i_global, j_global] += k[i_local, j_local] * (viscosity_kin + visc) * (1. / dx[e])

    LHS = Mdt + K
    cc = sp.zeros(nodes)
    LHS_copy = sp.copy(LHS)
    t_wall2 = viscosity_din * abs(du_dy[0])
    u_t2 = sp.sqrt((t_wall2 / rho_fld))
    for i in range(nodes):
        t_wall[i] = viscosity_kin * abs(du_dy[i])
        # u_t[i] = sp.sqrt((t_wall[i] / rho_fld))
        # y_plus[i] = (u_t[i] * x[i]) / viscosity_kin
        y_plus[i] = (u_t2 * x[i]) / viscosity_kin
        # u_plus[i] = u_last[i] / u_t[i]
        u_plus[i] = u_last[i] / u_t2
        u_plus_a[i] = (1./0.41) * sp.log(y_plus[i]) + 5.5
        cc[i] -= v0 * LHS_copy[i, 0]
        cc[i] -= vh * LHS_copy[i, -1]
        LHS[0, i] = 0
        LHS[i, 0] = 0
        LHS[-1, i] = 0
        LHS[i, -1] = 0

    LHS[0, 0] = 1
    LHS[-1, -1] = 1

    u_last[0] = v0
    u_last[-1] = vh
    RHS = sp.dot(Mdt, u_last) + sp.dot(M, Q) + cc
    RHS[0] = v0
    RHS[-1] = vh
    u = sp.linalg.solve(LHS, RHS)

    if t > 3 and sp.all(abs(u - u_last) < 10e-6):
        break
    u_last = sp.copy(u)


    # with open('turbulent.csv', mode='a') as turbulent:
    #     writer = csv.writer(turbulent, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     writer.writerow([viscosity_turb, u_plus, u_plus_a, y_plus])


# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------

with open('flowResultT.csv', mode='a') as flowResultT:
    writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(u_last)

with open('properties.csv', mode='w') as properties:
    writer = csv.writer(properties, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow([x, dt, tempo, viscosity_din, rho_fld, viscosity_kin, h, L])


plt.figure(1)
plt.plot(x, u)

# plt.figure(4)
# plt.plot(lc, x)

plt.figure(5)
plt.plot(x, du_dy)

plt.figure(3)
plt.semilogx(y_plus[0:nodes/2], u_plus[0:nodes/2])
plt.semilogx(y_plus[0:nodes/2], y_plus[0:nodes/2])
plt.semilogx(y_plus[0:nodes/2], u_plus_a[0:nodes/2], 'k-')
plt.xlim(1, 10**3)

plt.figure(2)
plt.plot(x, viscosity_turb)
plt.show()