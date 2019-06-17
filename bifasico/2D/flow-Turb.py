# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
from scipy.sparse.linalg import cg
import sys
import csv
import matrix as ma
import matplotlib.pyplot as plt

cwd = os.getcwd()

# --------------------------------------------------
#   Problem Parameters

dt = 100
tempo = 100000

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -0.012

# boundary
v0 = 0
vh = 0

h = 1.0
L = 8*h     # Duct length

fine = 1000
coarse = 500
d_fine = 0.15*h
d_coarse = 0.85*h

y_fine = sp.linspace(0, d_fine, fine, endpoint=False)
y_coarse = sp.linspace(d_fine, d_coarse, coarse, endpoint=False)
y_fine2 = sp.linspace(d_coarse, h, fine+1)
y = sp.append(y_fine, y_coarse)
y = sp.append(y, y_fine2)

# y = sp.linspace(0, h, 10000)

nodes = len(y)
elem = nodes-1

ien = sp.zeros((elem, 2), dtype='int32')
dy = sp.zeros(elem)
for i in range(0, elem):
    dy[i] = abs(y[i+1] - y[i])
    for j in [0, 1]:
        ien[i, j] = int(i + j)

k = sp.array([[1, -1], [-1, 1]])
m = sp.array([[2, 1], [1, 2]])
grad = sp.array([[-0.5, 0.5], [-0.5, 0.5]])
Q = sp.ones(nodes) * (-grad_p / rho_fld)
M = sp.zeros((nodes, nodes))
G = sp.zeros((nodes, nodes))


print "assembly"
for e in range(elem):
    for i_local in range(2):
        i_global = ien[e, i_local]
        for j_local in range(2):
            j_global = ien[e, j_local]
            M[i_global, j_global] += m[i_local, j_local] * (dy[e]) / 6.
            G[i_global, j_global] += grad[i_local, j_local]

Mdt = M / dt



# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)
du_dy = sp.zeros(nodes)
viscosity_turb = sp.zeros(nodes)

y_plus = sp.zeros(nodes)
u_plus = sp.zeros(nodes)
u_plus_a = sp.zeros(nodes)
lc = sp.zeros(nodes)


M_inv = sp.zeros(nodes)

for i in range(nodes):
    for j in range(nodes):
        M_inv[i] += M[i, j]
    if 0 <= y[i] <= 0.1 * h:
        lc[i] = 0.41 * y[i]
    elif h >= y[i] >= 0.9 * h:
        lc[i] = 0.41 * (h - y[i])
    else:
        lc[i] = 0.41 * 0.1 * h

    M_inv[i] = 1./M_inv[i]

print "M inv lumped"




# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------


A_const = 26
for t in range(tempo):
    print t, " / ", tempo
    # with open('flowResultT.csv', mode='a') as flowResultT:
    #     writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     writer.writerow(u_last)

    u[0] = 0
    u[-1] = 0

    t_wall = viscosity_din * (0.5/dy[0]) * (-3 * u[0] + 4 * u[1] - u[2])
    u_t = sp.sqrt(t_wall / rho_fld)
    y_plus = (u_t / viscosity_kin) * y
    u_plus = (1./u_t) * u
    print t_wall
    du_dy = sp.dot(G, u)
    for i in range(nodes):
        if 0 <= y[i] <= 0.1 * h:
            D = 1 - sp.exp(-y_plus[i]/A_const)
            lc[i] = 0.41 * y[i] * D
        elif h >= y[i] >= 0.9 * h:
            index = nodes - i - 1
            D = 1 - sp.exp(-y_plus[index]/A_const)
            lc[i] = 0.41 * (h - y[i]) * D

        u_plus_a[i] = (1. / 0.41) * sp.log(y_plus[i]) + 5.5
        viscosity_turb[i] = (lc[i]**2) * abs(du_dy[i])

    # viscosity_turb = sp.multiply(M_inv, viscosity_turb)  # TESTAR


    K = sp.zeros((nodes, nodes))
    for e in range(elem):
        v1 = ien[e, 0]
        v2 = ien[e, 1]
        visc = (viscosity_turb[v1] + viscosity_turb[v2]) * 0.5
        for i_local in range(2):
            i_global = ien[e, i_local]
            for j_local in range(2):
                j_global = ien[e, j_local]
                K[i_global, j_global] += k[i_local, j_local] * (viscosity_kin + visc) * (1. / dy[e])

    LHS = Mdt + K
    cc = sp.zeros(nodes)
    LHS_copy = sp.copy(LHS)
    for i in range(nodes):
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
    writer.writerow([y, dt, tempo, viscosity_din, rho_fld, viscosity_kin, h, L])


plt.figure(1)
plt.plot(y, u)

plt.figure(4)
plt.plot(lc, y)

plt.figure(5)
plt.plot(y, du_dy)

x = sp.linspace(0.1, 20, 10000)
plt.figure(3)
plt.semilogx(x, x, label='u=y')
plt.semilogx(y_plus[0:nodes/2], u_plus_a[0:nodes/2], 'k--', label='log law')
plt.semilogx(y_plus[0:nodes/2], u_plus[0:nodes/2], 'y.', label='numeric')
plt.xlim(0.1, 10**3)
plt.ylim(0, 25)
plt.legend()

plt.figure(2)
plt.plot(y, viscosity_turb)
plt.show()