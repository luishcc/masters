# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import sys
import csv
import matrix as ma
import matplotlib.pyplot as plt

cwd = os.getcwd()


# --------------------------------------------------
#   Problem Parameters

dt = 0.1
dt_inv = 1. / dt
tempo = 10000

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -12.

# boundary
v0 = 0
vh = 0


# --------------------------------------------------
#   Fluid Mesh and Matrices Generation

h = 1   # Duct height
L = 5*h     # Duct length
fine = 20
coarse = 10
d_fine = 0.15*h
d_coarse = 0.8*h

x_fine = sp.linspace(0, d_fine, fine, endpoint=False)
x_coarse = sp.linspace(d_fine, d_coarse, coarse, endpoint=False)
x_fine2 = sp.linspace(d_coarse, h, fine+1)
x = sp.append(x_fine, x_coarse)
x = sp.append(x, x_fine2)

# x = sp.linspace(0, h, 10)

nodes = len(x)
elem = nodes-1

ien = sp.zeros((elem, 2), dtype='int32')
dx = sp.zeros(elem)
for i in range(0, elem):
    dx[i] = abs(x[i+1] - x[i])
    for j in [0, 1]:
        ien[i, j] = int(i + j)


# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)


k = sp.array([[1, -1], [-1, 1]])
m = sp.array([[2, 1], [1, 2]])
b = sp.array([1, 1])
B = sp.zeros(nodes)
K = sp.zeros((nodes, nodes))
M = sp.zeros((nodes, nodes))

for e in range(elem):
    for i_local in range(2):
        i_global = ien[e, i_local]
        for j_local in range(2):
            j_global = ien[e, j_local]
            K[i_global, j_global] += k[i_local, j_local] * (1. / dx[e])
            M[i_global, j_global] += m[i_local, j_local] * (dx[e]) / 6.

Q = sp.ones(nodes) * (-grad_p / rho_fld)

Mdt = M * dt_inv
LHS = Mdt + viscosity_din * K

v0 = 0
vh = 0


cc = sp.zeros(nodes)

LHS_copy = sp.copy(LHS)
for i in range(nodes):
    cc[i] = cc[i] - v0 * LHS_copy[i, 0]
    cc[i] = cc[i] - vh * LHS_copy[i, -1]
    LHS[0, i] = 0
    LHS[i, 0] = 0
    LHS[-1, i] = 0
    LHS[i, -1] = 0


LHS[0, 0] = 1
LHS[-1, -1] = 1

# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)


# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------

with open('flowResultL.csv', mode='w') as flowResultL:
    writer = csv.writer(flowResultL, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)


for t in range(tempo):
    print t, " / ", tempo
    with open('flowResultL.csv', mode='a') as flowResultL:
        writer = csv.writer(flowResultL, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(u_last)

    RHS = sp.dot(Mdt, u_last) + sp.dot(M, Q) + cc
    RHS[0] = v0
    RHS[-1] = vh
    u = sp.linalg.solve(LHS, RHS)

    u_last = sp.copy(u)

# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------

with open('flowResult.csv', mode='a') as flowResult:
    writer = csv.writer(flowResult, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(u_last)

with open('properties.csv', mode='w') as properties:
    writer = csv.writer(properties, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow([x, dt, tempo, viscosity_din, rho_fld, viscosity_kin, h, L])
#
#
# plt.figure(1)
# plt.plot(x, u)
# plt.show()