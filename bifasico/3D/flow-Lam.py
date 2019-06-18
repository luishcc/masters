# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import sys
import time
import csv
import matrix as ma
import matplotlib.pyplot as plt

cwd = os.getcwd()

# Definição da malha
arquivo = "circulo"

print("Reading .msh file")
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros

# --------------------------------------------------
#   Problem Parameters

dt = 100
dt_inv = 1. / dt
tempo = 10000

nodes = len(x)
elenum = len(ien)
q = 100
Q = q * sp.ones(nodes)


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

x = sp.linspace(0, h, 1000)

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


def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    gx_local = sp.zeros((3, 3), dtype="float64")
    gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        Area = (a[0] + a[1] + a[2]) / 2.

        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)
                gx_local[i, j] = b[j] * (1/6.)
                gy_local[i, j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local] * (Area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global, gx_global, gy_global

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
rhs_part = sp.dot(M, Q) + cc


# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)


# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------

with open('flow-L.csv', mode='w') as flowResultL:
    writer = csv.writer(flowResultL, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(x)

for t in range(tempo):
    print t, " / ", tempo
    with open('flow-L.csv', mode='a') as flowResultL:
        writer = csv.writer(flowResultL, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(u_last)

    RHS = sp.dot(Mdt, u_last) + rhs_part
    RHS[0] = v0
    RHS[-1] = vh
    u = sp.linalg.solve(LHS, RHS)


    norm = u - u_last
    norm = sp.dot(norm, norm)
    if norm < 10e-8:
        break

    u_last = sp.copy(u)


# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------

with open('flow-L.csv', mode='a') as flowResult:
    writer = csv.writer(flowResult, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(u_last)

with open('prop-L.csv', mode='w') as properties:
    writer = csv.writer(properties, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow([dt, tempo, viscosity_din, rho_fld, viscosity_kin, h, L])

#
#
plt.figure(1)
plt.plot(x, u)
plt.show()