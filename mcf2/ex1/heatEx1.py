# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import os
import libmalha as lm
import matplotlib.pyplot as plt


cwd = os.getcwd()


# -------------------------------------------------------
# Exercicio
#
# div(grad(Phi)) = 0 , com dirichlet em toda fronteira 2D


# -------------------------------------------------------
#   Generating 2D mesh
# -------------------------------------------------------

Lx = 1.0
Ly = 1.0

nele_x = 3
nele_y = 3

x, y, ien = lm.lin2d(nele_x, Lx, nele_y, Ly)

ien = lm.ret_tri(ien)

nodes = len(x)
elements = len(ien)


# Boundary

node_dirichlet = []

for i in range(nodes):
    if x[i] == 0:
        f = y[i]
        node_dirichlet.append([i, f])
        continue
    if y[i] == 0:
        f = x[i]
        node_dirichlet.append([i, f])
        continue
    if x[i] == 1:
        f = y[i]**2+1
        node_dirichlet.append([i, f])
        continue
    if y[i] == 1:
        f = x[i]**2+1
        node_dirichlet.append([i, f])
        continue

# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")

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

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]


    return  k_global
# --------------------------------------
# Total Matrices (Heat equation)

K = fem_matrix(x, y, elements, nodes, ien)


# ---------------------------------------
# Condições de contorno
#  ---------------------------------------

dirichlet_len_total = len(node_dirichlet)

RHS_BC = sp.zeros(nodes)
LHS = sp.copy(K)

for i in range(dirichlet_len_total):
    index = int(node_dirichlet[i][0])
    value = node_dirichlet[i][1]
    for j in range(nodes):
        RHS_BC[j] -= value * K[j, index]
        if j != index:
            LHS[index, j] = 0
            LHS[j, index] = 0
        else:
            LHS[index, j] = 1

RHS = RHS_BC
for i in range(dirichlet_len_total):
    index = int(node_dirichlet[i][0])
    RHS[index] = node_dirichlet[i][1]

Phi = sp.linalg.solve(LHS, RHS)

# Salvar VTK
vtk_t = Io.InOut(x, y, ien, nodes, elements, Phi, None, None, None, None, None, None)
vtk_t.saveVTK(cwd+"/results", "Ex1")
