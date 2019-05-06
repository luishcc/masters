# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import semiLagrangean as sl
import libmalha as lm
import GMesh as gm
import os
import random

cwd = os.getcwd()

# Definição da malha
arquivo = "circulo"

print "Reading .msh file"
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros
nodes = len(x)
elenum = len(ien)
q = 1
Q = q * sp.ones(nodes)

dt = 0.1
tempo = 200

visc = 2
nu_total = sp.ones(nodes) * visc
nu_t = sp.zeros(nodes)

vel_a = sp.zeros(nodes)

for i in range(nodes):
    vel_a[i] = 0.125 - 0.5 * y[i]**2
    nu_t[i] = 0.5 * y[i]**2 
    
nu_total += nu_t

# Montagem de matrizes

print "Assembling Matrices"

def fem_matrix(_x, _y, _numele, _numnode, _ax, _ien):
    k_local = sp.zeros((3, 3))
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    a = sp.zeros(3)
    b = sp.zeros(3)
    c = sp.zeros(3)
    yy = sp.zeros(3)
    xx = sp.zeros(3)
    K_global = sp.zeros((_numnode, _numnode))
    M_global = sp.zeros((_numnode, _numnode))
    

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]
        Area = (a[0] + a[1] + a[2]) * 0.5

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (_ax * b[i] * b[j] + _ay * c[i] * c[j]) / (4 * Area)
                
        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                K_global[i_global, j_global] += k_local[i_local, j_local] * _ax
                M_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
      
    return  K_global, M_global

K, M, VB = fem_matrix(x, y, elenum, nodes, nu_total, ien)

MQ = sp.dot(M, Q)
Mdt = M/dt
KM = K + Mdt


    # Implementando condições de contorno de DIRICHLET e INICIAL

print "Defining Boundary and Initial conditions"

quantcc = len(malha.dirichlet_points)
cc = sp.zeros(nodes)
KK = sp.copy(KM)
Ti = sp.zeros(nodes)

for i in range(quantcc):
    index = int(malha.dirichlet_points[i][0]-1)
    value = malha.dirichlet_points[i][1]
    Ti[index] = malha.dirichlet_points[i, 1]
    for j in range(nodes):
        cc[j] -= value * KM[j, index]
        if j != index:
            KK[index, j] = 0
            KK[j, index] = 0
        else:
            KK[index, j] = 1


    # Solucao do sistema

Tfim = sp.zeros((nodes, tempo))
Tfim[0:nodes, 0] = sp.copy(Ti)

for i in range(1, tempo):
    B = MQ + sp.dot(Mdt, Tfim[0:nodes, i-1]) + cc
    for j in range(quantcc):
        index = int(malha.dirichlet_points[j][0]) - 1
        B[index] = malha.dirichlet_points[j][1]

    Tfim[0:nodes, i] = linalg.solve(KK, B)
    print "Solving System " + str((float(i)/tempo)*100) + "%"

mint = Tfim.min()
maxt = Tfim.max()

# Analitico

Tx_a = sp.zeros(nodes)
Ty_a = sp.zeros(nodes)

for i in range(nodes):
    Ty_a[i] = -15./48. + 24*malha.X[i] + 1.5*malha.Y[i]**2-malha.Y[i]**4

print "Generating .VTK file(s)"
for i in range(tempo):
    vtk = io.InOut(x, y, ien, len(x), len(ien), Tfim[:,i], Ty_a, Tx_a)
    vtk.saveVTK(cwd + "/results", arquivo+str(i))
print "Done"
