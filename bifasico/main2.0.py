# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import GMesh as gm
import os

cwd = os.getcwd()

####### DEV

# Escrito por Luís Cunha
# Última atualização : 09/01/2018

# Código para a solução da EDP d²T               dT
#                              --- * alfa + Q = --- + v . grad(T)
#                              dx²               dt
# através de elementos Finitos.


# Definição da malha
arquivo = "valid"

print "Reading .msh file"
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros
nodes = len(x)
elenum = len(ien)
q = 0
Q = q * sp.ones(nodes)
dt = 0.1
tempo = 200
alfax = 1
alfay = 1

velx = sp.zeros(nodes)
vely = sp.zeros(nodes)

for i in range(nodes):
    velx[i] = 0.125 - y[i]**2/2.0
    vely[i] = 0.0

# Montagem de matrizes

print "Assembling Matrices"

def fem_matrix(_x, _y, _numele, _numnode, _ax, _ay, _ien, _vx, _vy):
    BB_local = sp.zeros((3, 3))
    NN_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    NVB_local = sp.zeros((3, 3))
    vx = sp.zeros(3)
    vy = sp.zeros(3)
    a = sp.zeros(3)
    b = sp.zeros(3)
    c = sp.zeros(3)
    yy = sp.zeros(3)
    xx = sp.zeros(3)
    BB_global = sp.zeros((_numnode, _numnode))
    NN_global = sp.zeros((_numnode, _numnode))
    NVB_global = sp.zeros((_numnode, _numnode))

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]
            vx[i] = _vx[ien[elem, i]]
            vy[i] = _vy[ien[elem, i]]

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]
        Area = (a[0] + a[1] + a[2]) / 2.

        for i in range(3):
            for j in range(3):
                BB_local[i, j] = (_ax * b[i] * b[j] + _ay * c[i] * c[j]) / (4 * Area)
                NVB_local[i,j] = (vx[i] * b[j] + vy[i] * c[j]) * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                BB_global[i_global, j_global] += BB_local[i_local, j_local]
                NN_global[i_global, j_global] += NN_local[i_local, j_local]* (Area/12.)
                NVB_global[i_global, j_global] += NVB_local[i_local, j_local]


    return  BB_global, NN_global, NVB_global

K, M, VB = fem_matrix(x, y, elenum, nodes, alfax, alfay, ien, velx, vely)

MQ = sp.dot(M, Q)
Mdt = M/dt
KM = K + VB + Mdt


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

    # Implementando condições de contorno de Neumann

len_neu = len(malha.Boundary_Neumann)
for i in range(len_neu):
    node1 = malha.Boundary_Neumann[i][0]
    node2 = malha.Boundary_Neumann[i][1]
    for j in range(len(malha.neumann_points)):
        if node1 == int(malha.neumann_points[j,0]):
            index1 = j
        if node2 == int(malha.neumann_points[j,0]):
            index2 = j
    value1 = malha.neumann_points[index1][1]
    value2 = malha.neumann_points[index2][1]
    length = sp.sqrt((x[node1]-x[node2])**2 + (y[node1]-y[node2])**2)
    neu_bc1 = (length/2)*value1
    neu_bc2 = (length/2)*value2
    cc[node1-1] -= neu_bc1
    cc[node2-1] -= neu_bc2


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