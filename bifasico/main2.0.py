# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import semiLagrangean as sl
import libmalha as lm
import GMesh as gm
import os
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cwd = os.getcwd()

# Definição da malha
arquivo = "circulo"

print("Reading .msh file")
malha = gm.GMesh(arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN

# Parametros
nodes = len(x)
elenum = len(ien)
q = 100
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

# Geração de Particulas


n_particulas = 10
particulas = sp.zeros((n_particulas, 3))
particulas_vel = sp.zeros((n_particulas, 3))

for i in range(n_particulas):
    particulas[i, 0] = random.randint(-70,70)*0.01
    particulas[i, 1] = random.randint(-70,70)*0.01
    particulas[i, 2] = random.randint(1, 140)*0.01

# Montagem de matrizes

print("Assembling Matrices")


K, M = lm.fem_matrix(x, y, elenum, nodes, nu_total, ien)

MQ = sp.dot(M, Q)
Mdt = M/dt
KM = K + Mdt


    # Implementando condições de contorno de DIRICHLET e INICIAL

print("Defining Boundary and Initial conditions")

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

    # Analitico

    Tx_a = sp.zeros(nodes)
    Ty_a = sp.zeros(nodes)

    for i in range(nodes):
        Ty_a[i] = -15. / 48. + 24 * malha.X[i] + 1.5 * malha.Y[i] ** 2 - malha.Y[i] ** 4

    # Solucao do sistema


Tfim = sp.zeros((nodes, tempo))
Tfim[0:nodes, 0] = sp.copy(Ti)

def solve():
    for i in range(1, tempo):
        B = MQ + sp.dot(Mdt, Tfim[0:nodes, i-1]) + cc
        for j in range(quantcc):
            index = int(malha.dirichlet_points[j][0]) - 1
            B[index] = malha.dirichlet_points[j][1]

        Tfim[0:nodes, i] = linalg.solve(KK, B)
        print( "Solving System " + str((float(i)/tempo)*100) + "%")

        vtk = io.InOut(x, y, ien, len(x), len(ien), Tfim[:,i-1])
        vtk.saveVTK(cwd + "/results", arquivo+str(i-1))

    vtk = io.InOut(x, y, ien, len(x), len(ien), Tfim[:,tempo-1])
    vtk.saveVTK(cwd + "/results", arquivo+str(tempo-1))
    return

#solve()


# TROCAR Z POR Y
def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = sp.linspace(0, height_z, 50)
    theta = sp.linspace(0, 2*sp.pi, 50)
    theta_grid, z_grid=sp.meshgrid(theta, z)
    x_grid = radius*sp.cos(theta_grid) + center_x
    y_grid = radius*sp.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid


def plot(_pos, _name):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    Xc,Yc,Zc = data_for_cylinder_along_z(0.0,0.0, 1, 20)
    ax.plot_surface(Xc, Zc, Yc, alpha=0.2)
    ax.view_init(elev=26, azim=30)
    ax.set_axis_off()
    for i in range(n_particulas):
        ax.scatter(_pos[i,0], _pos[i,2], _pos[i,1])
    plt.show()
    fig.savefig(_name, format='pdf')
