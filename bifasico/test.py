# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as io
import semiLagrangean as sl
import libmalha as lm
import GMesh as gm
import os
from scipy import interpolate
import random
import forces as ff
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cwd = os.getcwd()

# Definição da malha
arquivo = "circulo"

print("Reading .msh file")
malha = gm.GMesh(arquivo + ".msh")


# Parametros
nodes = len(malha.X)
elenum = len(malha.IEN)

# Geração de Particulas


n_particulas = 10
particule_pos = sp.zeros((n_particulas, 3))
particule_vel = sp.zeros((n_particulas, 3))

for i in range(n_particulas):
    particule_pos[i, 0] = random.randint(-70, 70) * 0.01
    particule_pos[i, 1] = random.randint(-70, 70) * 0.01
    particule_pos[i, 2] = random.randint(1, 140) * 0.01


# Analitico


Ty_a = sp.sin(malha.X**2+malha.Y**2)

fluid_velocity_corner, fluid_velocity = ff.get_vel(particule_pos[0], Ty_a, malha.X, malha.Y, 0.02)

fluid_velocity_last = fluid_velocity - fluid_velocity * 0.1

print fluid_velocity_corner
f_lift = ff.lift(fluid_velocity, fluid_velocity_corner[0:2], 0, 2, 0.02, 2)
f_mass = ff.virtualMass(fluid_velocity, fluid_velocity_last, 2, 0.02, 0.05)


print f_lift, f_mass
