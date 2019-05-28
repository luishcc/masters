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


def get_vel(_pos, _fluid_vel, _x, _y, _rp):
    fluid_vel_part = sp.zeros(5)
    f = sp.interpolate.bisplrep(_x, _y, _fluid_vel)
    fluid_vel_part[0] = sp.interpolate.bisplev(_pos[0]-_rp, _pos[1], f)       # vx esquerda
    fluid_vel_part[1] = sp.interpolate.bisplev(_pos[0]+_rp, _pos[1], f)       # vx direita
    fluid_vel_part[2] = sp.interpolate.bisplev(_pos[0], _pos[1]+_rp, f)       # vy cima
    fluid_vel_part[3] = sp.interpolate.bisplev(_pos[0], _pos[1]-_rp, f)       # vy baixo
    fluid_vel_part_center = sp.interpolate.bisplev(_pos[0], _pos[1], f)     # v centro
    return fluid_vel_part, fluid_vel_part_center


# TROCAR Z POR Y
def data_for_cylinder_along_z(center_x, center_y, radius, height_z):
    z = sp.linspace(0, height_z, 50)
    theta = sp.linspace(0, 2 * sp.pi, 50)
    theta_grid, z_grid = sp.meshgrid(theta, z)
    x_grid = radius * sp.cos(theta_grid) + center_x
    y_grid = radius * sp.sin(theta_grid) + center_y
    return x_grid, y_grid, z_grid


def plot(_pos, _name):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    Xc, Yc, Zc = data_for_cylinder_along_z(0.0, 0.0, 1, 20)
    ax.plot_surface(Xc, Zc, Yc, alpha=0.2)
    ax.view_init(elev=26, azim=30)
    ax.set_axis_off()
    for i in range(n_particulas):
        ax.scatter(_pos[i, 0], _pos[i, 2], _pos[i, 1])
    plt.show()
    fig.savefig(_name, format='pdf')


def gravity(_rho, rp):
    f_g = (4 / 3.0) * sp.pi * (rp ** 3) * 9.81 * _rho
    return f_g


def drag(_pos, fluid_vel_part, ):
    return