# -*- coding: utf-8 -*-
import scipy as sp
import os
import random
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

cwd = os.getcwd()

# --------------------------------------------------
#   Plot functions

def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = sp.linspace(0, height_z, 50)
    theta = sp.linspace(0, 2*sp.pi, 50)
    theta_grid, z_grid=sp.meshgrid(theta, z)
    x_grid = radius*sp.cos(theta_grid) + center_x
    y_grid = radius*sp.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def plot(_x, _y, _z, _name):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    Xc,Yc,Zc = data_for_cylinder_along_z(0.0, 0.0, 0.5, L)
    ax.plot_surface(Xc, Zc, Yc, alpha=0.2)
    ax.view_init(elev=26, azim=30)
    # ax.set_axis_off()
    ax.set_ylim(-0.5, L+0.5)
    for i in range(len(_x)):
        ax.scatter(_x, _z, _y)
    # plt.show()
    fig.savefig(_name, format='png')

# --------------------------------------------------
#   Problem Parameters

dt = 0.05
dt_inv = 1. / dt
tempo = 1000
g = 9.81
h = 1.
L = 5*h

collision_coef = 0.3

# fluid
rho_fld = 1000
viscosity_din = 0.8e-03
viscosity_kin = viscosity_din / rho_fld

# particle
rho_part = 1000.0
radius = 0.0005
diameter = 2 * radius
volume = (4/3.) * sp.pi * radius**3

constant = dt / ((rho_part + 0.5 * rho_fld) * volume)

# --------------------------------------------------
#   particle Generation

n_particle = 5

posx_particle = sp.zeros(n_particle)
posy_particle = sp.zeros(n_particle)
posz_particle = sp.zeros(n_particle)

for i in range(n_particle):
    posx_particle[i] = random.randint(-35, 35)*0.01
    posy_particle[i] = random.randint(-35, 35)*0.01
    posz_particle[i] = random.randint(1, 10)*0.01


# particle velocity
vx = sp.zeros(n_particle)
vy = sp.zeros(n_particle)
vz = sp.zeros(n_particle)
vx_last = sp.zeros(n_particle)
vy_last = sp.zeros(n_particle)
vz_last = sp.zeros(n_particle)

# particle forces
f_g = sp.ones(n_particle) * (rho_part - rho_fld) * volume * g       # constant gravity

f_dragx = sp.zeros(n_particle)
f_dragy = sp.zeros(n_particle)
f_dragz = sp.zeros(n_particle)

f_liftx = sp.zeros(n_particle)
f_lifty = sp.zeros(n_particle)


# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------



for t in range(tempo):

    for i in range(0, n_particle):

       # get fluid velocity at particle position
        u = -2*(posx_particle[i]**2) - 2*(posy_particle[i]**2) + 2

        # relative velocity between fluid and particle
        relative_velx = -1 * vx[i]
        relative_vely = -1 - vy[i]
        relative_velz = u - vz[i]
        abs_relative_vel = sp.sqrt(relative_velx**2 + relative_vely**2 + relative_velz**2)

        # Compute forces

        f_dragx[i] = 3 * sp.pi * viscosity_din * diameter * relative_velx
        f_dragy[i] = 3 * sp.pi * viscosity_din * diameter * relative_vely
        f_dragz[i] = 3 * sp.pi * viscosity_din * diameter * relative_velz

        du_dy = -4 * posy_particle[i]
        du_dx = -4 * posx_particle[i]
        f_liftx[i] = 1.61 * sp.sqrt(viscosity_din * rho_fld * abs(du_dx)) * (diameter**2) * abs_relative_vel * du_dx
        f_lifty[i] = 1.61 * sp.sqrt(viscosity_din * rho_fld * abs(du_dy)) * (diameter**2) * abs_relative_vel * du_dy

        # Update particle velocity

        vy[i] = vy_last[i] + (f_lifty[i] - f_g[i] + f_dragy[i]) * constant
        vx[i] = vx_last[i] + (f_liftx[i] + f_dragx[i]) * constant
        vz[i] = vz_last[i] + f_dragz[i] * constant

        print f_g[i], f_lifty[i], f_liftx[i], f_dragy[i], f_dragx[i], f_dragz[i]

        # Update particle position
        posx_particle[i] = vx[i] * dt + posx_particle[i]
        posy_particle[i] = vy[i] * dt + posy_particle[i]
        posz_particle[i] = vz[i] * dt + posz_particle[i]

        # Check if particle is inside limits


        if posz_particle[i] > L:
            posz_particle[i] = posz_particle[i] - L

        vx_last[i] = vx[i]
        vy_last[i] = vy[i]
        vz_last[i] = vz[i]

    plot(posx_particle, posy_particle, posz_particle, 'result/test-'+str(t+1)+'.png')
    time.sleep(1.5)

# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------
