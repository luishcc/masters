# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import random
from scipy import interpolate
import matplotlib.pyplot as plt
import sys

cwd = os.getcwd()
# --------------------------------------------------
# Plot function
def defPlot(_xp, _yp, _n, _u, _yf, L, t):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    # ax2 = fig.add_subplot(1,2,2)
    for i in range(_n):
        ax.scatter(_xp[i], _yp[i])
    ax.plot(_u, _yf, 'b')
    ax.plot([0, L], [0, 0], 'k-')
    ax.plot([0, L], [1, 1], 'k-')
    ax.plot([0, L], [0.5, 0.5], 'r--')
    ax.plot([0, 0], [0, 1], 'k--')
    ax.plot([L, L], [0, 1], 'k--')
    ax.set_xlim(-1, L + 1)
    ax.set_ylim(-0.1, 1.1)
    plt.savefig('result/test' + str(t) + '.jpg', format='jpg')
    return


# --------------------------------------------------
#   Flow Input

flow = 'T'      # L - laminar; T - turbulent
u = sp.loadtxt('flow-' + flow + '.csv', delimiter=',')
u = u*5

fluid_properties = sp.loadtxt('prop-' + flow + '.csv', delimiter=',')
# Order ->  0 dt, 1 tempo, 2 viscosity_din, 3 rho_fld, 4 viscosity_kin, 5 h, 6 L

# --------------------------------------------------
#   Problem Parameters

dt = 0.01
dt_inv = 1. / dt
tempo = 1000
g = 9.81
h = fluid_properties[5]
L = fluid_properties[6]

collision_coef = 0.3

# fluid
rho_fld = fluid_properties[3]
viscosity_din = fluid_properties[2]
viscosity_kin = fluid_properties[4]

# particle
rho_part = 1300.0
radius = 0.0005
diameter = 2 * radius
volume = (4/3.) * sp.pi * radius**3

constant = dt / ((rho_part + 0.5 * rho_fld) * volume)

# --------------------------------------------------
#   particle Generation

n_particle = 5

posx_particle = sp.zeros(n_particle)
posy_particle = sp.zeros(n_particle)

for i in range(n_particle):
    posx_particle[i] = random.randint(1, 10)*0.1
    posy_particle[i] = random.randint(1, 99)*0.01
# posy_particle[0] = 0.3


fld_accx = sp.zeros(n_particle)


# particle velocity
vx = sp.zeros(n_particle)
vy = sp.zeros(n_particle)
vx_last = sp.zeros(n_particle)
vy_last = sp.zeros(n_particle)

# particle forces
f_g = sp.ones(n_particle) * (rho_part - rho_fld) * volume * g       # constant gravity
f_dragx = sp.zeros(n_particle)
f_dragy = sp.zeros(n_particle)
f_mass_partialx = sp.zeros(n_particle)                  # f_mass_partialy = 0 because fluid velocity uy = 0
f_lift = sp.zeros(n_particle)                           # Only act on y direction


# --------------------------------------------------
# --------------------------------------------------

# ---------   Beginning of Time Loop  --------------

#defPlot(posx_particle, posy_particle, n_particle, u, x, L, 0)      # Plot initial condition


u_last_interp = sp.interpolate.interp1d(u[0], u[1], fill_value=0, bounds_error=False)

u_interp_aux = sp.interpolate.interp1d(u[0], u[-1], fill_value=0, bounds_error=False)
u_last_interp_aux = sp.interpolate.interp1d(u[0], u[-2], fill_value=0, bounds_error=False)

for t in range(tempo):

    if t+3 > len(u):
        u_interp = u_interp_aux
        u_last_interp = u_last_interp_aux

    # Interpolate fluid velocity at n+1 and n time step
    u_interp = sp.interpolate.interp1d(u[0], u[t+2], fill_value=0, bounds_error=False)

    for i in range(0, n_particle):
        # get fluid velocity at particle position
        u_n0 = u_last_interp(posy_particle[i])
        u_n1 = u_interp(posy_particle[i])
        fld_accx[i] = (u_n1 - u_n0) * dt

        # relative velocity between fluid and particle
        relative_velx = u_n0 - vx[i]
        relative_vely = -1 * vy[i]
        abs_relative_vel = sp.sqrt(relative_velx**2 + relative_vely**2)

        # Compute forces
        if abs_relative_vel != 0 :
            Re_r = abs_relative_vel * diameter * (1./viscosity_kin)
            drag_coef = 27. / Re_r
            drag_factor = drag_coef * Re_r * (1/24.)
        else:
            drag_factor = 0
        f_dragx[i] = 3 * sp.pi * viscosity_din * diameter * relative_velx * drag_factor
        f_dragy[i] = 3 * sp.pi * viscosity_din * diameter * relative_vely * drag_factor

        f_mass_partialx[i] = 0.5 * rho_fld * volume * fld_accx[i]

        du_dy = (u_last_interp(posy_particle[i]+radius) - u_last_interp(posy_particle[i]-radius)) / diameter
        f_lift[i] = 1.61 * sp.sqrt(viscosity_din * rho_fld * abs(du_dy)) * (diameter**2) * abs_relative_vel * du_dy

        # Update particle velocity
        vy_f = sp.random.normal() * sp.sqrt(1.5) * 0.041 * u_n1
        vx_f = sp.random.normal() * sp.sqrt(1.5) * 0.041 * u_n1

        vy[i] = vy_last[i] + (f_lift[i] - f_g[i] + f_dragy[i]) * constant + vy_f
        vx[i] = vx_last[i] + (f_mass_partialx[i] + f_dragx[i]) * constant + vx_f


        print f_g[i],  f_lift[i], f_dragy[i], f_dragx[i], f_mass_partialx[i]

        # Update particle position
        posx_particle[i] = vx[i] * dt + posx_particle[i]
        posy_particle[i] = vy[i] * dt + posy_particle[i]

        # Check if particle is inside limits
        if posy_particle[i] < 0:
            aux = posy_particle[i] / vy[i]
            vy[i] = -1 * vy[i] * collision_coef
            posy_particle[i] = aux * vy[i]

        if posy_particle[i] > h:
            aux = (posy_particle[i]-h) / vy[i]
            vy[i] = -1 * vy[i] * collision_coef
            posy_particle[i] = h + aux * vy[i]


        if posx_particle[i] > L:
            posx_particle[i] = posx_particle[i] - L

        vx_last[i] = vx[i]
        vy_last[i] = vy[i]

    u_last_interp = u_interp

    defPlot(posx_particle, posy_particle, n_particle, u[t+1], u[0], L, t+1)

# ---------   End of Time Loop   -------------------

# --------------------------------------------------
# --------------------------------------------------


