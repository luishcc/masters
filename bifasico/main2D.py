# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import random
from scipy import interpolate
import matplotlib.pyplot as plt
import sys
import shutil

cwd = os.getcwd()
if os.path.isdir(cwd+'/result2D') == True:
    shutil.rmtree(cwd+'/result2D')

os.mkdir(cwd+'/result2D')


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
#   Problem Parameters

dt = 0.05
dt_inv = 1 / dt
tempo = 300
g = 9.81

collision_coef = 0.4

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -12.

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

# --------------------------------------------------
#   Fluid Mesh and Matrices Generation

h = 1   # Duct height
L = 5*h     # Duct length
fine = 20
coarse = 10
d_fine = 0.2*h
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




k = sp.array([[1, -1], [-1, 1]])
m = sp.array([[2, 1], [1, 2]])
b = sp.array([1, 1])
B = sp.zeros(nodes)
K = sp.zeros((nodes, nodes))
M = sp.zeros((nodes, nodes))

for elem in range(elem):
    for i_local in range(2):
        i_global = ien[elem, i_local]
        B[i_global] += b[i_local] * ((-grad_p / rho_fld) * dx[elem] * 0.5)
        for j_local in range(2):
            j_global = ien[elem, j_local]
            K[i_global, j_global] += k[i_local, j_local] / dx[elem]
            M[i_global, j_global] += m[i_local, j_local] * (dx[elem])/6


# --------------------------------------------------
#   Fluid Dirichlet Boundary Condition and System Definition

Mdt = M * dt_inv
LHS = Mdt + viscosity_kin * K

v0 = 0
vh = 0

LHS_copy = sp.copy(LHS)
for i in range(nodes):
    B[i] = B[i] - v0 * LHS_copy[i, 0]
    B[i] = B[i] - vh * LHS_copy[i, -1]
    LHS[0, i] = 0
    LHS[i, 0] = 0
    LHS[-1, i] = 0
    LHS[i, -1] = 0

B[0] = v0
B[-1] = vh
LHS[0, 0] = 1
LHS[-1, -1] = 1


# --------------------------------------------------
#   Allocating Fluid and Particle Variables

# fluid velocity (and acceleration in x direction)
u_last = sp.zeros(nodes)
u = sp.zeros(nodes)

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

# f_g = sp.zeros(n_particle)

# --------------------------------------------------
# --------------------------------------------------

# ---------   Beginning of Time Loop  --------------

defPlot(posx_particle, posy_particle, n_particle, u, x, L, 0)      # Plot initial condition

for t in range(tempo):
    u_last[0] = v0
    u_last[-1] = vh
    RHS = sp.dot(Mdt, u_last) + B
    RHS[0] = v0
    RHS[-1] = vh
    u = sp.linalg.solve(LHS, RHS)

    for i in range(nodes):
        u[i] = (6/(h**2)) * (x[i] * h - x[i]**2)

    # Interpolate fluid velocity at n+1 and n time step
    u_last_interp = sp.interpolate.interp1d(x, u_last, fill_value=0, bounds_error=False)
    u_interp = sp.interpolate.interp1d(x, u, fill_value=0, bounds_error=False)

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
        vy[i] = vy_last[i] + (f_lift[i] - f_g[i] + f_dragy[i]) * constant
        vx[i] = vx_last[i] + (f_mass_partialx[i] + f_dragx[i]) * constant

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

    u_last = sp.copy(u)

    defPlot(posx_particle, posy_particle, n_particle, u, x, L, t+1)

# ---------   End of Time Loop   -------------------

# --------------------------------------------------
# --------------------------------------------------

