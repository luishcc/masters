# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
import sys
import csv
import matplotlib.pyplot as plt

cwd = os.getcwd()

# --------------------------------------------------
#   Problem Parameters

dt = 10
tempo = 10000

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -12

# boundary
v0 = 0
vh = 0

h = 1.0
L = 8*h     # Duct length

fine = 4000
coarse = 2000
d_fine = 0.15*h
d_coarse = 0.85*h

y_fine = sp.linspace(0, d_fine, fine, endpoint=False)
y_coarse = sp.linspace(d_fine, d_coarse, coarse, endpoint=False)
y_fine2 = sp.linspace(d_coarse, h, fine+1)
y = sp.append(y_fine, y_coarse)
y = sp.append(y, y_fine2)

# y = sp.linspace(0, h, 1000)

nodes = len(y)
elem = nodes-1

dy = sp.zeros(elem)

Q = sp.ones(nodes) * (-grad_p / rho_fld)

row = sp.zeros(elem*4)
col = sp.zeros(elem*4)
k_d = sp.zeros(elem*4)
m_d = sp.zeros(elem*4)
M_lump = sp.zeros(nodes)
g_d = sp.zeros(elem*4)

print "assembly"
for e in range(0, elem):
    p = e * 4
    dy[e] = abs(y[e+1] - y[e])

    row[p] = e
    row[p+1] = e
    row[p+2] = e + 1
    row[p+3] = e + 1

    col[p] = e
    col[p+1] = e + 1
    col[p+2] = e
    col[p+3] = e + 1

    m_d[p] = 2 * (dy[e] / 6.)
    m_d[p+1] = 1 * (dy[e] / 6.)
    m_d[p+2] = 1 * (dy[e] / 6.)
    m_d[p+3] = 2 * (dy[e] / 6.)

    M_lump[e] += m_d[p] + m_d[p+1]
    M_lump[e+1] += m_d[p+2] + m_d[p+3]

    g_d[p] = -0.5
    g_d[p+1] = 0.5
    g_d[p+2] = -0.5
    g_d[p+3] = 0.5

M = sp.sparse.coo_matrix((m_d, (row, col)), shape=(nodes, nodes))
G = sp.sparse.coo_matrix((g_d, (row, col)), shape=(nodes, nodes))

G = G.tocsr()
M = M.tocsr()

Mdt = M / dt


# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)
du_dy = sp.zeros(nodes)
viscosity_turb = sp.zeros(nodes)

y_plus = sp.zeros(nodes)
u_plus_a = sp.zeros(nodes)
lc = sp.zeros(nodes)


M_inv = sp.zeros(nodes)

for i in range(nodes):
    if 0 <= y[i] <= 0.1 * h:
        lc[i] = 0.41 * y[i]
    elif h >= y[i] >= 0.9 * h:
        lc[i] = 0.41 * (h - y[i])
    else:
        lc[i] = 0.41 * 0.1 * h

    M_inv[i] = 1. / M_lump[i]

print "M inv lumped"


# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------


# with open('turbulent.csv', mode='w') as turbulent:
#     writer = csv.writer(turbulent, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    # writer.writerow(viscosity_turb, u_plus, u_plus_a, y_plus, lc)
    # writer.writerow(lc)


# with open('flow-T.csv', mode='w') as flowResultT:
#     writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     writer.writerow(y)

A_const = 26
for t in range(tempo):
    print t, " / ", tempo
    # with open('flow-T.csv', mode='a') as flowResultT:
    #     writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     writer.writerow(u_last)

    u[0] = 0
    u[-1] = 0

    t_wall = viscosity_din * (0.5/dy[0]) * (-3 * u[0] + 4 * u[1] - u[2])
    u_t = sp.sqrt(t_wall / rho_fld)
    y_plus = (u_t / viscosity_kin) * y
    u_plus = (1./u_t) * u
    print t_wall

    du_dy = G.dot(u)
    for i in xrange(nodes):
        if 0 <= y[i] <= 0.1 * h:
            D = 1 - sp.exp(-y_plus[i]/A_const)
            lc[i] = 0.41 * y[i] * D
        elif h >= y[i] >= 0.9 * h:
            index = nodes - i - 1
            D = 1 - sp.exp(-y_plus[index]/A_const)
            lc[i] = 0.41 * (h - y[i]) * D

        u_plus_a[i] = (1. / 0.41) * sp.log(y_plus[i]) + 5.5
        viscosity_turb[i] = (lc[i]**2) * abs(du_dy[i])

    viscosity_turb = sp.multiply(M_inv, viscosity_turb)

    for e in xrange(elem):
        p = 4 * e

        visc = (viscosity_turb[e] + viscosity_turb[e+1]) * 0.5

        k_d[p] = 1 * (1. / dy[e]) * (viscosity_kin + visc)
        k_d[p + 1] = -1 * (1. / dy[e]) * (viscosity_kin + visc)
        k_d[p + 2] = -1 * (1. / dy[e]) * (viscosity_kin + visc)
        k_d[p + 3] = 1 * (1. / dy[e]) * (viscosity_kin + visc)

    K = sparse.coo_matrix((k_d, (row, col)), shape=(nodes, nodes))
    K = K.tocsr()

    LHS = Mdt + K
    cc = sp.zeros(nodes)
    LHS_copy = sp.copy(LHS)

    # for i in range(nodes):
    #     LHS[0, i] = 0
    #     LHS[i, 0] = 0
    #     LHS[-1, i] = 0
    #     LHS[i, -1] = 0


    LHS[0, 0] = 1
    LHS[1, 0] = 0
    LHS[-1, -1] = 1
    LHS[-2, -1] = 1


    u_last[0] = v0
    u_last[-1] = vh

    RHS = Mdt.dot(u_last) + M.dot(Q)

    RHS[0] = v0
    RHS[-1] = vh

    u = spsolve(LHS, RHS)

    norm = u - u_last
    norm = sp.dot(norm, norm)
    print norm
    if norm < 10e-10:
        break

    u_last = sp.copy(u)

    # with open('turbulent.csv', mode='a') as turbulent:
    #     writer = csv.writer(turbulent, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        # writer.writerow([viscosity_turb, u_plus, u_plus_a, y_plus, lc])
        # writer.writerow(lc)

# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------

# with open('flow-T.csv', mode='a') as flowResultT:
#     writer = csv.writer(flowResultT, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     writer.writerow(u_last)
#
# with open('prop-T.csv', mode='w') as properties:
#     writer = csv.writer(properties, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     writer.writerow([dt, tempo, viscosity_din, rho_fld, viscosity_kin, h, L])


plt.figure(1)
plt.plot(y, u)
plt.title('Perfil de Velocidade')
plt.ylabel('u')
plt.xlabel('y')
plt.grid(True)

# plt.figure(4)
# plt.plot(lc, y)

# plt.figure(5)
# plt.plot(y, du_dy)

x = sp.linspace(0.1, 20, 10000)
plt.figure(3)
plt.semilogx(x, x, label='u+ = y+')
plt.semilogx(y_plus[0:nodes/2], u_plus_a[0:nodes/2], 'k--', label='Lei da parede')
plt.semilogx(y_plus[0:nodes/2], u_plus[0:nodes/2], 'y.', label=u'numérico')
plt.ylabel("u+")
plt.xlabel("y+")
plt.xlim(0.1, 10**3)
plt.ylim(0, 25)
plt.grid(True)
plt.legend()

plt.figure(2)
plt.plot(y, viscosity_turb)
plt.title('Viscosidade Turbulenta')
plt.xlabel("y")
plt.ylabel('nu_t')
plt.grid(True)
plt.show()
