# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import sys
import shutil
import csv

cwd = os.getcwd()
if os.path.isfile(cwd+'/flowResult.csv') == True:
    os.remove(cwd+'/flowResult.csv')


# --------------------------------------------------
#   Problem Parameters

dt = 0.01
dt_inv = 1 / dt
tempo = 1000

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -12.


# --------------------------------------------------
#   Fluid Mesh and Matrices Generation

h = 1   # Duct height
L = 5*h     # Duct length
fine = 20
coarse = 10
d_fine = 0.15*h
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
            M[i_global, j_global] += m[i_local, j_local] * (dx[elem])/6

Mdt = M * dt_inv

# --------------------------------------------------
#   Allocating Fluid  Variables

u_last = sp.zeros(nodes)
u = sp.zeros(nodes)
viscosity_turb = sp.zeros(elem)


# --------------------------------------------------
# --------------------------------------------------
# ---------   Beginning of Time Loop  --------------


for t in range(tempo):

    with open('flowResult.csv', mode='a') as flowResult:
        writer = csv.writer(flowResult, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow([t, u_last])

    for elem in range(elem):
        for i_local in range(2):
            i_global = ien[elem, i_local]
            for j_local in range(2):
                j_global = ien[elem, j_local]
                K[i_global, j_global] += k[i_local, j_local] * (viscosity_kin + viscosity_turb[elem]) / dx[elem]

    #  Dirichlet Boundary Condition and System Definition
    LHS = Mdt + K
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

    u_last[0] = v0
    u_last[-1] = vh
    RHS = sp.dot(Mdt, u_last) + B
    RHS[0] = v0
    RHS[-1] = vh
    u = sp.linalg.solve(LHS, RHS)
    u_last = sp.copy(u)


    # for i in range(nodes):
    #     u[i] = (6/(h**2)) * (x[i] * h - x[i]**2)



# ---------   End of Time Loop   -------------------
# --------------------------------------------------
# --------------------------------------------------

# with open('flowResult.csv') as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter=',')
#     line_count = 0
#     for row in csv_reader:
#         if line_count == 0:
#             print(f'Column names are {", ".join(row)}')
#             line_count += 1
#         else:
#             print(f'\t{row[0]} , {row[1][500]} ')
#             line_count += 1
#     print(f'Processed {line_count} lines.')
