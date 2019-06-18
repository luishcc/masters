# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import os
import random
from scipy import interpolate
import matplotlib.pyplot as plt
import sys
import shutil
import csv

cwd = os.getcwd()
if os.path.isdir(cwd+'/result') == True:
    shutil.rmtree(cwd+'/result')

os.mkdir(cwd+'/result')

tempo = 1001
nodes = 1000
u = sp.zeros((tempo, nodes))




# with open('flowResult.csv', 'r') as csvfile:
#     plots = csv.reader(csvfile, delimiter=',')
#     for row in plots:

matriz = sp.loadtxt('flow-T.csv', delimiter=',')



x = sp.linspace(0,1, len(matriz[-1]))
plt.figure(1)
plt.plot(matriz[0], matriz[-1])
plt.show()



exit()

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



