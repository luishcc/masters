# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
from scipy import interpolate
import matplotlib.pyplot as plt

cwd = os.getcwd()
)

os.mkd



# with open('flowResult.csv', 'r') as csvfile:
#     plots = csv.reader(csvfile, delimiter=',')
#     for row in plots:

matriz = sp.loadtxt('flow-T.csv', delimiter=',')


x = sp.linspace(0,1, len(matriz[-1]))
plt.figure(1)
plt.plot(matriz[0], matriz[-1])
plt.show()





