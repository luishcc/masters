# -*- coding: utf-8 -*-
import numpy as np

# Escrito por: Luís Cunha
# Última atualização: 14/08/2017
#
# Biblioteca para geração e manipulação de malhas computacionais
#
# Malhas 1D ([numero de elementos], [comprimento da malha]):
#       Retornam: Vetor de posições x
#                 Matriz de conectividade IEN
#
# 1) lin1d - Malha 1D linear
# 2) quad1d - Malha 1D quadrática
# 3) cube1d - Malha 1D cúbica
# 4) exp1d - Malha 1D exponencial
#
# Malhas 2D ([numero de elementos x], [comprimento da malha x],
#                                [numero de elementos y], [comprimento da malha y]):
#       Retornam: Vetor de posições x, y
#                 Matriz de conectividade IEN retângular
#
# 5) lin2d - Malha 2D linear
# 6) quad2d - Malha 2D quadrática
# 7) cube2d - Malha 2D cúbica
# 8) exp2d - Malha 2D exponencial
#
#   Outros:
#
# 9) ret_tri([matriz de conectividade retangular])
#        Converte matriz de conectividade retangular para triangular
#
# 10) combine1d([distribuição em x], [distribuição em y])
#        Combina dois vetores de posição em uma malha 2D

def lin1d(num_ele, L):
    dx = float(L) / num_ele
    x = np.zeros(num_ele+1)
    IEN = np.zeros((num_ele, 2), dtype="int32")
    for i in range(0, num_ele):
        for j in [0, 1]:
            IEN[i, j] = i+j+1
    for i in range(0, num_ele+1):
        x[i] = i * dx
    return x, IEN

def quad1d(num_ele, L):
    LL = np.sqrt(L)
    dx = float(LL) / num_ele
    y = np.zeros(num_ele + 1)
    x = np.zeros(num_ele + 1)
    IEN = np.zeros((num_ele, 2), dtype="int32")
    for i in range(0, num_ele):
        for j in [0, 1]:
            IEN[i, j] = i + j + 1
    for i in range(0, num_ele + 1):
        x[i] = i * dx
        y[i] = x[i]**2
    return y, IEN

def cube1d(num_ele, L):
    LL = np.cbrt(L)
    dx = float(LL) / num_ele
    y = np.zeros(num_ele + 1)
    x = np.zeros(num_ele + 1)
    IEN = np.zeros((num_ele, 2), dtype="int32")
    for i in range(0, num_ele):
        for j in [0, 1]:
            IEN[i, j] = i + j + 1
    for i in range(0, num_ele + 1):
        x[i] = i * dx
        y[i] = x[i]**3
    return y, IEN

def exp1d(num_ele, L):
    LL = np.log(L+1)
    dx = float(LL) / num_ele
    y = np.zeros(num_ele + 1)
    x = np.zeros(num_ele + 1)
    IEN = np.zeros((num_ele, 2), dtype="int32")
    for i in range(0, num_ele):
        for j in [0, 1]:
            IEN[i, j] = i + j + 1
    for i in range(0, num_ele + 1):
        x[i] = i * dx
        y[i] = np.exp(x[i])-1
    return y, IEN

def lin2d(num_x, Lx, num_y, Ly):
    dx = float(Lx) / num_x
    dy = float(Ly) / num_y
    x = np.zeros((num_x+1)*(num_y+1))
    y = np.zeros((num_x+1)*(num_y+1))
    IEN = np.zeros(((num_x)*(num_y), 4), dtype='int')
    k = 0
    for j in range(0, (num_x+1)*(num_y+1), num_x+1):
        for i in range(0, num_x+1):
            x[i+j] = i * dx
    for i in range(0, (num_x+1)*(num_y+1), num_x+1):
        for j in range(0, num_x+1):
            y[i+j] = k * dy
        k += 1
    for i in range(len(IEN)):
        IEN[i] = (i, i+1, i+(num_x+1)+1, i+(num_x+1))
        IEN[i] += i/(num_x)
    return x, y, IEN

def quad2d(num_x, Lx, num_y, Ly):
    LxL = np.sqrt(Lx)
    LyL = np.sqrt(Ly)
    dx = float(LxL) / num_x
    dy = float(LyL) / num_y
    x = np.zeros((num_x+1)*(num_y+1))
    y = np.zeros((num_x+1)*(num_y+1))
    x2 = np.zeros((num_x+1)*(num_y+1))
    y2 = np.zeros((num_x+1)*(num_y+1))
    IEN = np.zeros(((num_x)*(num_y), 4), dtype='int')
    k = 0
    for j in range(0, (num_x+1)*(num_y+1), num_x+1):
        for i in range(0, num_x+1):
            x2[i+j] = i * dx
            x[i+j] = x2[i+j]**2
    for i in range(0, (num_x+1)*(num_y+1), num_x+1):
        for j in range(0, num_x+1):
            y2[i+j] = k * dy
            y[i+j] = y2[i+j]**2
        k += 1
    for i in range(len(IEN)):
        IEN[i] = (i, i+1, i+(num_x+1)+1, i+(num_x+1))
        IEN[i] += i/(num_x)
    return x, y, IEN

def cube2d(num_x, Lx, num_y, Ly):
    LxL = np.cbrt(Lx)
    LyL = np.cbrt(Ly)
    dx = LxL / float(num_x)
    dy = LyL / float(num_y)
    x = np.zeros((num_x+1)*(num_y+1))
    y = np.zeros((num_x+1)*(num_y+1))
    x2 = np.zeros((num_x+1)*(num_y+1))
    y2 = np.zeros((num_x+1)*(num_y+1))
    IEN = np.zeros(((num_x)*(num_y), 4), dtype='int')
    k = 0
    for j in range(0, (num_x+1)*(num_y+1), num_x+1):
        for i in range(0, num_x+1):
            x2[i+j] = i * dx
            x[i+j] = x2[i+j]**3
    for i in range(0, (num_x+1)*(num_y+1), num_x+1):
        for j in range(0, num_x+1):
            y2[i+j] = k * dy
            y[i+j] = y2[i+j]**3
        k += 1
    for i in range(len(IEN)):
        IEN[i] = (i, i+1, i+(num_x+1)+1, i+(num_x+1))
        IEN[i] += i/(num_x)
    return x, y, IEN

def exp2d(num_x, Lx, num_y, Ly):
    LxL = np.log(Lx+1)
    LyL = np.log(Ly+1)
    dx = LxL / float(num_x)
    dy = LyL / float(num_y)
    x = np.zeros((num_x+1)*(num_y+1))
    y = np.zeros((num_x+1)*(num_y+1))
    x2 = np.zeros((num_x+1)*(num_y+1))
    y2 = np.zeros((num_x+1)*(num_y+1))
    IEN = np.zeros(((num_x)*(num_y), 4), dtype='int')
    k = 0
    for j in range(0, (num_x+1)*(num_y+1), num_x+1):
        for i in range(0, num_x+1):
            x2[i+j] = i * dx
            x[i+j] = np.exp(x2[i+j])-1
    for i in range(0, (num_x+1)*(num_y+1), num_x+1):
        for j in range(0, num_x+1):
            y2[i+j] = k * dy
            y[i+j] = np.exp(y2[i+j])-1
        k += 1
    for i in range(len(IEN)):
        IEN[i] = (i, i+1, i+(num_x+1)+1, i+(num_x+1))
        IEN[i] += i/(num_x)
    return x, y, IEN

def ret_tri(ien):
    num = len(ien)
    ient = np.zeros((num*2, 3), dtype="int32")
    ient[0:num] =  np.delete(ien, 3, axis=1)
    ient[num:len(ient)] = np.delete(ien, 1, axis=1)
    return ient

def combine1D(malha1, malha2):
    numx = len(malha1)
    numy = len(malha2)
    blankx = np.zeros(numx)
    blanky = np.zeros(numy)
    for i in range (numx):
        blankx[i] = malha1[i]
    xx = np.array(blankx)
    for i in range (numy):
        blanky[i] = malha2[i]
    yy = np.array(blanky)
    x, y = np.meshgrid(xx, yy)
    IEN = np.zeros(((numx)*(numy), 4), dtype='int')
    for i in range(len(IEN)):
        IEN[i] = (i, i+1, i+(numx+1)+1, i+(numx+1))
        IEN[i] += i/(numx)
    return x, y, IEN

