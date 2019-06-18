# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp

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

def neighbourElements(_np, _ien):
    result = [None] * _np
    for i in range(_np):
        result[i] = []
        for e in range(len(_ien)):
            for j in range(3):
                if _ien[e, j] == i:
                    result[i].append(e)

    return result



def neighbourTest(_np):
    result = [None] * _np
    for i in range(_np):
        result[i] = []
    #
    # result = np.empty((_np,), dtype=object)
    # for i in range(_np):
    #     result[i] = []

    for j in range(3):
        for i in range(_np):
            result[i].append(j)

    result[3].append(9)

    print(result)
    return result


def fem_matrix(_x, _y, _numele, _numnode, _ax, _ien):
    k_local = sp.zeros((3, 3))
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
    a = sp.zeros(3)
    b = sp.zeros(3)
    c = sp.zeros(3)
    yy = sp.zeros(3)
    xx = sp.zeros(3)
    K_global = sp.zeros((_numnode, _numnode))
    M_global = sp.zeros((_numnode, _numnode))
    visc = 0.

    for elem in range(_numele):
        for i in range(3):
            visc += _ax[i]
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

        visc = visc / 3.0

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]
        Area = (a[0] + a[1] + a[2]) * 0.5

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                K_global[i_global, j_global] += k_local[i_local, j_local] * visc
                M_global[i_global, j_global] += m_local[i_local, j_local] * (Area / 12.)

    return K_global, M_global
