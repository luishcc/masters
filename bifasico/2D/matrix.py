import scipy as sp


def turbulentMatrix(_ien, _dx, _k, _Mdt, _B, _viscosity_kin, _viscosity_turb, _BC):
    elem = len(ien)
    for elem in range(elem):
        for i_local in range(2):
            i_global = _ien[elem, i_local]
            for j_local in range(2):
                j_global = _ien[elem, j_local]
                K[i_global, j_global] += k[i_local, j_local] * (viscosity_kin + _viscosity_turb[elem]) / dx[elem]

    #  Dirichlet Boundary Condition and System Definition
    LHS = Mdt + K
    LHS_copy = sp.copy(LHS)

    for i in range(_nodes):
        B[i] = B[i] - _BC[0] * LHS_copy[i, 0]
        B[i] = B[i] - _BC[1] * LHS_copy[i, -1]
        LHS[0, i] = 0
        LHS[i, 0] = 0
        LHS[-1, i] = 0
        LHS[i, -1] = 0

    return LHS, B

def laminarMatrix(_ien, _dx, _dt, _viscosity_kin, _BC, _nodes, _grad):
    elem = len(_ien)

    k = sp.array([[1, -1], [-1, 1]])
    m = sp.array([[2, 1], [1, 2]])
    b = sp.array([1, 1])
    B = sp.zeros(_nodes)
    K = sp.zeros((_nodes, _nodes))
    M = sp.zeros((_nodes, _nodes))

    for elem in range(elem):
        for i_local in range(2):
            i_global = _ien[elem, i_local]
            B[i_global] += b[i_local] * (_grad * _dx[elem] * 0.5)
            for j_local in range(2):
                j_global = _ien[elem, j_local]
                M[i_global, j_global] += m[i_local, j_local] * (_dx[elem])/6
                K[i_global, j_global] += k[i_local, j_local] * _viscosity_kin / _dx[elem]

    #  Dirichlet Boundary Condition and System Definition
    Mdt = M * (1./_dt)
    LHS = Mdt + K
    LHS_copy = sp.copy(LHS)
    for i in range(_nodes):
        B[i] = B[i] - _BC[0] * LHS_copy[i, 0]
        B[i] = B[i] - _BC[1] * LHS_copy[i, -1]
        LHS[0, i] = 0
        LHS[i, 0] = 0
        LHS[-1, i] = 0
        LHS[i, -1] = 0

    return LHS, B, Mdt
