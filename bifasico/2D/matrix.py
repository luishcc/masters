import scipy as sp


def turbulentMatrix(_ien, _dx, _k, _mdt, _b, _viscosity_kin, _viscosity_turb, _bC):
    elem = len(_ien)
    K = sp.zeros((elem+1, elem+1))
    for e in range(elem):
        for i_local in range(2):
            i_global = _ien[e, i_local]
            for j_local in range(2):
                j_global = _ien[e, j_local]
                K[i_global, j_global] += _k[i_local, j_local] * (_viscosity_kin + _viscosity_turb[e]) / _dx[e]

    #  Dirichlet Boundary Condition and System Definition
    LHS = _mdt + K
    LHS_copy = sp.copy(LHS)

    for i in range(elem+1):
        _b[i] = _b[i] - _bC[0] * LHS_copy[i, 0]
        _b[i] = _b[i] - _bC[1] * LHS_copy[i, -1]
        LHS[0, i] = 0
        LHS[i, 0] = 0
        LHS[-1, i] = 0
        LHS[i, -1] = 0

    return LHS, _b
