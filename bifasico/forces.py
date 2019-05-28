def get_vel(_pos, _fluid_vel, _x, _y, _rp):
    fluid_vel_part = sp.zeros(5)
    f = sp.interpolate.bisplrep(_x, _y, _fluid_vel)
    fluid_vel_part[0] = sp.interpolate.bisplev(_pos[0]-_rp, _pos[1], f)       # vx esquerda
    fluid_vel_part[1] = sp.interpolate.bisplev(_pos[0]+_rp, _pos[1], f)       # vx direita
    fluid_vel_part[2] = sp.interpolate.bisplev(_pos[0], _pos[1]+_rp, f)       # vy cima
    fluid_vel_part[3] = sp.interpolate.bisplev(_pos[0], _pos[1]-_rp, f)       # vy baixo
    fluid_vel_part_center = sp.interpolate.bisplev(_pos[0], _pos[1], f)     # v centro
    return fluid_vel_part, fluid_vel_part_center

def gravity(_rho, rp):
    f_g = (4/3.0) * sp.pi * (rp**3) * 9.81 * _rho
    return f_g

def drag(_relative_vel, _visc, _rp):
    f_drag = 3 * sp.pi * _visc * 2 * _rp * _relative_vel
    return f_drag

def lift():
    return f_lift

def virtualMass():
    return f_mass