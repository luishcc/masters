import scipy as sp
from scipy import linalg
import InOut as io
import semiLagrangean as sl
import libmalha as lm
import GMesh as gm
import os
import forces
import random


def get_vel(_pos, _fluid_vel, _x, _y, _rp):
    fluid_vel_part = sp.zeros(5)
    f = sp.interpolate.bisplrep(_x, _y, _fluid_vel)
    fluid_vel_part[0] = sp.interpolate.bisplev(_pos[0]+_rp, _pos[1], f)       # vx direita
    fluid_vel_part[1] = sp.interpolate.bisplev(_pos[0]-_rp, _pos[1], f)       # vx esquerda
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

def lift(_fluid_vel, _fluid_vel_corners, _particle_velocity,
         _particle_density, _particle_radius, _fluid_viscosity):

    relative_vel = _fluid_vel - _particle_velocity
    a = 1.61 * sp.sqrt(_particle_density * _fluid_viscosity) * (_particle_radius*2)**2
    du_dy = (_fluid_vel_corners[0] - _fluid_vel_corners[1]) / (_particle_radius*2)
    f_lift = a * relative_vel * du_dy * sp.sqrt(abs(du_dy))
    return f_lift

def virtualMass(_fluid_vel, _fluid_vel_last, _fluid_density, _particle_radius, _dt):

    volume = (4/3.) * sp.pi * (_particle_radius ** 3)
    f_mass = _fluid_density * volume * 0.5 * ((_fluid_vel - _fluid_vel_last) / _dt)
    return f_mass

