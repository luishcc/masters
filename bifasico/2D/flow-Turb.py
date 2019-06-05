import scipy as sp
from scipy import linalg
import os
import sys
import csv
import matrix as ma
import matplotlib.pyplot as plt

cwd = os.getcwd()

# --------------------------------------------------
#   Problem Parameters

dt = 10000
tempo = 1000

# fluid
rho_fld = 1000.0
viscosity_din = 0.8e-3
viscosity_kin = viscosity_din / rho_fld
grad_p = -0.012

# boundary
v0 = 0
vh = 0

