#-*- coding: utf-8 -*-
#-------------------------------------------------------------
# %matplotlib inline
# from IPython import get_ipython;   
# get_ipython().magic('reset -sf')
# from __future__ import unicode_literals
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 0.1 #set the value globally
mpl.rcParams['lines.dashed_pattern'] = [6, 6]
mpl.rcParams['lines.scale_dashes'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import warnings
warnings.simplefilter('ignore')
from scipy import sparse
import numpy as np
import scipy.sparse as sps
from matplotlib import pyplot as plt
from scipy.sparse.linalg.dsolve import linsolve

A=np.loadtxt('linha-fortran.txt', delimiter=',')
A1=np.loadtxt('coluna-fortran.txt', delimiter=',')
 
 
fig1, ax = plt.subplots(1,1,figsize=(3.5,3), frameon=True)
ax.locator_params(axis = 'y',nbins=10)
ax.locator_params(axis = 'x',nbins=10)
ax.tick_params(direction='in', width='0.1')
plt.tick_params(axis='both', labelsize=12, pad=1)
plt.xlabel(r'N', fontsize=12,labelpad=2)
plt.ylabel(r'CPU time $[s]$',style='italic', fontsize=12, rotation='vertical',
           labelpad=12)
ax.grid(color='grey', linestyle='--', linewidth=0.1,
        drawstyle='steps')    
plt.plot(A[:,1],A[:,0],'r-', linewidth=0.75,label='row')
plt.plot(A1[:100,1],A1[:100,0],'b-', linewidth=0.75,label='column')
plt.legend(title='',loc=1, borderaxespad=0., fontsize=12)
legend = plt.legend(loc=2, fontsize=10, facecolor='white',frameon=True,
                   fancybox=False, framealpha=1, shadow=False, borderpad=0.2,
                    labelspacing=0.5)
legend.get_frame().set_linewidth(0.1)
legend.get_frame().set_edgecolor("black")
plt.title('Fortran 90')
ax.set_xscale('log')
plt.show()
fig1.savefig('N_vs_CPU_Fortran.pdf',
             format='pdf',dpi=300, bbox_inches='tight',
             transparent=True)
