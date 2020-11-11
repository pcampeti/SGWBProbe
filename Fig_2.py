#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:52:46 2020

@author: Paolo Campeti

This script reproduces Figure 2 in the paper. 

"""
import os.path as op
import numpy as np
import pylab as pl
from fgbuster.cosmology import _get_Cl_cmb
import seaborn as sns
import matplotlib as mpl

# seaborn settings
sns.set()
sns.set(style='whitegrid')

#mpl settings
mpl.rcParams['text.usetex'] = True
mpl.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})

# number of simulations
N_sim = 132
# min, max multipole ell
lmin, lmax = 2, 200

###############################################################################
# Plot Fig.2

ell_v = np.arange(lmin,lmax)
norm = ell_v*(ell_v+1)/2/np.pi

average_white = np.zeros((lmax))
fig, ax = pl.subplots()    
for i in range(0, N_sim):
    if i==0:
        ClBB_res_noise_white = np.load(op.join(op.dirname(__file__), 'files/LiteBIRD_fgs_residuals/ClBB_res_noise_white_'+str(i)+'.npy'))
        ax.loglog(norm*ClBB_res_noise_white[lmin:lmax], color='DarkOrange', linewidth=1.0, alpha=1.0, label='Realizations of Res. Foregrounds +Post-Comp.Sep Noise')
    if i != 58 and i != 110: #jump missing noise simulation
        ClBB_res_noise_white = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_fgs_residuals/ClBB_res_noise_white_'+str(i)+'.npy'))
        average_white += ClBB_res_noise_white[:lmax]
        ax.loglog(norm*ClBB_res_noise_white[lmin:lmax], color='DarkOrange', linewidth=1.0, alpha=1.0)
        
average_white = average_white / (132-2) # compute average residuals (remember to subtract the two missing simulations at denominator)

ax.loglog(norm*average_white[lmin:lmax], color='red', linewidth=2.0, alpha=1.0, label='Average Res. Foregrounds+Post-Comp.Sep. Noise')

ax.loglog(norm*_get_Cl_cmb(1.,0.)[2,lmin:lmax], c='black', label='Lensing')
ax.loglog(norm*_get_Cl_cmb(0.,0.01)[2,lmin:lmax], c='grey', ls='--', label='r = 0.01')
ax.loglog(norm*_get_Cl_cmb(0.,0.001)[2,lmin:lmax], c='grey', ls='-', label='r = 0.001')
ax.set_ylim(1e-6,1e0)
ax.legend(loc='upper left', fontsize=11)
ax.set_xlim([lmin, 200])
ax.set_xlabel('$\ell$', fontsize=15)
ax.set_ylabel('$\ell(\ell+1)C_\ell / 2\pi$', fontsize=15)
pl.savefig(op.join(op.dirname(__file__), 'figures/Fig_2.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
pl.show()


