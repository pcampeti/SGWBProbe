#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:00:48 2020

@author: paolo
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import matplotlib as mpl
import seaborn as sns

# seaborn settings
sns.set()
sns.set(style='whitegrid')

# mpl settings
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['figure.figsize'] = [5,3]
mpl.rcParams['text.usetex'] = True
mpl.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})

# SKA
SKA = np.load('/home/paolo/Codes/SGWBProbeComb/files/hc_SKA.npz')
SKA_freq = SKA['x']
SKA_hc = SKA['y']
SKA_strain = SKA_hc**2/SKA_freq

# Einstein Telescope
ET = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_ET.npz')
ET_freq = ET['x']
ET_strain = ET['y']

# Advanced LIGO
aLIGO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_aLIGO.npz')
aLIGO_freq = aLIGO['x']
aLIGO_strain = aLIGO['y']

# LISA
LISA = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_LISA_xcosmo.npz')
LISA_freq = LISA['x']
LISA_strain = LISA['y']

# muAres 
Ares = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_muAres_nofgs.npz')
Ares_freq = Ares['x']
Ares_strain = Ares['y']

# BBO 
BBO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_BBO_STAR.npz')
BBO_freq = BBO['x']
BBO_strain = BBO['y']

# DECIGO
DECIGO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DECIGO.npz')
DECIGO_freq = DECIGO['x']
DECIGO_strain = DECIGO['y']

# DO Optimal
DO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DO_Optimal.npz')
DO_freq = DO['x']
DO_strain = DO['y']

# DO Conservative
DO_cons = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DO_Conservative.npz')
DO_cons_freq = DO_cons['x']
DO_cons_strain = DO_cons['y']

# AEDGE
AEDGE = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_AEDGE.npz')
AEDGE_freq = AEDGE['x']
AEDGE_strain = AEDGE['y']

###############################################################################
# Plot strain curves
fig = plt.figure()
ax = plt.gca()
ax.plot(ET_freq, np.sqrt(ET_freq*(1/np.sqrt(3))*ET_strain), label='ET', linewidth='1.5', color=sns.xkcd_rgb["windows blue"])
ax.plot(aLIGO_freq, np.sqrt(aLIGO_freq*aLIGO_strain), label='aLIGO', linewidth='1.5', color=sns.xkcd_rgb["yellow"])
ax.plot(DECIGO_freq, np.sqrt(DECIGO_freq*DECIGO_strain), label='DECIGO', linewidth='1.5', color=sns.xkcd_rgb["amber"])
ax.plot(DO_freq, np.sqrt(DO_freq*DO_strain), label='DO Optimal', linewidth='1.5', color=sns.xkcd_rgb["greyish"])
ax.plot(DO_cons_freq, np.sqrt(DO_cons_freq*DO_cons_strain), label='DO Conservative', linewidth='1.5', color=sns.xkcd_rgb["faded green"])
ax.plot(BBO_freq, np.sqrt(BBO_freq*BBO_strain), label='BBO', linewidth='1.5', color=sns.xkcd_rgb["dusty purple"])
ax.plot(Ares_freq, np.sqrt(Ares_freq*Ares_strain), label='$\mu$Ares', linewidth='1.5', color=sns.xkcd_rgb["cyan"])
ax.plot(AEDGE_freq, np.sqrt(AEDGE_freq*AEDGE_strain), label='AEDGE', linewidth='1.5', color=sns.xkcd_rgb["pale red"])
ax.loglog(LISA_freq, np.sqrt(LISA_freq*LISA_strain), label='LISA', linewidth='1.5', color=sns.xkcd_rgb["black"])
ax.plot(SKA_freq, np.sqrt(SKA_freq*SKA_strain), label='SKA', linewidth='1.5', color=sns.xkcd_rgb["violet"])
ax.set_ylim([1e-25, 1e-12])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both', which = 'major', labelsize = 8.0)
ax.legend(fontsize=7, bbox_to_anchor=(1., 0.8))
plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_3.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
