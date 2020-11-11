#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 23:27:13 2020

@author: Paolo Campeti

This script reproduces Figure 1 in the paper. 
Uses methods imported from module sgwbprobecomb/SGWB_Signal.py.

"""
import os.path as op
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
# import our methods
from sgwbprobe.SGWB_Signal import Signal_GW

# seaborn settings 
sns.set()
sns.set(style='whitegrid')

# matplotlib settings
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['figure.figsize'] = [5,3]
mpl.rcParams['text.usetex'] = True
mpl.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})
axissize = 6
labelsize = 8
legendsize = 10
colornorm = colors.Normalize(vmin=0.0, vmax=5.0)
linesize = 2

# Generate primordial signals
class_axion_old = Signal_GW(r_vac=1e-5, r_star=400, k_p=1e15, sigma=9.1, axion=True, running=True)
class_axion1 = Signal_GW(r_vac=1e-5, r_star=50, k_p=1e6, sigma=4.8, axion=True, running=True)
class_axion2 = Signal_GW(r_vac=1e-5, r_star=0.15 , k_p=1e11, sigma=8, axion=True, running=True)
class_no_axion = Signal_GW(r_vac=0.01, axion=None, running=True)
class_no_axion_r0001 = Signal_GW(r_vac=0.001, axion=None, running=True)
class_no_axion_rBICEP = Signal_GW(r_vac=0.06, axion=None, running=True)

# wavenumber array
k = np.logspace(np.log10(1e-5), np.log10(1e20), 100000)

# Plot tensor power spectra models (Figure 1) 
fig = plt.figure()
plt.loglog(k, class_axion_old.total_spect(k),label='Axion-SU(2) $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=linesize, color=sns.xkcd_rgb["violet"])
plt.loglog(k, class_axion1.total_spect(k),label='Axion-SU(2) $r_{\star}=50$, $k_{p}=10^{6}$ $Mpc^{-1}$, $\sigma=4.8$',linewidth=linesize, color=sns.xkcd_rgb["blue"])
plt.plot(k, class_axion2.total_spect(k),label='Axion-SU(2) $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$',linewidth=linesize, color=sns.xkcd_rgb["light orange"])
plt.plot(k, class_no_axion_rBICEP.tensor_spect(k),label='Single-Field Slow-Roll r=0.06',linewidth=linesize, color=sns.xkcd_rgb["black"], linestyle='--')
plt.plot(k, class_no_axion.tensor_spect(k),label='Single-Field Slow-Roll r=0.01',linewidth=linesize, color=sns.xkcd_rgb["teal"])
plt.plot(k, class_no_axion_r0001.tensor_spect(k),label='Single-Field Slow-Roll r=0.001',linewidth=linesize, color=sns.xkcd_rgb["red"])
plt.xlabel(r'Wavenumber k $[Mpc^{-1}]$',fontsize = labelsize)
plt.ylabel(r'$\mathcal{P}_{T}(k)$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize = 7, loc='upper right')
axes = plt.gca()
axes.set_xlim([1e-4,1e20])
axes.set_ylim([5e-14,1e-1])
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_1.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()