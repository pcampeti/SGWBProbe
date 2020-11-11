#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 17:00:48 2020

@author: Paolo Campeti

This script reproduces Figures 4 - 7 in the paper. 
Uses methods imported from module sgwbprobecomb/foregrounds.py.

"""
import os.path as op
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import matplotlib as mpl
import seaborn as sns
# import our methods 
from sgwbprobecomb.foregrounds import GWD_fg, EGWD_fg, BBH_BNS_fg, MBHB_fg


# seaborn settings
sns.set()
sns.set(style='whitegrid')

#mpl settings
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['figure.figsize'] = [5,3]
mpl.rcParams['text.usetex'] = True
mpl.rc('font',**{'family':'serif','serif':['Times New Roman'],'size':14})

# Useful constants: seconds in a year
year_sec = 60*60*24*365 

# LISA 
LISA = np.load(op.join(op.dirname(__file__), 'files/S_h_LISA_xcosmo.npz'))
LISA_freq = LISA['x']
LISA_strain = LISA['y']

# muAres 
Ares = np.load(op.join(op.dirname(__file__), 'files/S_h_muAres_nofgs.npz'))
Ares_freq = Ares['x']
Ares_strain = Ares['y']

# BBO
BBO = np.load(op.join(op.dirname(__file__), 'files/S_h_BBO_STAR.npz'))
BBO_freq = BBO['x']
BBO_strain = BBO['y']

# DECIGO
DECIGO = np.load(op.join(op.dirname(__file__), 'files/S_h_DECIGO.npz'))
DECIGO_freq = DECIGO['x']
DECIGO_strain = DECIGO['y']

# DO Optimal
DO = np.load(op.join(op.dirname(__file__), 'files/S_h_DO_Optimal.npz'))
DO_freq = DO['x']
DO_strain = DO['y']

# DO Conservative
DO_cons = np.load(op.join(op.dirname(__file__), 'files/S_h_DO_Conservative.npz'))
DO_cons_freq = DO_cons['x']
DO_cons_strain = DO_cons['y']

# AEDGE
AEDGE = np.load(op.join(op.dirname(__file__), 'files/S_h_AEDGE.npz'))
AEDGE_freq = AEDGE['x']
AEDGE_strain = AEDGE['y']

# SKA
SKA = np.load(op.join(op.dirname(__file__), 'files/hc_SKA.npz'))
SKA_freq = SKA['x']
SKA_hc = SKA['y']

# AEDGE
ET = np.load(op.join(op.dirname(__file__), 'files/S_h_ET.npz'))
ET_freq = ET['x']
ET_strain = ET['y']

###############################################################################    
# Plot strain sensitivity for SKA (Fig.7, right panel)
fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-10), np.log10(1e3), 10000)

ax.loglog(SKA_freq, SKA_hc, linewidth='1.5', color=sns.xkcd_rgb["steel blue"], linestyle='--', label='SKA strain noise')
ax.plot(freq, MBHB_fg(freq), linewidth='2.0', color='purple', linestyle='-.', label='MBHB Foreground')
ax.set_xlim([1e-10, 1e-6])
ax.set_ylim([1e-17, 1e-12])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10)
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_7_right.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
# Plot strain sensitivity for LISA (Fig. 4, left panel)
fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-7), np.log10(1e3), 10000)

ax.loglog(freq, np.sqrt(freq*GWD_fg(freq, 3*year_sec)), label='Galactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["electric blue"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.plot(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.loglog(LISA_freq, np.sqrt(LISA_freq*LISA_strain), linewidth='1.5', color=sns.xkcd_rgb["black"], linestyle='--', label='LISA strain noise')
ax.set_xlim([1e-4, 1e-1])
ax.set_ylim([1e-23, 1e-18])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10)
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_4_left.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
###############################################################################
# Plot strain sensitivity for DO (Fig. 4, right panel)
fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-7), np.log10(1e3), 10000)

ax.loglog(freq, np.sqrt(freq*GWD_fg(freq, 3*year_sec)), label='Galactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["electric blue"])
ax.plot(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.plot(DO_freq, np.sqrt(DO_freq*DO_strain), linewidth='1.5', color=sns.xkcd_rgb["greyish"], zorder=4, linestyle='--', label='DO Optimal strain noise')
ax.plot(DO_cons_freq, np.sqrt(DO_cons_freq*DO_cons_strain), linewidth='1.5', color=sns.xkcd_rgb["dark cyan"], zorder=5, linestyle='--', label='DO Conservative strain noise')
ax.set_xlim([1e-3, 1e1])
ax.set_ylim([1e-24, 1e-18])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10, loc='upper right')
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_4_right.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
# Plot strain sensitivity for AEDGE (Fig. 5, left panel)
fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-7), np.log10(1e3), 10000)

ax.loglog(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.plot(AEDGE_freq, np.sqrt(AEDGE_freq*AEDGE_strain), linewidth='1.5', color=sns.xkcd_rgb["denim blue"], linestyle='--', label='AEDGE strain noise')
ax.set_xlim([5e-3, 5e0])
ax.set_ylim([5e-24, 1e-21])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10, loc='upper right')
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_5_left.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
# Plot strain sensitivity for muAres (Fig. 5, right panel)

fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-10), np.log10(1e3), 10000)

ax.loglog(freq, np.sqrt(freq*GWD_fg(freq, 4*year_sec)), label='Galactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["electric blue"])
ax.plot(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.loglog(Ares_freq, np.sqrt(Ares_freq*Ares_strain), linewidth='1.5', color=sns.xkcd_rgb["amber"], linestyle='--', label=r'$\mu$Ares strain noise')

freq_mbhb = np.logspace(np.log10(1e-10), np.log10(2.36e-5), 10000)

ax.plot(freq_mbhb, MBHB_fg(freq_mbhb), linewidth='2.0', color='purple', linestyle='-.', label='MBHB Foreground')
ax.set_xlim([1e-6, 1e-1])
ax.set_ylim([1e-23, 5e-16])

plt.xlabel(r'f $[Hz]$',fontsize = 7.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10, loc='upper right')
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_5_right.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
# Plot strain sensitivity for DECIGO (Fig. 6, left panel)

fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-7), np.log10(1e3), 10000)

ax.loglog(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.plot(DECIGO_freq, np.sqrt(DECIGO_freq*DECIGO_strain), linewidth='1.5', color=sns.xkcd_rgb["purplish blue"], linestyle='--', label='DECIGO strain noise')
ax.set_xlim([1e-4, 1e2])
ax.set_ylim([1e-25, 1e-19])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10, loc='upper right')
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_6_left.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
# Plot strain sensitivity for BBO (Fig. 6, right panel)

fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-5), np.log10(1e2), 1000)

ax.loglog(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.plot(freq, np.sqrt(freq*EGWD_fg(freq)), label='Extragalactic WD Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["scarlet"])
ax.plot(BBO_freq, np.sqrt(BBO_freq*BBO_strain), linewidth='1.5', color=sns.xkcd_rgb["dusty purple"], linestyle='--', label='BBO strain noise')
ax.set_xlim([1e-4, 1e1])
ax.set_ylim([1e-25, 1e-19])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10, loc='upper right')
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_6_right.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()


####################################################################################################################
# Plot strain sensitivity for ET (Fig. 7, left panel)

fig = plt.figure()
ax = plt.gca()

freq = np.logspace(np.log10(1e-7), np.log10(1e5), 10000)
ax.plot(freq, np.sqrt(freq*BBH_BNS_fg(freq)), label='BBH+BNS Foreground', linewidth='1.5', linestyle='-.', color=sns.xkcd_rgb["green"])
ax.loglog(ET_freq, np.sqrt(ET_freq*ET_strain), linewidth='1.5', color=sns.xkcd_rgb["denim blue"], linestyle='--', label='ET strain noise')
ax.set_xlim([1e-1, 1e4])
ax.set_ylim([1e-26, 1e-18])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Characteristic Strain $\sqrt{f S_{n}(f)}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=10)
plt.savefig(op.join(op.dirname(__file__), 'figures/Fig_7_left.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
