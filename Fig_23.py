#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:08:54 2020

@author: Paolo Campeti

Script to reproduce Fig. (23) of the paper.

"""
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns

# seaborn settings
sns.set()
sns.set(style='whitegrid')

# Useful constant
C = 299792458. #speed of light in m/s

# loading response functions

# LISA
Resp_LISA = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/Resp_LISA.npy')
freq_LISA = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/f_R_LISA.npy')
Resp_LISA = Resp_LISA/Resp_LISA[0] # normalize response to 1 for plotting purpose
Resp_LISA = np.abs(Resp_LISA) # take absolute value for plotting purpose

# DECIGO
Resp_DECIGO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/Resp_DECIGO.npy')
freq_DECIGO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/f_R_DECIGO.npy')
Resp_DECIGO = Resp_DECIGO/Resp_DECIGO[0] 
Resp_DECIGO = np.abs(Resp_DECIGO)

# DO 
Resp_DO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/Resp_DO.npy')
freq_DO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/f_R_DO.npy')
Resp_DO = Resp_DO/Resp_DO[0] 
Resp_DO = np.abs(Resp_DO) 

# muAres
Resp_Ares = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/Resp_muAres.npy')
freq_Ares = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/f_R_Ares.npy')
Resp_Ares = Resp_Ares/Resp_Ares[0] 
Resp_Ares = np.abs(Resp_Ares) 

#BBO
Resp_BBO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/Resp_BBO.npy')
freq_BBO = np.load('/home/paolo/Codes/GW_plot_fgs/FINAL_Interferometers/Responses/f_R_BBO.npy')
Resp_BBO = Resp_BBO/Resp_BBO[0] 
Resp_BBO = np.abs(Resp_BBO) 

# Plot Response Functions
fig = plt.figure()
ax = plt.gca()
ax.semilogx(freq_LISA, Resp_LISA, label='LISA', linewidth='2.0')
ax.plot(freq_DECIGO, Resp_DECIGO, label='DECIGO', linewidth='2.0')
ax.plot(freq_DO, Resp_DO, label='DO', linewidth='2.0')
ax.plot(freq_Ares, Resp_Ares, label=r'$\mu$Ares', linewidth='2.0')
ax.plot(freq_BBO, Resp_BBO, label='BBO', linewidth='2.0', color='black')
ax.set_xlim([1e-5, 5e3])
ax.set_ylim([-0.1, 1.1])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'Normalized Overlap Reduction Function $|\mathcal{R}_{IJ}|$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=12, loc='upper right', bbox_to_anchor=(1.3, 0.8))
plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_23.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
