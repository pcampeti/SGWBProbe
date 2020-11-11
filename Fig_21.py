#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:31:56 2020

@author: Paolo Campeti

This script reproduces Figure 21 in the paper. 
Uses methods imported from module sgwbprobecomb/SGWB_Signal.py, 
sgwbprobecomb/Binned_errors.py. and sgwbprobecomb/error_boxes.py.

"""
import os.path as op
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
# import our classes and methods
from sgwbprobe.SGWB_Signal import Signal_GW
from sgwbprobe.Binned_errors import Binned_GW
from sgwbprobe.error_boxes import make_error_boxes

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

# Constants
year_sec = 60*60*24*365 

# Load and unpack PTA and Interferometers instrumental strains

# Einstein Telescope
ET = np.load(op.join(op.dirname(__file__),'files/S_h_ET.npz'))
ET_freq = ET['x']
ET_strain = ET['y']
eff_ET = 1. # mission efficiency factor
ET_T_obs = 1 * year_sec * eff_ET

# LISA for Cosmologists
LISA_xcosmo = np.load(op.join(op.dirname(__file__),'files/S_h_LISA_xcosmo.npz'))
LISA_xcosmo_freq = LISA_xcosmo['x']
LISA_xcosmo_strain = LISA_xcosmo['y']
eff_LISA = 0.75
LISA_xcosmo_T_obs = 4 * year_sec * eff_LISA


# muAres without fgs
Ares_nofgs = np.load(op.join(op.dirname(__file__),'files/S_h_muAres_nofgs.npz'))
Ares_nofgs_freq = Ares_nofgs['x']
Ares_nofgs_strain = Ares_nofgs['y']
eff_Ares = 1.
Ares_nofgs_T_obs = 10 * year_sec * eff_Ares


# BBO STAR
BBO_STAR = np.load(op.join(op.dirname(__file__),'files/S_h_BBO_STAR.npz'))
BBO_STAR_freq = BBO_STAR['x']
BBO_STAR_strain = BBO_STAR['y']
eff_BBO = 1.
BBO_STAR_T_obs = 10 * year_sec * eff_BBO


# DECIGO
DECIGO = np.load(op.join(op.dirname(__file__),'files/S_h_DECIGO.npz'))
DECIGO_freq = DECIGO['x']
DECIGO_strain = DECIGO['y']
eff_DECIGO = 1.
DECIGO_T_obs = 10 * year_sec * eff_DECIGO


# DO Optimal
DO = np.load(op.join(op.dirname(__file__),'files/S_h_DO_Optimal.npz'))
DO_freq = DO['x']
DO_strain = DO['y']
eff_DO = 0.75
DO_T_obs = 4 * year_sec * eff_DO


# DO Conservative
DO_cons = np.load(op.join(op.dirname(__file__),'files/S_h_DO_Conservative.npz'))
DO_cons_freq = DO_cons['x']
DO_cons_strain = DO_cons['y']
eff_DO = 0.75
DO_cons_T_obs = 4 * year_sec * eff_DO

###############################################################################

# Generate primordial signals and wavenumber vector
 
class_axion1 = Signal_GW(r_vac=1e-5, r_star=400, k_p=1e15, sigma=9.1, axion=True, running=True)
class_axion2 = Signal_GW(r_vac=1e-5, r_star=0.15 , k_p=1e11, sigma=8, axion=True, running=True)
class_no_axion_r001 = Signal_GW(r_vac=0.01, axion=None, running=True)
class_no_axion = Signal_GW(r_vac=0.001, axion=None, running=True)

k = np.logspace(np.log10(1e-5), np.log10(1e16), 100000)

###############################################################################

#class for LISA

sens_curve_LISA = np.array(LISA_xcosmo_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_LISA = np.array(LISA_xcosmo_freq) * 6.5e14

class_binned = Binned_GW(
        name_exp = 'LISA',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=6.5e9,
                         N_bins_sens=1,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         interp=True,
                         sigma_L=1.0
                         )

binned_signal_whole, bins_mean_point_whole = class_binned.Omega_GW_binning()
xerr, yerr, bins_mean_point, binned_signal, binned_curve = class_binned.sens_curve_binning()

################################################################################

# BBO
omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
omega_gw_flat = class_no_axion.analytic_omega_WK(k)
k_BBO = np.array(BBO_STAR_freq) * 6.5e14 
sens_curve_BBO = np.array(BBO_STAR_strain) 

class_BBO = Binned_GW(
        name_exp = 'BBO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_flat,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=4,
                         T_obs=BBO_STAR_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0
                         )

xerr_BBO, yerr_BBO, bins_mean_point_BBO, binned_signal_BBO, binned_curve_BBO = class_BBO.sens_curve_binning()

###############################################################################

#class for LiteBIRD and r=0.001
Fisher = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_4.0_r0001.npy'))  
omega_gw_flat = class_no_axion.analytic_omega_WK(k)
power_spectrum = class_no_axion.tensor_spect(k)

class_binned_flat_CMB = Binned_GW(
        name_exp='LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         omega_gw=omega_gw_flat,
                         kmin_sens=1e-4,
                         N_bins_sens=2,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_flat, bins_mean_point_whole_flat = class_binned_flat_CMB.Omega_GW_binning()
xerr_flat, yerr_flat, bins_mean_point_flat, binned_signal_flat, binned_curve_flat = class_binned_flat_CMB.sens_curve_binning()

###############################################################################

#class for LiteBIRD and r=0.01

Fisher = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_0.9_r001.npy'))  
omega_gw_flat_r001 = class_no_axion_r001.analytic_omega_WK(k)
power_spectrum_r001 = class_no_axion_r001.tensor_spect(k)

class_binned_flat_CMB_r001 = Binned_GW(
        name_exp='LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         omega_gw=omega_gw_flat_r001,
                         kmin_sens=1e-4,
                         N_bins_sens=2,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum_r001,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_flat_r001, bins_mean_point_whole_flat_r001 = class_binned_flat_CMB_r001.Omega_GW_binning()

###############################################################################

# class for LiteBIRD and AX1 model (just to plot the binned signal, not for error bars)
Fisher_axion = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_0.9_r001.npy'))  
power_spectrum_axion = class_axion1.total_spect(k)

class_binned_axion_CMB = Binned_GW(
        name_exp='LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=0.9,
                         omega_gw=omega_gw,
                         kmin_sens=1e-4,
                         N_bins_sens=10,
                         CMB=True,
                         F=Fisher_axion,
                         tensor_spect=power_spectrum_axion,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_axion, bins_mean_point_whole_axion = class_binned_axion_CMB.Omega_GW_binning()

###############################################################################

# BBO Axion 2

omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
k_BBO2 = np.logspace(np.log10(1e-4), np.log10(1e3), 1000) * 6.5e14

class_BBO2 = Binned_GW(
        name_exp='BBO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=2,
                         T_obs=BBO_STAR_T_obs,
                         n_det = 2.,
                         sigma_L=1.0
                         )

binned_signal_axion2, bins_mean_point_axion2 = class_BBO2.Omega_GW_binning()

###############################################################################

#class for DECIGO

sens_curve_DECIGO = np.array(DECIGO_strain)  
k_decigo = class_axion1.freq_k_conv(DECIGO_freq)
class_DECIGO = Binned_GW(
        name_exp='DECIGO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_flat,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=3,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0
                         )

xerr_decigo, yerr_decigo, bins_mean_point_decigo, binned_signal_decigo, binned_curve_decigo = class_DECIGO.sens_curve_binning()
                
################################################################################
#class for muAres without foregrounds

sens_curve_MUARES_nofgs = np.array(Ares_nofgs_strain)
k_muares_nofgs = class_axion1.freq_k_conv(Ares_nofgs_freq)
class_MUARES_nofgs = Binned_GW(
              name_exp='muAres',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw_flat,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=3,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0
                         )

xerr_muares_nofgs, yerr_muares_nofgs, bins_mean_point_muares_nofgs, binned_signal_muares_nofgs, binned_curve_muares_nofgs = class_MUARES_nofgs.sens_curve_binning()
                
###############################################################################

#class for DO Optimal

sens_curve_DO = np.array(DO_strain)
k_DO = class_axion1.freq_k_conv(DO_freq)
class_DO = Binned_GW(      name_exp='DO_Opt',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_DO,
                         omega_gw=omega_gw,
                         k_sens=k_DO,
                         kmin_sens=k_DO[0],
                         N_bins_sens=7,
                         T_obs=DO_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0
                         )

xerr_DO, yerr_DO, bins_mean_point_DO, binned_signal_DO, binned_curve_DO = class_DO.sens_curve_binning()

###############################################################################

#class for DO Conservative

sens_curve_DO_cons = np.array(DO_cons_strain)
k_DO_cons = class_axion1.freq_k_conv(DO_cons_freq)
class_DO_cons = Binned_GW(
        name_exp='DO_Cons',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_DO_cons,
                         omega_gw=omega_gw,
                         k_sens=k_DO_cons,
                         kmin_sens=k_DO_cons[0],
                         N_bins_sens=7,
                         T_obs=DO_cons_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0
                         )

xerr_DO_cons, yerr_DO_cons, bins_mean_point_DO_cons, binned_signal_DO_cons, binned_curve_DO_cons = class_DO_cons.sens_curve_binning()

###############################################################################
# FOREGROUNDS BELOW

# class for muAres with foregrounds

class_MUARES = Binned_GW(
        name_exp='muAres_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw_flat,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=3,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1.0,
                         )

xerr_muares, yerr_muares, bins_mean_point_muares, binned_signal_muares, binned_curve_muares = class_MUARES.sens_curve_binning()

################################################################################
# BBO with fgs
class_BBO_fgs = Binned_GW(
        name_exp='BBO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=4.0,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_flat,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=4,
                         T_obs=BBO_STAR_T_obs,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1.0,
                         )

xerr_BBO_fgs, yerr_BBO_fgs, bins_mean_point_BBO_fgs, binned_signal_BBO_fgs, binned_curve_BBO_fgs = class_BBO_fgs.sens_curve_binning()

###############################################################################

#PLOT
fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=2, alpha=0.55, linestyle='--')

#plot for LiteBIRD flat spectrum r=0.001
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=6)
_ = make_error_boxes(ax, np.array(bins_mean_point_flat)/6.5e14, binned_signal_flat, xerr_flat/6.5e14, yerr_flat, facecolor='g', alpha=0.55, zorder=5)

#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat_r001)/6.5e14, binned_signal_whole_flat_r001, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=4, linestyle='--')

plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=1, label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$', 
           alpha=0.55, linestyle='--')

_ = make_error_boxes(ax, np.array(bins_mean_point_BBO)/6.5e14, binned_signal_BBO, xerr_BBO/6.5e14, yerr_BBO, facecolor='red', alpha=0.7, zorder=5)

plt.axvline(x=0.2, color='grey', linestyle='--', linewidth=1, zorder=45)

plt.text(5e-19, 1e-13, r'$\bf LiteBIRD$', fontsize=10, color='green', zorder=89)
plt.text(4e-3, 5e-15, r'$\bf BBO$', fontsize=10, color='red', zorder=89)

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()

ax = 1e-21
bx = 1e2
ay = 1e-21
by = 1e-8

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig(op.join(op.dirname(__file__),'figures/Fig_21.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()