#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:31:56 2020

@author: Paolo Campeti

This script reproduces Figure 14-16 in the paper. 
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

#seaborn settings
sns.set()
sns.set(style='whitegrid')

#matplotlib settings
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
Ares_T_obs = Ares_nofgs_T_obs  

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


# AEDGE
AEDGE = np.load(op.join(op.dirname(__file__),'files/S_h_AEDGE.npz'))
AEDGE_freq = AEDGE['x']
AEDGE_strain = AEDGE['y']
eff_AEDGE = 1.#0.6
AEDGE_T_obs = 10 * year_sec * eff_AEDGE

###############################################################################

# Generate primordial signals and wavenumber vector
class_axion1 = Signal_GW(r_vac=1e-5, r_star=400, k_p=1e15, sigma=9.1, axion=True, running=True)
class_axion2 = Signal_GW(r_vac=1e-5, r_star=0.15 , k_p=1e11, sigma=8, axion=True, running=True)
class_no_axion = Signal_GW(r_vac=0.01, axion=None, running=True)
class_no_axion_r0001 = Signal_GW(r_vac=0.001, axion=None, running=True)

# wavenumber array
k = np.logspace(np.log10(1e-5), np.log10(1e16), 100000)

###############################################################################

#class for AEDGE

sens_curve_AEDGE = np.array(AEDGE_strain)
omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
k_AEDGE = np.array(AEDGE_freq) * 6.5e14

class_binned_AEDGE = Binned_GW(
                         name_exp='AEDGE',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_AEDGE,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_AEDGE,
                         kmin_sens=k_AEDGE[0],
                         N_bins_sens=3,
                         T_obs=AEDGE_T_obs,
                         n_det = 1.,
                         interp=True,
                         sigma_L=0.1,
                         )

xerr_AEDGE, yerr_AEDGE, bins_mean_point_AEDGE, binned_signal_AEDGE, binned_curve_AEDGE = class_binned_AEDGE.sens_curve_binning()

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
                         delta_log_k=2.,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=5,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         interp=True,
                         fgs=True,
                         sigma_L=0.1,
                         )

binned_signal_whole, bins_mean_point_whole = class_binned.Omega_GW_binning()
xerr, yerr, bins_mean_point, binned_signal, binned_curve = class_binned.sens_curve_binning()
################################################################################

# BBO
omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
k_BBO = np.array(BBO_STAR_freq) * 6.5e14 
sens_curve_BBO = np.array(BBO_STAR_strain) 

class_BBO = Binned_GW(
         name_exp = 'BBO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=8,
                         T_obs=BBO_STAR_T_obs,
                         n_det = 2.,
                         interp=True,
                         sigma_L=1.0,
                         )

binned_signal_axion2, bins_mean_point_axion2 = class_BBO.Omega_GW_binning()
xerr_BBO, yerr_BBO, bins_mean_point_BBO, binned_signal_BBO, binned_curve_BBO = class_BBO.sens_curve_binning()

###############################################################################

#class for LiteBIRD and r=0.01 (used only to plot the r=0.01 binned signal,
#not for error bars)

Fisher = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_0.9_r001.npy'))  
omega_gw_flat = class_no_axion.analytic_omega_WK(k)
power_spectrum = class_no_axion.tensor_spect(k)

class_binned_flat_CMB = Binned_GW(
         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         omega_gw=omega_gw_flat,
                         kmin_sens=1e-4,
                         N_bins_sens=3,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum,
                         sigma_L=1.0,
                         )                 
                      
binned_signal_whole_flat, bins_mean_point_whole_flat = class_binned_flat_CMB.Omega_GW_binning()

###############################################################################

# class for LiteBIRD and Axion model r_vac=1e-5, r_star=835, k_p=1e13, sigma=9 
#(used only to plot the r=0.01 binned signal, not for error bars)

Fisher_axion = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_0.9_r001.npy'))  
power_spectrum_axion = class_axion1.total_spect(k)

class_binned_axion_CMB = Binned_GW(
         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         omega_gw=omega_gw,
                         kmin_sens=1e-4,
                         N_bins_sens=3,
                         CMB=True,
                         F=Fisher_axion,
                         tensor_spect=power_spectrum_axion,
                         sigma_L=1.0,
                         )                 
                      
binned_signal_whole_axion, bins_mean_point_whole_axion = class_binned_axion_CMB.Omega_GW_binning()

###############################################################################
# class for LiteBIRD and Axion model 2 
Fisher_axion2 = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_2.0_AX2.npy'))  
power_spectrum_axion2 = class_axion2.total_spect(k)

class_binned_axion2_CMB = Binned_GW(
         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         omega_gw=omega_gw_axion2,
                         kmin_sens=1e-4,
                         N_bins_sens=3,
                         CMB=True,
                         F=Fisher_axion2,
                         tensor_spect=power_spectrum_axion2,
                         sigma_L=1.0,
                         )                 
                      
binned_signal_whole_axion2, bins_mean_point_whole_axion2 = class_binned_axion2_CMB.Omega_GW_binning()
xerr_axion2, yerr_axion2, bins_mean_point_axion2_LiteBIRD, binned_signal_axion2_LiteBIRD, binned_curve_axion2_LiteBIRD = class_binned_axion2_CMB.sens_curve_binning()

###############################################################################

#class for DECIGO
sens_curve_DECIGO = np.array(DECIGO_strain)  
k_decigo = class_axion1.freq_k_conv(DECIGO_freq)
class_DECIGO = Binned_GW(
                         name_exp = 'DECIGO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2., 
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=6,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0,
                         )

xerr_decigo, yerr_decigo, bins_mean_point_decigo, binned_signal_decigo, binned_curve_decigo = class_DECIGO.sens_curve_binning()
                                
###############################################################################

#class for muAres without foregrounds

sens_curve_MUARES_nofgs = np.array(Ares_nofgs_strain)
k_muares_nofgs = class_axion1.freq_k_conv(Ares_nofgs_freq)
class_MUARES_nofgs = Binned_GW(
                         name_exp = 'muAres',                         
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=7,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0,
                         )

xerr_muares_nofgs, yerr_muares_nofgs, bins_mean_point_muares_nofgs, binned_signal_muares_nofgs, binned_curve_muares_nofgs = class_MUARES_nofgs.sens_curve_binning()
                
###############################################################################

#class for DO Optimal

sens_curve_DO = np.array(DO_strain)
k_DO = class_axion1.freq_k_conv(DO_freq)
class_DO = Binned_GW(name_exp = 'DO_Optimal',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_DO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_DO,
                         kmin_sens=k_DO[0],
                         N_bins_sens=3,
                         T_obs=DO_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0,
                         )

xerr_DO, yerr_DO, bins_mean_point_DO, binned_signal_DO, binned_curve_DO = class_DO.sens_curve_binning()

###############################################################################

#class for DO Conservative

sens_curve_DO_cons = np.array(DO_cons_strain)
k_DO_cons = class_axion1.freq_k_conv(DO_cons_freq)
class_DO_cons = Binned_GW(
        name_exp = 'DO_Conservative',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_DO_cons,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_DO_cons,
                         kmin_sens=k_DO_cons[0],
                         N_bins_sens=3,
                         T_obs=DO_cons_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0,
                         )

xerr_DO_cons, yerr_DO_cons, bins_mean_point_DO_cons, binned_signal_DO_cons, binned_curve_DO_cons = class_DO_cons.sens_curve_binning()

###############################################################################

# class for LiteBIRD and r=0.001 (used only to plot the r=0.001 binned signal,
# not for error bars)
Fisher = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_0.9_r001.npy'))  
omega_gw_flat_r0001 = class_no_axion_r0001.analytic_omega_WK(k)
power_spectrum_r0001 = class_no_axion_r0001.tensor_spect(k)

class_binned_flat_CMB_r0001 = Binned_GW(
         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.0,
                         omega_gw=omega_gw_flat_r0001,
                         kmin_sens=1e-4,
                         N_bins_sens=1,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum_r0001,
                         sigma_L=1.0,
                         )                 
                      
binned_signal_whole_flat_r0001, bins_mean_point_whole_flat_r0001 = class_binned_flat_CMB_r0001.Omega_GW_binning()

###############################################################################
#FOREGROUNDS BELOW

# BBO with fgs

class_BBO_fgs = Binned_GW(
                         name_exp='BBO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.0,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=8,
                         T_obs=BBO_STAR_T_obs,
                         n_det = 2.,
                         interp=True,
                         fgs=True,
                         sigma_L=1e-9,
                         )

xerr_BBO_fgs, yerr_BBO_fgs, bins_mean_point_BBO_fgs, binned_signal_BBO_fgs, binned_curve_BBO_fgs = class_BBO_fgs.sens_curve_binning()

################################################################################
#class for DECIGO with fgs
class_DECIGO_fgs = Binned_GW(
                         name_exp = 'DECIGO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.0,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=6,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1e-3,
                         )

xerr_decigo_fgs, yerr_decigo_fgs, bins_mean_point_decigo_fgs, binned_signal_decigo_fgs, binned_curve_decigo_fgs = class_DECIGO_fgs.sens_curve_binning()

###############################################################################
# class for muAres with foregrounds

class_MUARES = Binned_GW(
        name_exp = 'muAres_two_fgs_spectral',                        
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=7,
                         T_obs=Ares_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1e-3,
                         )

xerr_muares, yerr_muares, bins_mean_point_muares, binned_signal_muares, binned_curve_muares = class_MUARES.sens_curve_binning()
###############################################################################

#class for DECIGO with fgs
class_DECIGO_spectral = Binned_GW(
                         name_exp = 'DECIGO_spectral_shape',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=2.0,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=6,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1.,
                         )

xerr_decigo_spectral, yerr_decigo_spectral, bins_mean_point_decigo_spectral, binned_signal_decigo_spectral, binned_curve_decigo_spectral = class_DECIGO_spectral.sens_curve_binning()

###############################################################################


#PLOT
fig = plt.figure()
ax = plt.gca()


#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18, alpha=0.55, linestyle='--')

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')

#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=19, alpha=0.55, linestyle='--')

#plot for LiteBIRD axion2 spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion2_LiteBIRD)/6.5e14, binned_signal_axion2_LiteBIRD, xerr_axion2/6.5e14, yerr_axion2, facecolor='g', alpha=0.55, zorder=4)

#plot for muAres without fgs axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_muares_nofgs)/6.5e14, binned_signal_muares_nofgs, xerr_muares_nofgs/6.5e14, yerr_muares_nofgs, facecolor=sns.xkcd_rgb["amber"], alpha=0.7, zorder=4)
_ = make_error_boxes(ax, np.array(bins_mean_point_muares)/6.5e14, binned_signal_muares, xerr_muares/6.5e14, yerr_muares, facecolor=sns.xkcd_rgb["amber"], alpha=0.3, zorder=3)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=20, label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$')


plt.text(9e-17, 1.5e-16, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(7e-6, 1e-12, r'$\bf \mu Ares$', fontsize=11, color=sns.xkcd_rgb["amber"])

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()

ax = 1e-21
bx = 1e2
ay = 1e-19
by = 1e-10

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)



plt.savefig(op.join(op.dirname(__file__),'figures/Fig_14.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################

#PLOT
fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18, alpha=0.55, linestyle='--')

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')


#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=19, alpha=0.55, linestyle='--')

#plot for LiteBIRD axion2 spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion2_LiteBIRD)/6.5e14, binned_signal_axion2_LiteBIRD, xerr_axion2/6.5e14, yerr_axion2, facecolor='g', alpha=0.55, zorder=4)

#plot for DECIGO axion spectrum 
_ = make_error_boxes(ax, np.array(bins_mean_point_decigo)/6.5e14, binned_signal_decigo, xerr_decigo/6.5e14, yerr_decigo, facecolor='blue', alpha=0.7, zorder=4)
_ = make_error_boxes(ax, np.array(bins_mean_point_decigo_fgs)/6.5e14, binned_signal_decigo_fgs, xerr_decigo_fgs/6.5e14, yerr_decigo_fgs, facecolor='blue', alpha=0.4, zorder=4)
_ = make_error_boxes(ax, np.array(bins_mean_point_decigo_spectral)/6.5e14, binned_signal_decigo_spectral, xerr_decigo_spectral/6.5e14, yerr_decigo_spectral, facecolor='blue', alpha=0.2, zorder=2)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=20, label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$')

plt.axvline(x=0.2, color='grey', linestyle='--', linewidth=1, zorder=45)


plt.text(9e-17, 1.5e-16, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(3e-3, 9e-15, r'$\bf DECIGO$', fontsize=10, color='blue')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()

ax = 1e-21
bx = 1e2
ay = 1e-19
by = 1e-10

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig(op.join(op.dirname(__file__),'figures/Fig_15.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
#PLOT
fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18, alpha=0.55, linestyle='--')

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')


#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=19, alpha=0.55, linestyle='--')

#plot for LiteBIRD axion2 spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion2_LiteBIRD)/6.5e14, binned_signal_axion2_LiteBIRD, xerr_axion2/6.5e14, yerr_axion2, facecolor='g', alpha=0.55, zorder=4)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=20, label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$')
_ = make_error_boxes(ax, np.array(bins_mean_point_BBO)/6.5e14, binned_signal_BBO, xerr_BBO/6.5e14, yerr_BBO, facecolor='red', alpha=0.7, zorder=4)

plt.axvline(x=0.2, color='grey', linestyle='--', linewidth=1, zorder=45)


plt.text(9e-17, 1.5e-16, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(5e-3, 2.4e-15, r'$\bf BBO$', fontsize=10, color='red')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()

ax = 1e-21
bx = 1e2
ay = 1e-19
by = 1e-10

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)
plt.savefig(op.join(op.dirname(__file__),'figures/Fig_16.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()