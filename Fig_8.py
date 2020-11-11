#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:31:56 2020

@author: Paolo Campeti

This script reproduces Figure 8 in the paper. 
Uses methods imported from module sgwbprobecomb/SGWB_Signal.py and 
sgwbprobecomb/Binned_errors.py.

"""
import os.path as op
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
# import our classes and methods
from sgwbprobecomb.SGWB_Signal import Signal_GW
from sgwbprobecomb.Binned_errors import Binned_GW

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

# Useful constant: seconds in a year
year_sec = 60*60*24*365 

# Load and unpack PTA and Interferometers instrumental strains

# SKA
SKA_file = np.load(op.join(op.dirname(__file__), 'files/hc_SKA.npz'))
SKA_freq = SKA_file['x']
SKA_hc = SKA_file['y']
SKA_strain = SKA_hc**2/SKA_freq
eff_SKA = 1. # mission efficiency factor
SKA_T_obs = 10 * year_sec * eff_SKA


# Einstein Telescope
ET = np.load(op.join(op.dirname(__file__),'files/S_h_ET.npz'))
ET_freq = ET['x']
ET_strain = ET['y']
eff_ET = 1. # mission efficiency factor
ET_T_obs = 1 * year_sec * eff_ET

# Advanced LIGO
aLIGO = np.load(op.join(op.dirname(__file__),'files/S_h_aLIGO.npz'))
aLIGO_freq = aLIGO['x']
aLIGO_strain = aLIGO['y']
eff_aLIGO = 1. # mission efficiency factor
aLIGO_T_obs = 4 * year_sec * eff_aLIGO


# LISA 
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

# AEDGE
AEDGE = np.load(op.join(op.dirname(__file__),'files/S_h_AEDGE.npz'))
AEDGE_freq = AEDGE['x']
AEDGE_strain = AEDGE['y']
eff_AEDGE = 0.6
AEDGE_T_obs = 5 * year_sec * eff_AEDGE

###############################################################################
# Generate primordial signals and wavenumber vector
 
class_axion1 = Signal_GW(r_vac=1e-5, r_star=835, k_p=1e13, sigma=9, axion=True)
class_axion2 = Signal_GW(r_vac=1e-5, r_star=0.15 , k_p=1e11, sigma=8, axion=True)
class_no_axion = Signal_GW(r_vac=0.01, axion=None)
class_no_axion_r0001 = Signal_GW(r_vac=0.001, axion=None)

# Wavenumber array
k = np.logspace(np.log10(1e-5), np.log10(1e20), 100000)

###############################################################################

#class for SKA

sens_curve_SKA = np.array(SKA_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_SKA = np.array(SKA_freq) * 6.5e14

class_binned_SKA = Binned_GW(
        name_exp='SKA',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_SKA,
                         omega_gw=omega_gw,
                         k_sens=k_SKA,
                         kmin_sens=k_SKA[0],
                         N_bins_sens=5,
                         T_obs=SKA_T_obs,
                         n_det=1.,
                         sigma_L=1.0,
                         )

xerr_SKA, yerr_SKA, bins_mean_point_SKA, binned_signal_SKA, binned_curve_SKA = class_binned_SKA.sens_curve_binning()
###############################################################################

#class for AEDGE

sens_curve_AEDGE = np.array(AEDGE_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_AEDGE = np.array(AEDGE_freq) * 6.5e14

class_binned_AEDGE = Binned_GW(name_exp='AEDGE',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_AEDGE,
                         omega_gw=omega_gw,
                         k_sens=k_AEDGE,
                         kmin_sens=k_AEDGE[0],
                         N_bins_sens=4,
                         T_obs=AEDGE_T_obs,
                         n_det = 1.,
                         interp=True,
                         sigma_L=1.0,
                         )

xerr_AEDGE, yerr_AEDGE, bins_mean_point_AEDGE, binned_signal_AEDGE, binned_curve_AEDGE = class_binned_AEDGE.sens_curve_binning()

###############################################################################

#class for ET

sens_curve_ET = np.array(ET_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_ET = np.array(ET_freq) * 6.5e14

class_binned_ET = Binned_GW(name_exp='ET',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_ET,
                         omega_gw=omega_gw,
                         k_sens=k_ET,
                         kmin_sens=1.5*6.5e14 ,
                         N_bins_sens=5,
                         T_obs=ET_T_obs,
                         interp=True,
                         n_det = 3.,
                         sigma_L=1.0,
                         )

xerr_ET, yerr_ET, bins_mean_point_ET, binned_signal_ET, binned_curve_ET = class_binned_ET.sens_curve_binning()

###############################################################################
#class for aLIGO

sens_curve_aLIGO = np.array(aLIGO_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_aLIGO = np.array(aLIGO_freq) * 6.5e14

class_binned_aLIGO = Binned_GW(name_exp='Adv_LIGO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_aLIGO,
                         omega_gw=omega_gw,
                         k_sens=k_aLIGO,
                         kmin_sens=k_aLIGO[0],
                         N_bins_sens=5,
                         T_obs=aLIGO_T_obs,
                         n_det = 1.,
                         sigma_L=1.,
                         )

xerr_aLIGO, yerr_aLIGO, bins_mean_point_aLIGO, binned_signal_aLIGO, binned_curve_aLIGO = class_binned_aLIGO.sens_curve_binning()
###############################################################################
#class for LISA

sens_curve_LISA = np.array(LISA_xcosmo_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_LISA = np.array(LISA_xcosmo_freq) * 6.5e14

class_binned = Binned_GW(name_exp='LISA',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=8,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         interp=True,
                         sigma_L=1.0,
                         )

binned_signal_whole, bins_mean_point_whole = class_binned.Omega_GW_binning()
xerr, yerr, bins_mean_point, binned_signal, binned_curve = class_binned.sens_curve_binning()

################################################################################

# BBO

omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
k_BBO = np.array(BBO_STAR_freq) * 6.5e14 
sens_curve_BBO = np.array(BBO_STAR_strain) 

class_BBO = Binned_GW(name_exp='BBO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=12,
                         T_obs=BBO_STAR_T_obs,
                         n_det = 2.,
                         interp=True,
                         sigma_L=1.0,
                         )

binned_signal_axion2, bins_mean_point_axion2 = class_BBO.Omega_GW_binning()
xerr_BBO, yerr_BBO, bins_mean_point_BBO, binned_signal_BBO, binned_curve_BBO = class_BBO.sens_curve_binning()

###############################################################################
#class for LiteBIRD and r=0

Fisher = np.load(op.join(op.dirname(__file__),'files/LiteBIRD_Fisher_matrices/Fisher_1.3_r0.npy'))  
omega_gw_flat_r0001 = class_no_axion_r0001.analytic_omega_WK(k)
power_spectrum_r0001 = class_no_axion_r0001.tensor_spect(k)

class_binned_flat_CMB_0001 = Binned_GW(name_exp='LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         omega_gw=omega_gw_flat_r0001,
                         kmin_sens=1e-4,
                         N_bins_sens=7,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum_r0001,
                         sigma_L=1.0,
                         )                 
                      
xerr_flat_0001, yerr_flat_0001, bins_mean_point_flat_0001, binned_signal_flat_0001, binned_curve_flat_0001 = class_binned_flat_CMB_0001.sens_curve_binning()

################################################################################

#class for DECIGO

sens_curve_DECIGO = np.array(DECIGO_strain)  
k_decigo = class_axion1.freq_k_conv(DECIGO_freq)
class_DECIGO = Binned_GW(name_exp='DECIGO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=11,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0,
                         )

xerr_decigo, yerr_decigo, bins_mean_point_decigo, binned_signal_decigo, binned_curve_decigo = class_DECIGO.sens_curve_binning()
                
################################################################################

#class for muAres without foregrounds

sens_curve_MUARES_nofgs = np.array(Ares_nofgs_strain)
k_muares_nofgs = class_axion1.freq_k_conv(Ares_nofgs_freq)
class_MUARES_nofgs = Binned_GW(name_exp='muAres',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=12,
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
class_DO = Binned_GW(name_exp='DO_Opt',
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
                         sigma_L=1.0,
                         )

xerr_DO, yerr_DO, bins_mean_point_DO, binned_signal_DO, binned_curve_DO = class_DO.sens_curve_binning()

###############################################################################

#class for DO Conservative

sens_curve_DO_cons = np.array(DO_cons_strain)
k_DO_cons = class_axion1.freq_k_conv(DO_cons_freq)
class_DO_cons = Binned_GW(name_exp='DO_Cons',
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
                         sigma_L=1.0,
                         )

xerr_DO_cons, yerr_DO_cons, bins_mean_point_DO_cons, binned_signal_DO_cons, binned_curve_DO_cons = class_DO_cons.sens_curve_binning()


###############################################################################
###############################################################################
##########       FOREGROUNDS BELOW!!                                ###########
###############################################################################
###############################################################################

#class for SKA for fgs

class_binned_SKA_fgs = Binned_GW(name_exp='SKA_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_SKA,
                         omega_gw=omega_gw,
                         k_sens=k_SKA,
                         kmin_sens=k_SKA[0],
                         N_bins_sens=5,
                         T_obs=SKA_T_obs,
                         n_det=1.,
                         fgs=True,
                         sigma_L=1.0,
                         )

xerr_SKA_fgs, yerr_SKA_fgs, bins_mean_point_SKA_fgs, binned_signal_SKA_fgs, binned_curve_SKA_fgs = class_binned_SKA_fgs.sens_curve_binning()

###############################################################################
# Class for ET with fgs

class_binned_ET_fgs = Binned_GW(name_exp='ET_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_ET,
                         omega_gw=omega_gw,
                         k_sens=k_ET,
                         kmin_sens=1.5*6.5e14 ,
                         N_bins_sens=5,
                         T_obs=ET_T_obs,
                         n_det = 3.,
                         fgs=True,
                         interp=True,
                         sigma_L=1.0,
                         )

xerr_ET_fgs, yerr_ET_fgs, bins_mean_point_ET_fgs, binned_signal_ET_fgs, binned_curve_ET_fgs = class_binned_ET_fgs.sens_curve_binning()


###############################################################################

#class for AEDGE with fgs

class_binned_AEDGE_fgs = Binned_GW(name_exp='AEDGE_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_AEDGE,
                         omega_gw=omega_gw,
                         k_sens=k_AEDGE,
                         kmin_sens=k_AEDGE[0],
                         N_bins_sens=4,
                         T_obs=AEDGE_T_obs,
                         n_det = 1.,
                         interp=True,
                         fgs=True,
                         sigma_L=0.1,
                         )

xerr_AEDGE_fgs, yerr_AEDGE_fgs, bins_mean_point_AEDGE_fgs, binned_signal_AEDGE_fgs, binned_curve_AEDGE_fgs = class_binned_AEDGE_fgs.sens_curve_binning()

###############################################################################

#class for LISA with fgs

class_binned_fgs = Binned_GW(name_exp='LISA_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=8,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         interp=True,
                         fgs=True,
                         sigma_L=0.1,
                         )

xerr_fgs, yerr_fgs, bins_mean_point_fgs, binned_signal_fgs, binned_curve_fgs = class_binned_fgs.sens_curve_binning()

################################################################################
#class for DECIGO with fgs

class_DECIGO_fgs = Binned_GW(name_exp='DECIGO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=11,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1e-3,
                         )

xerr_decigo_fgs, yerr_decigo_fgs, bins_mean_point_decigo_fgs, binned_signal_decigo_fgs, binned_curve_decigo_fgs = class_DECIGO_fgs.sens_curve_binning()
##############################################################################

#class for DECIGO only spectral shape

class_DECIGO_spectral = Binned_GW(name_exp='DECIGO_spectral_shape',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=11,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1.0,
                         )

xerr_decigo_spectral, yerr_decigo_spectral, bins_mean_point_decigo_spectral, binned_signal_decigo_spectral, binned_curve_decigo_spectral = class_DECIGO_spectral.sens_curve_binning()
                
###############################################################################

#class for DO Optimal with fgs

class_DO_fgs = Binned_GW(name_exp='DO_Optimal_with_fgs',
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
                         fgs=True,
                         sigma_L=0.1,
                         )

xerr_DO_fgs, yerr_DO_fgs, bins_mean_point_DO_fgs, binned_signal_DO_fgs, binned_curve_DO_fgs = class_DO_fgs.sens_curve_binning()

###############################################################################

#class for DO Conservative with fgs
class_DO_cons_fgs = Binned_GW(name_exp='DO_Cons_withfgs',
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
                         fgs=True,
                         sigma_L=0.1,
                         )

xerr_DO_cons_fgs, yerr_DO_cons_fgs, bins_mean_point_DO_cons_fgs, binned_signal_DO_cons_fgs, binned_curve_DO_cons_fgs = class_DO_cons_fgs.sens_curve_binning()

###############################################################################
# class for muAres with foregrounds
class_MUARES = Binned_GW(name_exp='muAres_two_fgs_spectral',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.3,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=12,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1e-3,
                         fgs=True,
                         )

xerr_muares, yerr_muares, bins_mean_point_muares, binned_signal_muares, binned_curve_muares = class_MUARES.sens_curve_binning()
                
###############################################################################
# Plot
fig = plt.figure()
ax = plt.gca()
ax.loglog(np.array(bins_mean_point_BBO)/6.5e14, binned_curve_BBO[:len(bins_mean_point_BBO)], linewidth='1.5', color=sns.xkcd_rgb["dusty purple"], linestyle='--')
ax.loglog(np.array(bins_mean_point_decigo)/6.5e14, binned_curve_decigo[:len(bins_mean_point_decigo)], linewidth='1.5', color=sns.xkcd_rgb["amber"], linestyle='--')
ax.loglog(np.array(bins_mean_point_aLIGO)/6.5e14, binned_curve_aLIGO[:len(bins_mean_point_aLIGO)], label='aLIGO', linewidth='1.5', color=sns.xkcd_rgb["yellow"], linestyle='--')
ax.loglog(np.array(bins_mean_point_DO)/6.5e14, binned_curve_DO[:len(bins_mean_point_DO)], linewidth='1.5', color=sns.xkcd_rgb["greyish"], linestyle='--')
ax.loglog(np.array(bins_mean_point_DO_cons)/6.5e14, binned_curve_DO_cons[:len(bins_mean_point_DO_cons)], linewidth='1.5', color=sns.xkcd_rgb["faded green"], linestyle='--')
ax.loglog(np.array(bins_mean_point_ET)/6.5e14, binned_curve_ET[:len(bins_mean_point_ET)], linewidth='1.5', color=sns.xkcd_rgb["steel blue"], linestyle='--')
ax.loglog(np.array(bins_mean_point_muares_nofgs)/6.5e14, binned_curve_muares_nofgs[:len(bins_mean_point_muares_nofgs)], linewidth='1.5', color=sns.xkcd_rgb["cyan"], linestyle='--')
ax.loglog(np.array(bins_mean_point)/6.5e14, binned_curve[:len(bins_mean_point)], linewidth='1.5', color=sns.xkcd_rgb["black"], linestyle='--')
ax.loglog(np.array(bins_mean_point_flat_0001)/6.5e14, binned_curve_flat_0001[:len(bins_mean_point_flat_0001)], label='LiteBIRD r=0', linewidth='1.5', color=sns.xkcd_rgb["scarlet"], linestyle='-')
ax.loglog(np.array(bins_mean_point_AEDGE)/6.5e14, binned_curve_AEDGE[:len(bins_mean_point_AEDGE)], linewidth='1.5', color=sns.xkcd_rgb["pale red"], linestyle='--')
ax.loglog(np.array(bins_mean_point_SKA)/6.5e14, binned_curve_SKA[:len(bins_mean_point_SKA)], linewidth='1.5', color=sns.xkcd_rgb["violet"], linestyle='--')

ax.loglog(np.array(bins_mean_point_decigo_fgs)/6.5e14, binned_curve_decigo_fgs[:len(bins_mean_point_decigo_fgs)], label='DECIGO', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["amber"])
ax.loglog(np.array(bins_mean_point_decigo_spectral)/6.5e14, binned_curve_decigo_spectral[:len(bins_mean_point_decigo_spectral)], linestyle='-.', linewidth='1.5', color=sns.xkcd_rgb["amber"])
ax.loglog(np.array(bins_mean_point_DO_fgs)/6.5e14, binned_curve_DO_fgs[:len(bins_mean_point_DO_fgs)], label='DO Optimal', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["greyish"])
ax.loglog(np.array(bins_mean_point_DO_cons_fgs)/6.5e14, binned_curve_DO_cons_fgs[:len(bins_mean_point_DO_cons_fgs)],label='DO Conservative', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["faded green"])
ax.loglog(np.array(bins_mean_point_muares)/6.5e14, binned_curve_muares[:len(bins_mean_point_muares)], label=r'$\mu$Ares', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["cyan"])
ax.loglog(np.array(bins_mean_point_fgs)/6.5e14, binned_curve_fgs[:len(bins_mean_point_fgs)], label='LISA', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["black"])
ax.loglog(np.array(bins_mean_point_AEDGE_fgs)/6.5e14, binned_curve_AEDGE_fgs[:len(bins_mean_point_AEDGE_fgs)], label='AEDGE', linestyle='-', linewidth='1.5', color=sns.xkcd_rgb["pale red"])
ax.loglog(np.array(bins_mean_point_ET_fgs)/6.5e14, binned_curve_ET_fgs[:len(bins_mean_point_ET_fgs)], label='ET', linewidth='1.5', linestyle='-', color=sns.xkcd_rgb["steel blue"])
ax.loglog(np.array(bins_mean_point_SKA_fgs)/6.5e14, binned_curve_SKA_fgs[:len(bins_mean_point_SKA_fgs)], label='SKA', linewidth='1.5', linestyle='-', color=sns.xkcd_rgb["violet"])

ax.set_xlim([1e-19, 1e4])
ax.set_ylim([1e-19, 1e-6])
plt.xlabel(r'f $[Hz]$',fontsize = 10.0)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = 10.0)
plt.tick_params(axis = 'both',which = 'major', labelsize = 10.0)
ax.legend(fontsize=8, bbox_to_anchor=(1, 1.0))
plt.savefig(op.join(op.dirname(__file__),'figures/Fig_8.pdf'), format='pdf', dpi=1000, bbox_inches='tight')
plt.show()