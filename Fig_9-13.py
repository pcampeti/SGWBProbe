"""
@author: Paolo Campeti

"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
# import our classes and methods
from sgwbprobecomb.SGWB_Signal import Signal_GW
from sgwbprobecomb.Binned_errors import Binned_GW
from sgwbprobecomb.error_boxes import make_error_boxes

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
SKA_file = np.load('/home/paolo/Codes/SGWBProbeComb/files/hc_SKA.npz')
SKA_freq = SKA_file['x']
SKA_hc = SKA_file['y']
SKA_strain = SKA_hc**2/SKA_freq
eff_SKA = 1. # mission efficiency factor
SKA_T_obs = 10 * year_sec * eff_SKA


# Einstein Telescope
ET = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_ET.npz')
ET_freq = ET['x']
ET_strain = ET['y']
eff_ET = 1. # mission efficiency factor
ET_T_obs = 1 * year_sec * eff_ET


# LISA for Cosmologists
LISA_xcosmo = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_LISA_xcosmo.npz')
LISA_xcosmo_freq = LISA_xcosmo['x']
LISA_xcosmo_strain = LISA_xcosmo['y']
eff_LISA = 0.75
LISA_xcosmo_T_obs = 4 * year_sec * eff_LISA


# muAres without fgs
Ares_nofgs = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_muAres_nofgs.npz')
Ares_nofgs_freq = Ares_nofgs['x']
Ares_nofgs_strain = Ares_nofgs['y']
eff_Ares = 1.
Ares_nofgs_T_obs = 10 * year_sec * eff_Ares
Ares_R12 = np.load('/home/paolo/Codes/SGWBProbeComb/files/Responses/Resp_muAres.npy')
Ares_f_R = np.load('/home/paolo/Codes/SGWBProbeComb/files/Responses/f_R_Ares.npy')


# BBO 
BBO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_BBO_STAR.npz')
BBO_freq = BBO['x']
BBO_strain = BBO['y']
eff_BBO = 1.
BBO_T_obs = 10 * year_sec * eff_BBO


# DECIGO
DECIGO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DECIGO.npz')
DECIGO_freq = DECIGO['x']
DECIGO_strain = DECIGO['y']/3
eff_DECIGO = 1.
DECIGO_T_obs = 10 * year_sec * eff_DECIGO


# DO Optimal
DO = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DO_Optimal.npz')
DO_freq = DO['x']
DO_strain = DO['y']
eff_DO = 0.75
DO_T_obs = 4 * year_sec * eff_DO
DO_R12 = np.load('/home/paolo/Codes/SGWBProbeComb/files/Responses/Resp_DO.npy')
DO_f_R = np.load('/home/paolo/Codes/SGWBProbeComb/files/Responses/f_R_DO.npy')


# DO Conservative
DO_cons = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_DO_Conservative.npz')
DO_cons_freq = DO_cons['x']
DO_cons_strain = DO_cons['y']
eff_DO = 0.75
DO_cons_T_obs = 4 * year_sec * eff_DO

# AEDGE
AEDGE = np.load('/home/paolo/Codes/SGWBProbeComb/files/S_h_AEDGE.npz')
AEDGE_freq = AEDGE['x']
AEDGE_strain = AEDGE['y']
eff_AEDGE = 0.6
AEDGE_T_obs = 5 * year_sec * eff_AEDGE

###############################################################################
# Generate primordial signals using our class Signal_GW
class_axion_SKA = Signal_GW(r_vac=1e-5, r_star=800, k_p=1e8, sigma=5.1, axion=True, running=True)
class_axion1 = Signal_GW(r_vac=1e-5, r_star=400, k_p=1e13, sigma=8.1, axion=True, running=True)
class_axion2 = Signal_GW(r_vac=1e-5, r_star=0.15 , k_p=1e11, sigma=8, axion=True, running=True)
class_no_axion = Signal_GW(r_vac=0.01, axion=None, running=True)
class_no_axion_r0001 = Signal_GW(r_vac=0.001, axion=None, running=True)
class_no_axion_rBICEP = Signal_GW(r_vac=0.06, axion=None, running=True)

# wavenumber array
k = np.logspace(np.log10(1e-5), np.log10(1e20), 100000)

###############################################################################
# class for SKA
sens_curve_SKA = np.array(SKA_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_SKA = np.array(SKA_freq) * 6.5e14

class_binned_SKA = Binned_GW(
                         name_exp='SKA', 
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_SKA,
                         omega_gw=omega_gw,
                         k_sens=k_SKA,
                         kmin_sens=k_SKA[0],
                         N_bins_sens=5,
                         T_obs=SKA_T_obs,
                         n_det=1.,
                         interp=True,
                         sigma_L=1.0
                         )

xerr_SKA, yerr_SKA, bins_mean_point_SKA, binned_signal_SKA, binned_curve_SKA = class_binned_SKA.sens_curve_binning()
###############################################################################
# class for AEDGE

sens_curve_AEDGE = np.array(AEDGE_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_AEDGE = np.array(AEDGE_freq) * 6.5e14

class_binned_AEDGE = Binned_GW(
                         name_exp='AEDGE',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_AEDGE,
                         omega_gw=omega_gw,
                         k_sens=k_AEDGE,
                         kmin_sens=k_AEDGE[0],
                         N_bins_sens=4,
                         T_obs=AEDGE_T_obs,
                         n_det = 1.,
                         interp=True,
                         sigma_L=0.1
                         )

xerr_AEDGE, yerr_AEDGE, bins_mean_point_AEDGE, binned_signal_AEDGE, binned_curve_AEDGE = class_binned_AEDGE.sens_curve_binning()


###############################################################################
#class for ET

sens_curve_ET = np.array(ET_strain)
omega_gw = class_axion1.analytic_omega_WK(k)
k_ET = np.array(ET_freq) * 6.5e14

class_binned_ET = Binned_GW(
                         name_exp='ET',   
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_ET,
                         omega_gw=omega_gw,
                         k_sens=k_ET,
                         kmin_sens=1.5*6.5e14 ,
                         N_bins_sens=5,
                         T_obs=ET_T_obs,
                         n_det = 3.,
                         interp=True,
                         sigma_L=1.0
                         )

xerr_ET, yerr_ET, bins_mean_point_ET, binned_signal_ET, binned_curve_ET = class_binned_ET.sens_curve_binning()


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
                         delta_log_k=1.2,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=7,
                         T_obs=LISA_xcosmo_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0
                         )

binned_signal_whole, bins_mean_point_whole = class_binned.Omega_GW_binning()
xerr, yerr, bins_mean_point, binned_signal, binned_curve = class_binned.sens_curve_binning()

################################################################################
# BBO

omega_gw_axion2 = class_axion2.analytic_omega_WK(k)
k_BBO = np.array(BBO_freq) * 6.5e14 
sens_curve_BBO = np.array(BBO_strain) 

class_BBO = Binned_GW(
                        name_exp = 'BBO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=10,
                         T_obs=BBO_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0
                         )

binned_signal_axion2, bins_mean_point_axion2 = class_BBO.Omega_GW_binning()
xerr_BBO, yerr_BBO, bins_mean_point_BBO, binned_signal_BBO, binned_curve_BBO = class_BBO.sens_curve_binning()

###############################################################################
#class for LiteBIRD and r=0.01

Fisher = np.load('/home/paolo/Codes/SGWBProbeComb/files/LiteBIRD_Fisher_matrices/Fisher_0.01.npy')  
omega_gw_flat = class_no_axion.analytic_omega_WK(k)
power_spectrum = class_no_axion.tensor_spect(k)

class_binned_flat_CMB = Binned_GW(
                         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         omega_gw=omega_gw_flat,
                         kmin_sens=1e-4,
                         N_bins_sens=5,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_flat, bins_mean_point_whole_flat = class_binned_flat_CMB.Omega_GW_binning()
xerr_flat, yerr_flat, bins_mean_point_flat, binned_signal_flat, binned_curve_flat = class_binned_flat_CMB.sens_curve_binning()

###############################################################################

# class for LiteBIRD and Axion model r_vac=1e-5, r_star=835, k_p=1e13, sigma=9

Fisher_axion = np.load('/home/paolo/Codes/SGWBProbeComb/files/LiteBIRD_Fisher_matrices/Fisher_1.2_AXION1.npy')  
power_spectrum_axion = class_axion1.total_spect(k)

class_binned_axion_CMB = Binned_GW(
                         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         omega_gw=omega_gw,
                         kmin_sens=1e-4,
                         N_bins_sens=6,
                         CMB=True,
                         F=Fisher_axion,
                         tensor_spect=power_spectrum_axion,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_axion, bins_mean_point_whole_axion = class_binned_axion_CMB.Omega_GW_binning()
xerr_axion, yerr_axion, bins_mean_point_axion, binned_signal_axion, binned_curve_axion = class_binned_axion_CMB.sens_curve_binning()

################################################################################

#class for DECIGO

sens_curve_DECIGO = np.array(DECIGO_strain)  
k_decigo = class_axion1.freq_k_conv(DECIGO_freq)
class_DECIGO = Binned_GW(
                         name_exp = 'DECIGO',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw_axion2,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=9,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0
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
                         delta_log_k=1.2,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=12,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         sigma_L=1.0,
                         R_auto=Ares_R12[0],
                         R_12=Ares_R12,
                         f_R=Ares_f_R,
                         cosmic_var=True
                         )

xerr_muares_nofgs, yerr_muares_nofgs, bins_mean_point_muares_nofgs, binned_signal_muares_nofgs, binned_curve_muares_nofgs = class_MUARES_nofgs.sens_curve_binning()
                
###############################################################################

#class for DO Optimal

sens_curve_DO = np.array(DO_strain)
k_DO = class_axion1.freq_k_conv(DO_freq)
class_DO = Binned_GW(
                         name_exp = 'DO_Optimal',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DO,
                         omega_gw=omega_gw,
                         k_sens=k_DO,
                         kmin_sens=k_DO[0],
                         N_bins_sens=7,
                         T_obs=DO_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0,
                         R_auto=DO_R12[0],
                         R_12=DO_R12,
                         f_R=DO_f_R,
                         cosmic_var=None
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
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DO_cons,
                         omega_gw=omega_gw,
                         k_sens=k_DO_cons,
                         kmin_sens=k_DO_cons[0],
                         N_bins_sens=7,
                         T_obs=DO_cons_T_obs,
                         interp=True,
                         n_det = 1.,
                         sigma_L=1.0,
                         R_auto=DO_R12[0],
                         R_12=DO_R12,
                         f_R=DO_f_R,
                         cosmic_var=None
                         )

xerr_DO_cons, yerr_DO_cons, bins_mean_point_DO_cons, binned_signal_DO_cons, binned_curve_DO_cons = class_DO_cons.sens_curve_binning()

###############################################################################

#class for LiteBIRD and r=0.001 (I use it only for the r=1e-3 signal)

Fisher = np.load('/home/paolo/Codes/SGWBProbeComb/files/LiteBIRD_Fisher_matrices/Fisher_0.01.npy')  
omega_gw_flat_r0001 = class_no_axion_r0001.analytic_omega_WK(k)
power_spectrum_r0001 = class_no_axion_r0001.tensor_spect(k)
omega_gw_BICEP = class_no_axion_rBICEP.analytic_omega_WK(k)

class_binned_flat_CMB_r0001 = Binned_GW(
                        name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         omega_gw=omega_gw_flat_r0001,
                         kmin_sens=1e-4,
                         N_bins_sens=8,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum_r0001,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_flat_r0001, bins_mean_point_whole_flat_r0001 = class_binned_flat_CMB_r0001.Omega_GW_binning()

###############################################################################

#class for LiteBIRD and r=0.06 

Fisher = np.load('/home/paolo/Codes/SGWBProbeComb/files/LiteBIRD_Fisher_matrices/Fisher_0.01.npy')  
omega_gw_BICEP = class_no_axion_rBICEP.analytic_omega_WK(k)
power_spectrum_rBICEP = class_no_axion_rBICEP.tensor_spect(k)

class_binned_flat_CMB_rBICEP = Binned_GW(
                         name_exp = 'LiteBIRD',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         omega_gw=omega_gw_BICEP,
                         kmin_sens=1e-4,
                         N_bins_sens=8,
                         CMB=True,
                         F=Fisher,
                         tensor_spect=power_spectrum_rBICEP,
                         sigma_L=1.0
                         )                 
                      
binned_signal_whole_BICEP, bins_mean_point_whole_BICEP = class_binned_flat_CMB_rBICEP.Omega_GW_binning()

# FOREGROUNDS BELOW

#class for SKA for fgs

class_binned_SKA_fgs = Binned_GW(
                         name_exp='SKA_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_SKA,
                         omega_gw=omega_gw,
                         k_sens=k_SKA,
                         kmin_sens=k_SKA[0],
                         N_bins_sens=5,
                         T_obs=SKA_T_obs,
                         n_det=1.,
                         fgs=True ,
                         interp=True,
                         sigma_L=0.1
                         )

xerr_SKA_fgs, yerr_SKA_fgs, bins_mean_point_SKA_fgs, binned_signal_SKA_fgs, binned_curve_SKA_fgs = class_binned_SKA_fgs.sens_curve_binning()

###############################################################################

## Class for ET with fgs

class_binned_ET_fgs = Binned_GW(
                         name_exp='ET_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_ET,
                         omega_gw=omega_gw,
                         k_sens=k_ET,
                         kmin_sens=1.5*6.5e14 ,#k_ET[0],
                         N_bins_sens=5,
                         T_obs=ET_T_obs,
                         n_det = 3.,
                         fgs=True,
                         interp=True,
                         sigma_L=1.
                         )

xerr_ET_fgs, yerr_ET_fgs, bins_mean_point_ET_fgs, binned_signal_ET_fgs, binned_curve_ET_fgs = class_binned_ET_fgs.sens_curve_binning()


###############################################################################

#class for AEDGE with fgs

omega_gw = class_axion1.analytic_omega_WK(k)

class_binned_AEDGE_fgs = Binned_GW(
                         name_exp = 'AEDGE_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_AEDGE,
                         omega_gw=omega_gw,
                         k_sens=k_AEDGE,
                         kmin_sens=k_AEDGE[0],
                         N_bins_sens=4,
                         T_obs=AEDGE_T_obs,
                         n_det = 1.,
                         fgs=True,
                         interp=True,
                         sigma_L=0.1
                         )

xerr_AEDGE_fgs, yerr_AEDGE_fgs, bins_mean_point_AEDGE_fgs, binned_signal_AEDGE_fgs, binned_curve_AEDGE_fgs = class_binned_AEDGE_fgs.sens_curve_binning()

###############################################################################

#class for LISA with fgs

omega_gw = class_axion1.analytic_omega_WK(k)

class_binned_fgs = Binned_GW(
                         name_exp='LISA_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_LISA,
                         omega_gw=omega_gw,
                         k_sens=k_LISA,
                         kmin_sens=1.21303790e+10,
                         N_bins_sens=7,
                         T_obs=LISA_xcosmo_T_obs,
                         n_det = 1.,
                         fgs=True,
                         interp=True,
                         sigma_L=0.1
                         )

xerr_fgs, yerr_fgs, bins_mean_point_fgs, binned_signal_fgs, binned_curve_fgs = class_binned_fgs.sens_curve_binning()

################################################################################

# BBO with fgs

class_BBO_fgs = Binned_GW(
                         name_exp='BBO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_BBO,
                         omega_gw=omega_gw,
                         k_sens=k_BBO,
                         kmin_sens=k_BBO[0],
                         N_bins_sens=10,
                         T_obs=BBO_T_obs,
                         n_det = 2.,
                         fgs=True,
                         interp=True,
                         sigma_L=1.0
                         )

xerr_BBO_fgs, yerr_BBO_fgs, bins_mean_point_BBO_fgs, binned_signal_BBO_fgs, binned_curve_BBO_fgs = class_BBO_fgs.sens_curve_binning()

################################################################################
#class for DECIGO with fgs

class_DECIGO_fgs = Binned_GW(
        name_exp = 'DECIGO_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DECIGO,
                         omega_gw=omega_gw,
                         k_sens=k_decigo,
                         kmin_sens=k_decigo[0],
                         N_bins_sens=9,
                         T_obs = DECIGO_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=1.
                         )

xerr_decigo_fgs, yerr_decigo_fgs, bins_mean_point_decigo_fgs, binned_signal_decigo_fgs, binned_curve_decigo_fgs = class_DECIGO_fgs.sens_curve_binning()
                
###############################################################################

#class for DO Optimal with fgs

class_DO_fgs = Binned_GW(
        name_exp = 'DO_Optimal_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DO,
                         omega_gw=omega_gw,
                         k_sens=k_DO,
                         kmin_sens=k_DO[0],
                         N_bins_sens=7,
                         T_obs=DO_T_obs,
                         interp=True,
                         n_det = 1.,
                         fgs=True,
                         sigma_L=0.1
                         )

xerr_DO_fgs, yerr_DO_fgs, bins_mean_point_DO_fgs, binned_signal_DO_fgs, binned_curve_DO_fgs = class_DO_fgs.sens_curve_binning()

###############################################################################

#class for DO Conservative with fgs

class_DO_cons_fgs = Binned_GW(
        name_exp = 'DO_Conservative_with_fgs',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_DO_cons,
                         omega_gw=omega_gw,
                         k_sens=k_DO_cons,
                         kmin_sens=k_DO_cons[0],
                         N_bins_sens=7,
                         T_obs=DO_cons_T_obs,
                         interp=True,
                         n_det = 1.,
                         fgs=True,
                         sigma_L=1e-3
                         )

xerr_DO_cons_fgs, yerr_DO_cons_fgs, bins_mean_point_DO_cons_fgs, binned_signal_DO_cons_fgs, binned_curve_DO_cons_fgs = class_DO_cons_fgs.sens_curve_binning()

###############################################################################

#class for muAres with foregrounds

class_MUARES = Binned_GW(
                         name_exp='muAres_two_fgs_spectral',
                         kmin=1e-4,
                         k=k,
                         N_bins=80,
                         delta_log_k=1.2,
                         sens_curve=sens_curve_MUARES_nofgs,
                         omega_gw=omega_gw,
                         k_sens=k_muares_nofgs,
                         kmin_sens=k_muares_nofgs[0],
                         N_bins_sens=12,
                         T_obs=Ares_nofgs_T_obs,
                         interp=True,
                         n_det = 2.,
                         fgs=True,
                         sigma_L=0.1,
                         cosmic_var=True,
                         R_auto=Ares_R12[0],
                         R_12=Ares_R12,
                         f_R=Ares_f_R
                         )

xerr_muares, yerr_muares, bins_mean_point_muares, binned_signal_muares, binned_curve_muares = class_MUARES.sens_curve_binning()
                
###############################################################################
#PLOT

fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)
_ = make_error_boxes(ax, np.array(bins_mean_point)/6.5e14, binned_signal, xerr/6.5e14, yerr, facecolor='b', alpha=0.7, zorder=10)
_ = make_error_boxes(ax, np.array(bins_mean_point_fgs)/6.5e14, binned_signal_fgs, xerr_fgs/6.5e14, yerr_fgs, facecolor='b', alpha=0.4, zorder=9)

## ET axion
_ = make_error_boxes(ax, np.array(bins_mean_point_ET)/6.5e14, binned_signal_ET, xerr_ET/6.5e14, yerr_ET, facecolor='purple', alpha=0.7, zorder=10)
_ = make_error_boxes(ax, np.array(bins_mean_point_ET_fgs)/6.5e14, binned_signal_ET_fgs, xerr_ET_fgs/6.5e14, yerr_ET_fgs, facecolor='purple', alpha=0.4, zorder=9)

# r=0.06 nT=-r/8 signal
plt.loglog(np.array(bins_mean_point_whole_BICEP)/6.5e14, binned_signal_whole_BICEP, color='violet',label=r'Primordial Signal $r=0.06$',linewidth=1.0, zorder=1, linestyle='--')

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')


#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=1,
           alpha=0.55, linestyle='--')

#plot for LiteBIRD axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion)/6.5e14, binned_signal_axion, xerr_axion/6.5e14, yerr_axion, facecolor='g', alpha=0.55, zorder=4)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=2,
           label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$', alpha=0.55, linestyle='--')


#plot for SKA
_ = make_error_boxes(ax, np.array(bins_mean_point_SKA)/6.5e14, binned_signal_SKA, xerr_SKA/6.5e14, yerr_SKA, facecolor='orange', alpha=0.7, zorder=1)
_ = make_error_boxes(ax, np.array(bins_mean_point_SKA_fgs)/6.5e14, binned_signal_SKA_fgs, xerr_SKA_fgs/6.5e14, yerr_SKA_fgs, facecolor='orange', alpha=0.3, zorder=1)


plt.text(5e-19, 3e-11, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(2e-12, 5e-12, r'$\bf SKA$', fontsize=10, color='orange')
plt.text(3e-4, 2e-9, r'$\bf LISA$', fontsize=10, color='blue')
plt.text(3e0, 1.5e-11, r'$\bf ET$', fontsize=10, color='purple')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()
axes.set_ylim([1e-25,1e-2])

ax = 1e-20
bx = 1e3
ay = 1e-21
by = 1e-4

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_9.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
###############################################################################

#PLOT

fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')

#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=1,
           alpha=0.55, linestyle='--')

# r=0.06 nT=-r/8 signal
plt.loglog(np.array(bins_mean_point_whole_BICEP)/6.5e14, binned_signal_whole_BICEP, color='violet',label=r'Primordial Signal $r=0.06$',linewidth=1.0, zorder=1, linestyle='--')

#plot for LiteBIRD axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion)/6.5e14, binned_signal_axion, xerr_axion/6.5e14, yerr_axion, facecolor='g', alpha=0.55, zorder=4)

#plot for muAres without fgs axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_muares_nofgs)/6.5e14, binned_signal_muares_nofgs, xerr_muares_nofgs/6.5e14, yerr_muares_nofgs, facecolor=sns.xkcd_rgb["amber"], alpha=0.7, zorder=9)
_ = make_error_boxes(ax, np.array(bins_mean_point_muares)/6.5e14, binned_signal_muares, xerr_muares/6.5e14, yerr_muares, facecolor=sns.xkcd_rgb["amber"], alpha=0.3, zorder=8)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=2,
           label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$', alpha=0.55, linestyle='--')

plt.text(1e-18, 4e-12, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(5e-6, 3e-11, r'$\bf \mu Ares$', fontsize=10, color=sns.xkcd_rgb["amber"])

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()
axes.set_ylim([1e-25,1e-2])

ax = 1e-20
bx = 5e1
ay = 1e-18
by = 1e-7

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_13.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################

#PLOT

fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')


#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=1,
           alpha=0.55, linestyle='--')

#plot for LiteBIRD axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion)/6.5e14, binned_signal_axion, xerr_axion/6.5e14, yerr_axion, facecolor='g', alpha=0.55, zorder=4)

# r=0.06 nT=-r/8 signal
plt.loglog(np.array(bins_mean_point_whole_BICEP)/6.5e14, binned_signal_whole_BICEP, color='violet',label=r'Primordial Signal $r=0.06$',linewidth=1.0, zorder=1, linestyle='--')

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=2,
           label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$', alpha=0.55, linestyle='--')

#plot for DO Optimal axion spectrum 
_ = make_error_boxes(ax, np.array(bins_mean_point_DO)/6.5e14, binned_signal_DO, xerr_DO/6.5e14, yerr_DO, facecolor='black', alpha=0.7, zorder=30)
_ = make_error_boxes(ax, np.array(bins_mean_point_DO_fgs)/6.5e14, binned_signal_DO_fgs, xerr_DO_fgs/6.5e14, yerr_DO_fgs, facecolor='black', alpha=0.4, zorder=29)

plt.text(5e-19, 3e-11, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(2e-2, 1e-10, r'$\bf DO$', fontsize=10, color='black')
plt.text(2e-3, 1e-11, r'$\bf Optimal$', fontsize=10, color='black')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()
axes.set_ylim([1e-25,1e-2])

ax = 1e-20
bx = 5e1
ay = 1e-21
by = 1e-4

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_11.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################

#PLOT

fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')

#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=1,
           alpha=0.55, linestyle='--')

#plot for LiteBIRD axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion)/6.5e14, binned_signal_axion, xerr_axion/6.5e14, yerr_axion, facecolor='g', alpha=0.55, zorder=4)

#plot for BBO axion spectrum 
plt.loglog(np.array(bins_mean_point_axion2)/6.5e14, binned_signal_axion2, color='orange', linewidth=1.0, zorder=2,
           label='Axion Signal $r_{\star}=0.15$, $k_{p}=10^{11}$ $Mpc^{-1}$, $\sigma=8$', alpha=0.55, linestyle='--')

# r=0.06 nT=-r/8 signal
plt.loglog(np.array(bins_mean_point_whole_BICEP)/6.5e14, binned_signal_whole_BICEP, color='violet',label=r'Primordial Signal $r=0.06$',linewidth=1.0, zorder=1, linestyle='--')

#plot for DO Conservative axion spectrum 
_ = make_error_boxes(ax, np.array(bins_mean_point_DO_cons)/6.5e14, binned_signal_DO_cons, xerr_DO_cons/6.5e14, yerr_DO_cons, facecolor='grey', alpha=0.7, zorder=29)
_ = make_error_boxes(ax, np.array(bins_mean_point_DO_cons_fgs)/6.5e14, binned_signal_DO_cons_fgs, xerr_DO_cons_fgs/6.5e14, yerr_DO_cons_fgs, facecolor='grey', alpha=0.4, zorder=29)

plt.text(5e-19, 3e-11, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(2e-3, 1e-9, r'$\bf DO$', fontsize=10, color='grey')
plt.text(2e-5, 1e-10, r'$\bf Conservative$', fontsize=10, color='grey')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()
axes.set_ylim([1e-25,1e-2])

ax = 1e-20
bx = 5e1
ay = 1e-21
by = 1e-4

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_10.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()

###############################################################################
#PLOT

fig = plt.figure()
ax = plt.gca()

#plot for LISA axion
plt.loglog(np.array(bins_mean_point_whole)/6.5e14, binned_signal_whole, color='blue',label=r'Axion Signal $r_{\star}=400$, $k_{p}=10^{15}$ $Mpc^{-1}$, $\sigma=9.1$',linewidth=1.0, zorder=18)

# AEDGE axion
_ = make_error_boxes(ax, np.array(bins_mean_point_AEDGE)/6.5e14, binned_signal_AEDGE, xerr_AEDGE/6.5e14, yerr_AEDGE, facecolor='#17becf', alpha=0.7, zorder=10)
_ = make_error_boxes(ax, np.array(bins_mean_point_AEDGE_fgs)/6.5e14, binned_signal_AEDGE_fgs, xerr_AEDGE_fgs/6.5e14, yerr_AEDGE_fgs, facecolor='#17becf', alpha=0.4, zorder=9)

# r=0.001 signal
plt.loglog(np.array(bins_mean_point_whole_flat_r0001)/6.5e14, binned_signal_whole_flat_r0001, color='green',label=r'Primordial Signal $r=0.001$',linewidth=1.0, zorder=1, linestyle='--')

#plot for LiteBIRD flat spectrum r=0.01
plt.loglog(np.array(bins_mean_point_whole_flat)/6.5e14, binned_signal_whole_flat, color='red',label=r'Primordial Signal $r=0.01$',linewidth=1.0, zorder=1,
           alpha=0.55, linestyle='--')

#plot for LiteBIRD axion spectrum
_ = make_error_boxes(ax, np.array(bins_mean_point_axion)/6.5e14, binned_signal_axion, xerr_axion/6.5e14, yerr_axion, facecolor='g', alpha=0.55, zorder=4)

# r=0.06 nT=-r/8 signal
plt.loglog(np.array(bins_mean_point_whole_BICEP)/6.5e14, binned_signal_whole_BICEP, color='violet',label=r'Primordial Signal $r=0.06$',linewidth=1.0, zorder=1, linestyle='--')

plt.text(1e-14, 1e-13, r'$\bf LiteBIRD$', fontsize=10, color='green')
plt.text(1e-3, 5e-13, r'$\bf AEDGE$', fontsize=10, color='#17becf')

plt.xlabel(r'f $[Hz]$',fontsize = labelsize)
plt.ylabel(r'$h^{2} \Omega_{GW}$',fontsize = labelsize)
plt.tick_params(axis = 'both',which = 'major', labelsize = axissize)
plt.legend(fontsize=6, loc='upper left')#, bbox_to_anchor=(1, 0.5))
axes = plt.gca()
axes.set_ylim([1e-25,1e-2])

ax = 1e-20
bx = 5e2
ay = 5e-17
by = 1e-11

plot1 = plt.subplot(111)
plt.xscale('log')
plt.yscale('log')
plt.xlim(ax, bx)
plt.ylim(ay, by)

plt.savefig('/home/paolo/Codes/SGWBProbeComb/figures/Fig_12.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()