#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:41:59 2020

@author: Paolo Campeti

This module contains the function g_of_k used to compute the number of 
effective degrees of freedom $g_{*rho}$ and $g_{*s}$ as a function of 
wavenumber k. 
Uses Eq.(2.15) of Saikawa & Shirai 2018 (Ref.[5]) to relate 
the frequency of GWs to the temperature T_hc at which the corresponding 
mode re-enters the horizon. Uses also functions tabulated and provided within 
Ref.[5].
"""
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from sgwbprobecomb.interpolate import log_interp1d


def g_of_k(k):
    '''
    Effective degrees of freedom $g_{*rho}$ and $g_{*s}$ as a function of 
    wavenumber k. Uses Eq.(2.15) of Saikawa & Shirai 2018 (Ref.[5]) to relate 
    the frequency of GWs to the temperature T_hc at which the corresponding 
    mode re-enters the horizon. Uses functions tabulated and provided within 
    Ref.[5].
    '''
    
    freq = k/6.5e14
    
    # Load tabulated functions provided by Ref.[5].
    filetoread = '/home/paolo/Codes/SGWBProbeComb/files/effective_degrees.dat' 
    dataset = np.genfromtxt(filetoread, names=None, skip_header=0) #mind the skip_header
    T = dataset[:,0] #[GeV]
    grho = dataset[:,1] 
    gs =  dataset[:,3]
    
    func_grho = log_interp1d(T, grho) # interpolating function
    func_gs = log_interp1d(T, gs) # interpolating function
    T_interp = np.logspace(np.log10(T[0]), np.log10(T[-1]), 100000) #interpolated temperature
    grho_interp = func_grho(T_interp) # interpolated function
    gs_interp = func_gs(T_interp) # interpolated function
    
    # compute array of f0 today corresponding to a certain array of T_hc
    f = 2.65 * (grho_interp/105.25388)**(1/2) * (gs_interp/105.25245)**(-1/3) * (T_interp/1e8)# [Hz], T_hc in GeV
    
    func_grho_f = log_interp1d(f, grho_interp) # interpolating function g_{*}(f)  instead of g_{*}(T_in) at horizon crossing
    func_gs_f = log_interp1d(f, gs_interp) # interpolating function g_{*}(f)  instead of g_{*}(T_in) at horizon crossing
    
    # piecewise function to accomodate also value outside the range of tabulated effectuve degrees of freedom 
    grho_k = np.piecewise(freq, [freq < f[0],  (f[0] <= freq) * (freq <= f[-1]), freq > f[-1] ], [3.383, func_grho_f, 105.25388])
    gs_k = np.piecewise(freq, [freq < f[0],  (f[0] <= freq) * (freq <= f[-1]), freq > f[-1] ], [3.931, func_gs_f, 105.25245])
    
    return grho_k, gs_k, func_grho, func_gs
