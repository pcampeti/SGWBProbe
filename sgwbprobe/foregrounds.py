#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:42:16 2020

@author: Paolo Campeti

This module contains the functions used to compute the analytical 
approximations for the GWD, EGWD, MBHB and BBH+BNS foreground component for the
interferometers and PTA (see Section 4.2.1).
    
"""
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import astropy.units as u


def GWD_fg(f, T_obs):
    """
    Galactic WD foreground (see Eq.(4.15) and Ref.[133]).
    """
    # seconds in a year
    year_sec = 60*60*24*365 

    # subtraction parameters
    A = 9e-45
    a = np.array([0.133,0.171,0.165,0.138])
    b = np.array([243,292,299,-221])
    k = np.array([482,1020,611,521])
    g = np.array([917,1680,1340,1680])
    f_k = np.array([0.00258,0.00215,0.00173,0.00113])

    # mission duration time
    if T_obs < 1.*year_sec:
        index = 0
    elif T_obs >= 1.*year_sec and T_obs < 2.*year_sec:
        index = 1
    elif T_obs >= 2.*year_sec and T_obs < 4.*year_sec:
        index = 2
    else:
        index = 3
     
    # strain 
    S_c_f = (
            A * np.exp(-(f**a[index]) + (b[index] * f * np.sin(k[index] * f)))
            * (f**(-7/3)) * (1 + np.tanh(g[index] * (f_k[index] - f)))
            ) 
    return S_c_f

def EGWD_fg(f):
    """
    Extragalactic WD foreground (see Eq.(4.16) and Ref.[131]).
    """
    A = 4.2e-47
    res = np.zeros((len(f)))
    for i,freq in enumerate(f):    
        if freq >=3e-3:
        # strain     
            res[i] = A * freq**(-7/3) * np.exp(-2*(freq/5e-2)**2) 
        else:
            res[i] = np.NaN
    return np.array(res)


def MBHB_fg(f):
    '''
    MBHB foreground (see Eq.(4.19) and Ref.[134]).
    '''
    h0 = 0.69e-15 
    f0 = 4.27e-8
    gamma = -1.08
    hc = h0 * (f/f0)**(-2/3) * (1 + f/f0)**gamma
    return hc


def BBH_BNS_fg(f):
    '''
    BBH+BNS foreground (see Eq.(4.18) and Ref.[90]).
    '''
    Omega_star = 8.9e-10
    f_star = 25
    Omega = Omega_star * (f/f_star)**(2/3)
    h = 0.67
    H_0 = 100 * h * u.kilometer/u.second/u.megaparsec
    H_0 = H_0.to_value(u.hertz) 
    fact = (3*H_0**2)/(4*np.pi**2)        
    res = fact * Omega * f**(-3)
    return res