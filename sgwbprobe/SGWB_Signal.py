#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:23:21 2020

@author: Paolo Campeti

This module contains the class Signal_GW needed to compute the energy density 
for the SU(2) Axion model of inflation and the standard signle-field slow-roll 
one.   

"""

import numpy as np
import scipy
import warnings
warnings.filterwarnings("ignore")
import scipy.interpolate
from sgwbprobe.effective_degrees import g_of_k


class Signal_GW:
    """
    Container of methods useful to compute the energy density of gravitational 
    waves for the single-field slow-roll and axion-SU(2) models described in 
    Section 2 of the paper.
    
    Parameters
    ----------
    r_vac: float.
        Tensor-to-scalar ratio for quantum vacuum fluctuations (named simply r
        in the paper).
    n_T: float (optional).
        Tensor spectral index. If None is calculated from the
        inflationary consistency relation.
    r_star: float (optional).
        Parameter of the axion-SU(2) spectrum (see Sec.2.2).        
    k_p: float (optional).
        Parameter of the axion-SU(2) spectrum (see Sec.2.2).
    sigma: float (optional).
        Parameter of the axion-SU(2) spectrum (see Sec.2.2).
    axion: Boolean type (optional), defaults to None.
        If True computes Omega_GW for the axion-SU(2) model, otherwise for the
        standard single-field slow-roll model.
    k_star: float (optional).
        Pivot scale of the tensor power spectrum. Default value 0.05.
    running: Boolean type (optional).
        If True includes the running in the tensor power spectrum according to
        the inflatonary consistency relation (see Sec.2.1).
        
    """
    def __init__(self, r_vac, nT=None, r_star=None, k_p=None, sigma=None, 
                 axion=None, k_star=0.05, running=None):

        self.A_S = 2.1e-9
        self.ns = 0.9649 
        self.k_star = k_star
        self.r_vac = r_vac 
        self.A_T = self.r_vac * self.A_S
        self.r_star = r_star
        self.k_p = k_p
        self.sigma = sigma
        self.tau_0 = 1.41745198407190182479E+04#1.41745e4
        self.tau_eq = 420 #Mpc^-1
        self.h = 0.6736
        self.C = 2.99792458e8 #speed of light
        self.H0 = (self.h *1e5)/self.C #H0 in Mpc^-1 
        self.k_eq = 1/self.tau_eq   
        self.axion = axion
        self.running = running
        
        if nT is not None:
            self.nT = nT
        else:    
            self.nT = -self.r_vac/8
        
    def total_spect(self, k):
        '''
        Returns the total tensor spectrum: vacuum + sources. 
        '''
        P_h_total = self.sourced_spect(k) + self.tensor_spect(k)
        return P_h_total

    
    def sourced_spect(self, k):
        '''
        Returns the sourced spectrum of Eq.(2.3)-(2.4) in the paper.
        '''
        P_s = self.scalar_spect(k)
        r_star = self.r_star
        index = -(1./(2.*self.sigma**2)) * (np.log(k/self.k_p))**2
        res = r_star * P_s * np.exp(index)
        return res
 
     
    def scalar_spect(self, k):
        '''
        Returns the scalar vacuum spectrum.
        '''
        P_xi = self.A_S * (k/self.k_star)**(self.ns-1)
        return P_xi
    
    def tensor_spect(self, k):
        '''
        Returns the tensor vacuum spectrum.
        '''
        if self.running is not None:
            P_h = self.A_T * (k/self.k_star)**(-self.r_vac/8 + self.r_vac/16 * ((self.ns-1)+ self.r_vac/8)*np.log(k/self.k_star))
        else:
            P_h = self.A_T * (k/self.k_star)**(self.nT)
        return P_h

        
    def freq_k_conv(self, f):
        '''
        Returns the wavenumber in Mpc^-1 corresponding to the freq f in Hz.
        '''
        k = 6.5e14*f
        return k
        
    def analytic_omega_WK(self, k):
        '''
        Returns the piecewise analytical approximation for the Omega_GW signal 
        today of Eq.(2.9) in the paper (as in Watanabe & Komatsu 2006). 
        Includes also the effect of change in the effective degrees of
        freedom in the early Universe as described in Section 2.3.
        
        Returns
        -------
        res: numpy array.
            Array containing the Omega_GW signal today as a function of
            wavenumber k.
        '''
        def A(self, k):
            '''
            Eq.(B9) of Watanabe & Komatsu 2006.
            '''
            x = k*self.tau_eq            
            A = 3/(2*x) - np.cos(2*x)/(2*x) + np.sin(2*x)/(x**2)
            return A
        
        def B(self, k):
            '''
            Eq.(B10) of Watanabe & Komatsu 2006.
            '''
            x = k*self.tau_eq
            B = -1 + 1/(x**2) - np.cos(2*x)/(x**2) - np.sin(2*x)/(2*x)
            return B
        
        def Omega_intermediate(self, k):
            '''
            Eq.(20) of Watanabe & Komatsu 2006.
            '''        
            x = k*self.tau_0
            teq_over_t02 = (self.tau_eq / self.tau_0)**2
            j2 = scipy.special.spherical_jn(2, x)
            y2 = -(1/x) * ((3/(x**2) - 1) * np.cos(x) + (3/x)*np.sin(x))
            Ak = A(self, k)
            Bk = B(self, k)
            h2 = self.h**2
            k2 = k**2
            den = 1/(12*self.H0**2)
            if self.axion:
                spect = self.total_spect(k)
            else:
                spect = self.tensor_spect(k)
            res = h2 * spect * k2 * teq_over_t02 * den * (Ak*j2 + Bk*y2)**2
            return res
   
        def Omega_super(self, k):
            '''
            Eq.(21) of Watanabe & Komatsu 2006. 
            '''
            x = k*self.tau_0
            j2 = scipy.special.spherical_jn(2, x)
            h2 = self.h**2
            k2 = k**2
            den = 1/(12*self.H0**2)
            if self.axion:
                spect = self.total_spect(k)
            else:
                spect = self.tensor_spect(k)
            res = h2 * spect * k2 * den * (3*j2/x)**2
            return res
        
        
        Omega = np.zeros((len(k)))
        for counter, k_val in enumerate(k): 
            if k_val < self.k_eq:
                Omega[counter] = Omega_super(self, k_val)
                                
            else:
                Omega[counter] = Omega_intermediate(self, k_val) 

        # Effect of change in the effective degrees of freedom in the early 
        # Universe as described in Section 2.3.
        grho_k, gs_k, _, _ = g_of_k(k)
        suppression = (grho_k/3.383) * (gs_k/3.931)**(-4/3)
        
        res = Omega * suppression
        
        return res  
