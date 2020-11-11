#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:36:59 2020

@author: Paolo Campeti

This module contains the class Binned_GW needed to compute the binned 
sensitivity given the sensitivity curve for a GW experiment.

"""
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from scipy import stats
import astropy.units as u
import scipy.integrate as integrate
from sgwbprobe.interpolate import log_interp1d, log_interp1d_fill

class Binned_GW:
    ''' Container of methods used to compute the binned sensitivity curve 
        given the sensitivity curve of a GW observatory.
    '''
    
    def __init__(
                 self,
                 name_exp,
                 kmin,
                 k, 
                 N_bins, 
                 delta_log_k,
                 omega_gw,
                 kmin_sens,
                 N_bins_sens, 
                 T_obs=None, 
                 tensor_spect=None,
                 k_sens=None, 
                 sens_curve=None,
                 CMB=None, 
                 F=None,
                 A_S=2.1e-9, 
                 interp=None,
                 n_det=None,
                 fgs=None,
                 sigma_L=None,
                 cosmic_var=None,
                 f_R=None,
                 R_auto=None,
                 R_12=None
                 ):
        
        '''
        Parameters
        ----------
        name_exp: string. 
            Name of the current experiment. 
        kmin: float.
            Minimum wavenumber k of the signal  
        k: numpy array.
            array of wavenumbers for the signal. 
        N_bins: integer.
            Number of bins in which the signal is binned along the whole k 
            range, only for plotting purposes. N_bins=80 is usually enough to
            cover the range from 10^-18 to 10^4 Hz, you should increase it if
            the you are using smaller bins or decrease it for larger bins.
        delta_log_k: float, the logarithm bin step.
        omega_gw: numpy array.
            The GW signal (gravitational wave energy density in GWs today in 
            the Universe) computed from the Signal_GW class.
        kmin_sens: float.
            The minimum k of the sensitivity range of the experiment.
        N_bins_sens: integer. 
            The number of bins in the experiment sensitivity range.
        T_obs: float. 
            The mission observation time in seconds.
        tensor_spect: numpy array (optional). 
            Input tensor power spectrum for the model analyzed (Eqs.(2.1) 
            and (2.5)).
        k_sens: numpy array.
            Array of wavenumbers k representing the band-width of a 
            given experiment.
        sens_curve: numpy  array. 
            The experiment instrumental strain sensitivity S_h from Eq.(4.13).         
        CMB : Boolean type (optional). 
            True if we are computing a CMB sensitivity curve, False otherwise.
        F: numpy ndarray with shape (N_bins_sens, N_bins_sens) (optional). 
            The CMB sensitivity Fisher matrix (optional, necessary only if we 
            want to compute the sensitivity for a CMB experiment).
        A_S: float (optional).
            Amplitude of the scalar perturbations spectrum, default to 
            A_S=2.1e-9.
        interp: Boolean type (optional). 
            True if you want to interpolate along the instrument bandwidth the 
            input instrumental strain sensitivity 
            curve (sens_curve).
        n_det : integer. 
            Number of detectors in the cross-correlation for interferometers 
            (see Eq.(4.10)), default to n_det=1.
        fgs: Boolean type (optional). 
            True if you want error bars including foregrounds residuals.
        sigma_L: float. 
            Fractional uncertainty on the amplitude of the BBH+BNS foreground 
            given by an external experiment (see Sec.4.2.2).
        cosmic_var: Boolean type (optional). 
            True if you want to include the cosmic variance in the 
            interferometer SNR (see Ref.[5]).
        f_R: numpy array (optional).
            Frequency band for the response function necessary only if you want
            also the cosmic variance for the interferometers. 
        R_auto: numpy array.
            The frequency response R_II for the interferometer.
        R_12: numpy array. 
            The frequency response R_IJ with I/=J for the interferometer.
        
        '''
                
        self.name_exp = name_exp
        self.kmin = kmin
        self.delta_log_k = delta_log_k
        self.k = k
        self.N_bins = N_bins
        self.sens_curve = sens_curve
        self.omega_gw = omega_gw
        self.k_sens = k_sens
        self.kmin_sens = kmin_sens
        self.N_bins_sens = N_bins_sens
        self.CMB = CMB
        self.F = F
        self.tensor_spect = tensor_spect
        self.A_S = A_S
        self.interp = interp
        self.T_obs = T_obs
        self.n_det = n_det
        self.fgs = fgs
        self.sigma_L = sigma_L
        self.f_R = f_R
        self.R_auto = R_auto
        self.R_12 = R_12
        self.cosmic_var = cosmic_var
        
    def sens_curve_binning(self): 
        '''
        Ths method bins the sensitivity curve and the GW signal and returns the
        following quantities.
        
        Returns
        -------
        xerr: numpy ndarray. 
            The half-width of the error rectangle for each bin.
        yerr: numpy ndarray.
            The half-height of the error rectangle for each bin.
        bins_mean_point: numpy array. 
            The x coordinate of the error rectangle center.
        binned_signal: numpy array.
            The GW signal binned in the experiment sensitivity range.
        binned_curve: numpy array.
            The binned sensitivity curve of the experiment.  
        
        '''
        print('Forecast for ' + self.name_exp)
        k_i = self.initialize_k_i_sens()
        #append kmax value to the wavenumbers array
        bins_array = np.append(
                               k_i, 
                               self.kmin_sens 
                               * np.exp(self.delta_log_k)**len(k_i)
                               ) 
        #binning the signal
        bin_means_GW, _, _ = stats.binned_statistic(
                                                    self.k, 
                                                    self.omega_gw, 
                                                    statistic='mean',
                                                    bins=bins_array
                                                    )
        #binned Omega_GW signal in the sensitivity range
        binned_signal = bin_means_GW         
        bins_mean_point = [
                           np.exp((np.log(x) + np.log(bins_array[i - 1]))/2)
                           for i, x in enumerate(bins_array)
                           ][1:]
        
        signal = np.insert(binned_signal, 0, binned_signal[0], axis=0)
        signal = np.append(signal, binned_signal[-1])
        k_signal = np.insert(bins_mean_point, 0, bins_array[0], axis=0)
        k_signal = np.append(k_signal, bins_array[-1])
               
        
        if self.CMB: 
            ''' 
            Compute the binned sensitivity for a CMB experiment, as in 
            Eq.(3.12).
            '''
            sigma = self.CMB_Sensitivity_Calculator()
            transfer = self.A_S * self.omega_gw / self.tensor_spect 
            # binning the tensor spectrum
            bin_means_transfer, _, _ = stats.binned_statistic(
                                                           self.k, 
                                                           transfer, 
                                                           statistic='mean',
                                                           bins=bins_array
                                                           )
            
            binned_curve =  sigma[:self.N_bins_sens] * bin_means_transfer 
                    
        
        else:
            ''' 
            Compute the binned sensitivity for interferometers or PTA.
            '''
            if self.interp:
                '''
                This interpolates the sensitivity curve if needed.
                '''
                f = log_interp1d(self.k_sens, self.sens_curve)
                k_interp = np.logspace(np.log10(self.k_sens[0]), 
                                       np.log10(self.k_sens[-1]), 
                                       100000)
                interp_sens = f(k_interp)
                # binned sensitivity curve
                binned_curve = self.time_int(k_signal,
                                             signal,
                                             k_interp, 
                                             bins_array, 
                                             interp_sens, 
                                             self.n_det
                                             )  

            else:
                '''
                If there is no need to interpolate the sensitivity curve
                '''
                # binned sensitivity curve
                binned_curve = self.time_int(k_signal,
                                             signal,
                                             self.k_sens,
                                             bins_array, 
                                             self.sens_curve,
                                             self.n_det
                                             )  
        
                
        # bins half-width in a list        
        bins_diff = [
                    (self.kmin_sens*np.exp(self.delta_log_k/2)**(i) 
                     - self.kmin_sens*np.exp(self.delta_log_k/2)**(i-1)) 
                     for i in range(1, 2*len(k_i)+1)
                     ] 
        # x error (half-width of bin on left and right of center point) array 
        # for the make_error_boxes function
        xerr = np.vstack((bins_diff[::2], bins_diff[1::2])) 
        # y error (above and below center point) array for 
        # the make_error_boxes function
        yerr = np.vstack((binned_curve, binned_curve)) 
         
        return xerr, yerr, bins_mean_point, binned_signal, binned_curve


    def time_int(self, k_signal, signal, k, bins_array, sens_curve, n_det=1.):
        '''
        Computes the binned sensitivity to Omega_GW (including a factor h^2).
        
        Returns
        -------
        res: numpy array.
            The binned sensitivity to Omega_GW.
        ''' 
        sigma_L = self.sigma_L
        SNR0=1.
        T = self.T_obs
        n = n_det            
        f  = k / 6.5e14
        f_signal  = k_signal / 6.5e14
        fbins_array = bins_array / 6.5e14
        h = 0.6736
        H_0 = 100 * h * u.kilometer/u.second/u.megaparsec
        H_0 = H_0.to_value(u.hertz) 
        fact = (4*np.pi**2) / (3*H_0**2)
        S_s = (1./fact) * f_signal**(-3) * signal /h**2 
        S_s_interp = log_interp1d(f_signal, S_s)
        S_n_interp = log_interp1d(f, sens_curve)
            

        def M(x):
            if self.cosmic_var:
                '''
                If True includes cosmic variance for the interferometers. 
                '''
                R_auto = self.R_auto
                R_12_interp = log_interp1d_fill(self.f_R, self.R_12, a=self.R_12[0], b=1e-4)
            
                if self.name_exp=='LISA' or self.name_exp=='LISA_with_fgs':
                    '''
                    Special case for TDI variables (valid for LISA and DO).
                    '''
                    S_n = S_n_interp(x) * R_12_interp(x) * np.sqrt(2)
                else:    
                    S_n = S_n_interp(x) * R_12_interp(x)
                
                res = S_n_interp(x)**2 * (
                                          1 + 2*S_s_interp(x)*R_auto/S_n + 
                                          2*S_s_interp(x)**2 *
                                          (R_auto**2 + R_12_interp(x)**2)/S_n**2 
                                          )
            else:
                res = S_n_interp(x)**2
            return res
        
        def S_BBH(f):
            '''
            BBH+BNS foreground from Eq.(4.18) in the paper.
            '''
            Omega_star = 8.9e-10
            f_star = 25
            Omega = Omega_star * (f/f_star)**(2/3)
            res = Omega * f**(-3) / fact
            return res
        
        def S_MBHB(f):
            '''
            MBHB foreground from Eq.(4.19) in the paper.
            '''
            h0 = 0.69e-15 
            f0 = 4.27e-8
            gamma = -1.08
            hc = h0 * (f/f0)**(-2/3) * (1 + f/f0)**gamma
            res = hc**2 / f
            return res

        
        if self.fgs:
            '''
            Computes binned error bars for interferometers and PTA including
            foreground residuals.
            '''
            if (self.name_exp == 'ET_with_fgs' or
                self.name_exp == 'DECIGO_spectral_shape' or 
                self.name_exp == 'muAres_spectral_shape'):
                '''
                Computes error bars in the case sigma_fg = infinity in the 
                filter in Eq.(4.27) for ET, DECIGO and muAres, that is 
                subtraction using only the spectral shape and not external 
                experiments.
                '''
                
                S_fg = S_BBH
                
                def integrand_A(x):
                    res = 1./(x**6 * M(x))
                    return res
    
                def integrand_B(x):
                    res = S_fg(x)/(x**3 * M(x))
                    return res
                
                def integrand_C(x):
                    res = S_fg(x)**2/M(x)
                    return res
  
                res = np.zeros((len(f)))    
                for index in range(len(bins_array[:-1])):
        
                    A = integrate.quad(
                                        integrand_A, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    B = integrate.quad(
                                        integrand_B, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
        
                    C = integrate.quad(
                                        integrand_C, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    W = B/(n*T*C)
                
                    def integrand_D(x):
                        res = 1./M(x) * (np.abs(1/x**3 - n*T*W*S_fg(x)))**2 
                        return res
                
                    def integrand_E(x):
                        res = W * (S_fg(x)**2)/M(x)
                        return res

                    def integrand_F(x):
                        res = W * S_fg(x)/(x**3 * M(x))
                        return res

                    D = integrate.quad(
                                       integrand_D, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    E = integrate.quad(
                                       integrand_E, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    F = integrate.quad(
                                       integrand_F, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                                
                    num = (T/n)*D 
                    den = T**2 * (A -n*T*F)**2
                    res[index] = h**2 * SNR0 * fact * np.sqrt(num/den)

            elif self.name_exp == 'SKA_with_fgs':
                '''
                Case sigma=infinity in the filter Eq.(4.27) for SKA, that is
                subtraction using only the spectral shape of the MBHB 
                foreground, without any external information.
                '''
                
                S_fg = S_MBHB  
                
                def integrand_A(x):
                    res = 1./(x**6 * M(x))
                    return res
    
                def integrand_B(x):
                    res = S_fg(x)/(x**3 * M(x))
                    return res
                
                def integrand_C(x):
                    res = S_fg(x)**2/M(x)
                    return res
  
                res = np.zeros((len(f)))    
                for index in range(len(bins_array[:-1])):
        
                    A = integrate.quad(
                                        integrand_A, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    B = integrate.quad(
                                        integrand_B, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
        
                    C = integrate.quad(
                                        integrand_C, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    W = B/(n*T*C)
                
                    def integrand_D(x):
                        res = 1./M(x) * (np.abs(1/x**3 - n*T*W*S_fg(x)))**2 
                        return res
                
                    def integrand_E(x):
                        res = W * (S_fg(x)**2)/M(x)
                        return res

                    def integrand_F(x):
                        res = W * S_fg(x)/(x**3 * M(x))
                        return res

                    D = integrate.quad(
                                       integrand_D, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    E = integrate.quad(
                                       integrand_E, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    F = integrate.quad(
                                       integrand_F, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                                
                    num = (T/n)*D 
                    den = T**2 * (A -n*T*F)**2
                    res[index] = h**2 * SNR0 * fact * np.sqrt(num/den)
                
            
            elif self.name_exp == 'muAres_two_fgs_spectral':           
                '''
                Implements the two-fgs filter for muAres MBHB and BBH+BNS 
                foregrounds in the case the MBHB fg is subtracted only 
                exploiting its spectral dependence (that is sigma2=infinity).
                '''
                S_fg1 = S_BBH
                S_fg2 = S_MBHB
                sigma1 = 1e-3 # external sigma_fg for the BBH+BNS foregrounds
                
                def integrand_F1_F2(x):
                    res = S_fg1(x)*S_fg2(x)/M(x)
                    return res
    
                def integrand_F1_F1(x):
                    res = S_fg1(x)**2/M(x)
                    return res
                
                def integrand_F2_F2(x):
                    res = S_fg2(x)**2/M(x)
                    return res
      
                def integrand_fg1(x):
                    res = S_fg1(x)/(x**3 * M(x))
                    return res
                
                def integrand_fg2(x):
                    res = S_fg2(x)/(x**3 * M(x))
                    return res
  
                res = np.zeros((len(f)))    
                for index in range(len(bins_array[:-1])):
        
                    F1_F2 = integrate.quad(
                                        integrand_F1_F2, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    F1_F1 = integrate.quad(
                                        integrand_F1_F1, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
        
                    F2_F2 = integrate.quad(
                                        integrand_F2_F2, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                
                
                    integral_fg1 = integrate.quad(
                                        integrand_fg1, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    integral_fg2 = integrate.quad(
                                        integrand_fg2, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
        
                    
                    pref_1 = T**(-2)*sigma1**(-2) + 2/T * F1_F1
                    pref_2 =  2/T * F2_F2
                    W = (2/T) * 1/(pref_1 * pref_2 - 4/(T**2) * F1_F2**2)
                    
                
                    def integrand_abs(x):
                        res = ((1./M(x)) * (np.abs(1/x**3 - W*(S_fg1(x)*(pref_2*integral_fg1 - 2/T * F1_F2 * integral_fg2) 
                        + S_fg2(x)*(pref_1*integral_fg2 - 2/T * F1_F2 * integral_fg1))))**2) 
                        return res
                
                    def integrand_sigma1(x):
                        res = (S_fg1(x) / M(x)) * (1/x**3 - W*(S_fg1(x)*(pref_2*integral_fg1 - 2/T * F1_F2 * integral_fg2) 
                        + S_fg2(x)*(pref_1*integral_fg2 - 2/T * F1_F2 * integral_fg1)))
                        return res
                    
                    def integrand_mu(x):
                        res =  1/(x**3 * M(x)) * (1/x**3 - W*(S_fg1(x)*(pref_2*integral_fg1 - 2/T * F1_F2 * integral_fg2) 
                        + S_fg2(x)*(pref_1*integral_fg2 - 2/T * F1_F2 * integral_fg1)))
                        return res


                    integral_abs = integrate.quad(
                                       integrand_abs, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    integral_sigma1 = integrate.quad(
                                       integrand_sigma1, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                                    
                    integral_mu = integrate.quad(
                                       integrand_mu, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                
                    num = (T/n)*integral_abs + sigma1**2 * T**2 * (integral_sigma1)**2 
                    den = T**2 * (integral_mu)**2
    
                    res[index] = h**2 * SNR0 * fact * np.sqrt(num/den)

            
            else:
                '''
                Computes binned error bars using the filter in Eq.(4.27) 
                including foregrounds for all other experiments. 
                '''
                
                S_fg = S_BBH
                
                def integrand_A(x):
                    res = 1./(x**6 * M(x))
                    return res
    
                def integrand_B(x):
                    res = S_fg(x)/(x**3 * M(x))
                    return res
                
                def integrand_C(x):
                    res = S_fg(x)**2/(M(x))
                    return res
  
                res = np.zeros((len(f)))    
                for index in range(len(bins_array[:-1])):
        
                    A = integrate.quad(
                                        integrand_A, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    B = integrate.quad(
                                        integrand_B, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
        
                    C = integrate.quad(
                                        integrand_C, 
                                        fbins_array[index], 
                                        fbins_array[index+1]
                                        )[0]
                    
                    W = B/(sigma_L**(-2) + n*T*C)
                
                    def integrand_D(x):
                        res = (1./(M(x))) * (np.abs(1/x**3 - n*T*W*S_fg(x)))**2 
                        return res
                
                    def integrand_E(x):
                        res = W * (S_fg(x)**2)/(M(x))
                        return res

                    def integrand_F(x):
                        res = W * S_fg(x)/(x**3 * M(x))
                        return res

                    D = integrate.quad(
                                       integrand_D, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    E = integrate.quad(
                                       integrand_E, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                    F = integrate.quad(
                                       integrand_F, 
                                       fbins_array[index], 
                                       fbins_array[index+1]
                                       )[0]
                
                
                    num = (T/n)*D + sigma_L**2 * T**2 * (B-n*T*E)**2
                    den = T**2 * (A -n*T*F)**2
                                
                    res[index] = h**2 * SNR0 * fact * np.sqrt(num/den)
                
        else:
            '''
            Binned error bars NOT including foregrounds.
            '''
            def integrand(x):
                return 1./(x**6 * M(x)) 
        
            res = np.zeros((len(f)))    
            for index in range(len(bins_array[:-1])):
                res[index] = h**2 * SNR0 * (n * T * (1./fact)**2  
                                            * integrate.quad(
                                                             integrand, 
                                                             fbins_array[index], 
                                                             fbins_array[index+1]
                                                             )[0]
                                            )**(-1/2)
                                            
        return res

    def Omega_GW_binning(self):
        '''
        This method bins the GW signal in wavenumber k.
        '''
        k_i = self.initialize_k_i()
        #append kmax value to the wavenumbers array
        bins_array = np.append(
                               k_i, 
                               self.kmin*np.exp(self.delta_log_k)**len(k_i)
                               ) 
        bin_means, _, _ = stats.binned_statistic(
                                                 self.k, 
                                                 self.omega_gw, 
                                                 statistic='mean',
                                                 bins=bins_array
                                                 )
        binned_signal = bin_means #binned Omega_GW
        #x coordinate of center bin point for make_error_boxes function
        bins_mean_point = [np.exp((np.log(x) + np.log(bins_array[i - 1]))/2) 
                           for i, x in enumerate(bins_array)
                           ][1:] 
        return binned_signal, bins_mean_point


    def initialize_k_i(self):
        k_i = np.zeros((self.N_bins))
        for i in range(self.N_bins):
            k_i[i] = self.kmin*np.exp(self.delta_log_k)**i
        return k_i

    
    def initialize_k_i_sens(self):
        k_i = np.zeros((self.N_bins_sens))
        for i in range(self.N_bins_sens):
            k_i[i] = self.kmin_sens*np.exp(self.delta_log_k)**i
        return k_i
       
    def CMB_Sensitivity_Calculator(self):
        '''
        Computes binned errors for CMB experiments from the input Fisher 
        matrix F.
        
        Returns
        -------
        sigma: numpy array.
            The binned error from Fisher matrix, as in Eq.(3.11).
        '''
        F = self.F
        C = np.linalg.inv(F)
        sigma_square = C.diagonal() #sigma^2
        sigma = np.sqrt(sigma_square) #sigma
        #collapse (1,int(nn)+1) shape array to (int(nn)+1) shape array
        sigma = np.ravel(sigma) 
        return sigma
