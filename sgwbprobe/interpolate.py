#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:40:58 2020

@author: Paolo Campeti

This module contains the functions log_interp1d and log_interp1d_fill which can 
be used to interpolate in logarithimc space.
"""
import numpy as np
import scipy
import warnings
warnings.filterwarnings("ignore")
import scipy.interpolate


def log_interp1d(xx, yy, kind='linear'):
    '''
    Interpolator used in the code above.
    '''
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

def log_interp1d_fill(xx, yy, a, b, kind='linear'):
    '''
    Interpolator used in the code above.
    '''
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, fill_value=(np.log10(a),np.log10(b)), kind=kind, bounds_error=False)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp
