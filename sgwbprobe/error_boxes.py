#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:39:30 2020

@author: Paolo Campeti

This module contains the method make_error_boxes used to plot the error 
rectangles.

"""
import warnings
warnings.filterwarnings("ignore")
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def make_error_boxes(ax, xdata, ydata, xerror, yerror, facecolor='r',
                     edgecolor='black', alpha=0.5, zorder=10):
    '''
    Function used to create the error bars on Figures 9-21 in the paper.
    '''
    # Create list for all the error patches
    errorboxes = []
    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(
                         errorboxes, 
                         facecolor=facecolor,
                         alpha=alpha,
                         edgecolor=edgecolor,
                         zorder=zorder
                         )
    # Add collection to axes
    ax.add_collection(pc)
    return pc
