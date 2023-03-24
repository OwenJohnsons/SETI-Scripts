"""
Code Purpose: This module contains functions for EIRP calculations
Last Major Update: 16/01/2023
Author: Owen A. Johnson
Note: The code that this was based on was written by Dr. Vishal Gajjar. A lot of the snippets of code were taken from his code and modified to fit my needs.
"""

import numpy as np 

def obs_len(SNR, SEFD, distance, EIRP): 
    dist_conv = 3.24078e-20 # 1 parsec = 3.24078e-20 light years
    distance = dist_conv*distance # convert distance to light years
    lytom = 9.46073e15 # 1 light year = 9.46073e15 meters
    obslen = 0.5*pow((SNR*SEFD*4*np.pi*pow((distance*lytom),2)/EIRP),2)
    return obslen
    
def obsEIRP(SNR, SEFD, distance, Obslen, chbw, dc):
    '''
    SNR = Signal to Noise Ratio
    SEFD = System Equivalent Flux Density
    distance = distance to target in parsecs
    chbw = channel bandwidth in MHz
    dc = duty cycle
    '''
    dist_conv = 3.26156 # 1 parsec = 3.26156 light years
    distance = dist_conv*distance # convert distance to light years
    lytom = 9.46073e15 # 1 light year = 9.46073e15 meters
    EIRP = SNR*(SEFD/chbw)*4*np.pi*pow((distance*lytom),2)*np.sqrt(chbw/(2*Obslen*dc))
    return EIRP

def obsEIRPtrans(SNR, SEFD, distance):
    '''
    EIRP for transient sources (broadband)
    SNR = Signal to Noise Ratio
    SEFD = System Equivalent Flux Density
    distance = distance to target in parsecs
    '''
    dist_conv = 3.26156 # 1 parsec = 3.26156 light years
    distance = dist_conv*distance # convert distance to light years
    lytom = 9.46073e15 # 1 light year = 9.46073e15 meters
    EIRP = SNR*SEFD*4*np.pi*pow((distance*lytom),2)/np.sqrt(2.0*1.0e-3*80e6)
    return EIRP

def EIRP_appending(empty_array, SNR, SEFD, distance, Obslen, chbw, dc):
    empty_array.append(np.log10(obsEIRP(SNR, SEFD, distance, Obslen, chbw, dc)))
    return empty_array

def calc_SEFD(A, Tsys, eff=1.0):
    """ Calculate SEFD
    Tsys = system temperature
    A = collecting area
    Ae = effective collecting area
    eff = aperture efficency (0.0 to 1.0)
    """
    kb = 1.3806488e3  # 1.38064852e-23    Boltzmann constant
    Ae = A*eff
    return 2 * Tsys * kb / Ae

def calc_DishArea(d):
    """ Compute dish area
    d = dish diameter
    """
    return np.pi * (d/2)**2
    