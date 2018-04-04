#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modeling the Universe - Problemset 7 
Problem 1
@author: Giuliana Noto (gn2244)
"""
import math
import numpy as np
import kplr
import my_transit #Professor Bryan's script such as for PSet3 

# Find the target KOI.
client = kplr.API()
koi = client.koi(97.01)

# Get a list of light curve datasets.
lcs = koi.get_light_curves(short_cadence=False)

# Open the first dataset and read it
f = lcs[0].open()
hdu_data = f[1].data
time = hdu_data["time"]  # get the time of each observation
flux = hdu_data["sap_flux"] # get the flux
flux_err = hdu_data["sap_flux_err"] # get the error in the flux
f.close()

#Defining functions to integrate, and intensity 
def I(r):
    '''A Limb-darkening function'''
    mu = (1 - (r**2))**(0.5)
    return(1 - (1 - (mu**(0.5))))

def func1(r, p, z):
    return(I(r) * (1 - my_transit.delta(p,r,abs(z))) * 2 * r)

def func2(r, p, z):
    return(I(r) * 2 * r)

#Separate flux,time,flux error for 124<t<125
extracted_data = np.where((125>time) & (time>124))
time     = time[extracted_data]
flux     = flux[extracted_data]
flux_err = flux_err[extracted_data]

for i in range(5):
    #Compute flux mean,standard dev 
    flux_mean = np.mean(flux)
    flux_std  = np.std(flux)
    #Exclude transit - remove points 2sigma away from mean
    boolean_values = (abs(flux-flux_mean)/flux_std<2)
    flux = flux[boolean_values]

normalized_flux = flux/flux_mean 




