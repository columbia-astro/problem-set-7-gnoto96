#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modeling the Universe - Problemset 7 
@author: Giuliana Noto (gn2244)
"""
import numpy as np
import matplotlib.pyplot as plt 
import my_transit #Professor Bryan's script such as for PSet3 
import scipy
plt.style.use('ggplot')

'''*********************
Problem 1
Model transit to data (normalized)
**********************'''

time, flux, flux_err = np.loadtxt("KOI97.01_1.out", unpack=True)

#Separate flux,time,flux error for 124<t<125
extracted_data = np.where((125>time) & (time>124))
time     = time[extracted_data]
flux     = flux[extracted_data]
flux_err = flux_err[extracted_data]

#Defining functions to integrate, and intensity 
def I(r):
    '''A Limb-darkening function'''
    mu = (1 - (r**2))**(0.5)
    return(1 - (1 - (mu**(0.5))))

def func1(r, p ,z):
    return(I(r) * (1 - my_transit.delta(p,r,abs(z))) * 2 * r)

def func2(r, p, z):
    return(I(r) * 2 * r)

'''Exclude Transit'''
#keep original flux & flux_err values
orig_flux = flux
orig_fluxerr = flux_err

for i in range(5):
    #Compute flux mean,standard dev 
    flux_mean = np.mean(flux)
    flux_std  = np.std(flux)
    #to exclude transit - remove points 2sigma away from mean
    boolean_values = (abs(flux-flux_mean)/flux_std<2)
    flux = flux[boolean_values]
    flux_err = flux_err[boolean_values]

norm_flux = orig_flux/flux_mean #normalize flux by mean 
norm_flux_err  = orig_fluxerr/flux_mean #normalize std by mean 

'''Calculate obscured to unobscured flux ratio (Problem 3 reference)'''    
def model_flux(t,tau,t0,p):
    flux_ratio = []
    for i in t:
        z = (i - t0)/tau
        func1_integral, error1 = scipy.integrate.quad(func1,0,1,args=(p,z))
        func2_integral, error2 = scipy.integrate.quad(func2,0,1,args=(p,z))
        unobs_obsflux = func1_integral/func2_integral 
        flux_ratio.append(unobs_obsflux)
    return(flux_ratio)

'''Compute chi^2'''
def chisq(flux,flux_model,flux_error):
    chisq = sum(((flux - flux_model)/flux_error)**2)
    return(chisq)

#Values given- tau=0.1,t0=124.51,p=0.0775
flux_ratio = model_flux(time,0.1,124.51,0.0775)
chi_2 = chisq(norm_flux,flux_ratio,norm_flux_err)
print("Chi Squared:",chi_2)

'''*********************
Problem 2
Plot data, find if chi^2 value is a good fit 
**********************'''

#Plot Transit Curve
plt.errorbar(time, norm_flux, norm_flux_err,label="Kepler Data") #data 
plt.plot(time,flux_ratio,label="Predicted Transit") #simulated fit 
plt.title("Transit Curve - Problem 2")
plt.ylabel('Flux',size=12)
plt.xlabel('Time',size=12)
plt.legend()
plt.savefig("Problem2.png",dpi=200)
#
#Plot p-value
def degfreedom(N,M):
    v = N - M
    return(v)

dof = degfreedom(len(norm_flux),3)
p_value = scipy.special.gammaincc(dof/2,chi_2/2)
print("P-value:",p_value)

'''*********************
Problem 3
Vary tau in calculating chi squared
**********************'''

#Calculate model and respective chi squared 
fluxmodel_values = 0
lowest_tau = 0
chivalue = 10000
tau_list   = []
chisq_list = []
for tau in np.arange(0.08,0.13,0.001):
    fluxv = model_flux(time,tau,124.51,0.0775)
    chi_2 = chisq(norm_flux,fluxv,norm_flux_err)
    tau_list.append(tau)
    chisq_list.append(chi_2)
    if(chi_2<chivalue):
        chivalue = chi_2
        lowest_tau = tau 
        fluxmodel_values = fluxv

print("Min Chi Squared value:",chivalue)
print("Corresponding min tau:",lowest_tau)

#Plot Results
fig = plt.figure()
ax = plt.subplot(111)
ax.errorbar(time, norm_flux, norm_flux_err,label="Kepler Data") #data 
ax.plot(time,fluxmodel_values,label="Predicted Transit") #simulated fit 
plt.title("Transit Curve - Problem 3")
plt.ylabel('Flux',size=12)
plt.xlabel('Time',size=12)
plt.legend()
plt.savefig("Problem2.png",dpi=200)

#Calculate p-value
dof = degfreedom(len(fluxmodel_values),3)
p_value = scipy.special.gammaincc(dof/2,chivalue/2)
print("P-value, lowest chi squared:",p_value)

#Plot chi-squared plot
fig1 = plt.figure()
ax1 = plt.subplot(111)
ax1.plot(tau_list,chisq_list) 
plt.title("Chi-Squared vs. Tau Plot")
plt.ylabel("Chi Squared",size=12)
plt.xlabel("Tau",size=12)
plt.legend()
plt.savefig("Problem2_Chi2.png",dpi=200)

#chi squared min value plus 1, to find sigma 1 tau 
chisq_newtau = chivalue + 1
ax1.axhline(chisq_newtau)



    








