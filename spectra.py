#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:23:26 2021

@author: janssens
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import netCDF4 as nc
from skimage.measure import block_reduce
import gc
sys.path.insert(1, '/home/janssens/scripts/pp3d/')
from functions import *

# Run specifics
itmin = 63#23
itmax = 67
di    = 1
izmin = 39
izmax = 40

klp = 4

lp = '/scratch-shared/janssens/bomex100_e12'
ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')
ilp = np.loadtxt(lp+'/lscale.inp.001')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])

dx = xf[1] - xf[0]
dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain

plttime = np.arange(itmin, itmax, di)

def plot_spectrum(k1d, spec, lab, plttime):
    fig = plt.figure(); ax = plt.gca()
    for i in range(len(plttime)):
        col = plt.cm.cubehelix(i/len(plttime))
        ax.loglog(k1d,spec[i,izpl,:],c=col,label='t=%.2f'%time[plttime[i]])
    # ax.set_ylim((1e-4,1e2))
    ax.set_ylabel(lab)
    ax.set_xlabel(r"Wavenumber [1/m]")
    ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)
    
    ax2 = ax.twiny()
    fig.subplots_adjust(bottom=0.22)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.spines['bottom'].set_position(('axes',-0.22))
    ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
    ax2.set_xscale('log')
    ax2.set_xlabel('Wavelength [m]')
    plt.show()

#%% Compute spectra at a given height

zflim = zf[izmin:izmax]
N = xf.size; N2 = N//2
spec_qt = np.zeros((len(plttime),len(zflim),N2))
spec_thl = np.zeros((len(plttime),len(zflim),N2))
spec_thlv = np.zeros((len(plttime),len(zflim),N2))
spec_w = np.zeros((len(plttime),len(zflim),N2))
spec_ql = np.zeros((len(plttime),len(zflim),N2))
spec_wqt = np.zeros((len(plttime),len(zflim),N2))
spec_wthl = np.zeros((len(plttime),len(zflim),N2))
spec_wthlv = np.zeros((len(plttime),len(zflim),N2))

for i in range(len(plttime)):
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    # u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    # v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    
    thl_av = np.mean(thl,axis=(1,2))
    thlv = thl + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt

    thl = thl - thl_av
    thlv = thlv - np.mean(thlv,axis=(1,2))
    qt = qt - np.mean(qt,axis=(1,2))
    ql = ql - np.mean(qt,axis=(1,2))
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    del wh
    
    for iz in range(len(zflim)):
        k1d,spec_qt[i,iz,:] = compute_spectrum(qt[iz,:,:], dx)
        k1d,spec_thl[i,iz,:] = compute_spectrum(thl[iz,:,:], dx)
        k1d,spec_thlv[i,iz,:] = compute_spectrum(thlv[iz,:,:], dx)
        k1d,spec_w[i,iz,:] = compute_spectrum(wf[iz,:,:], dx)
        k1d,spec_ql[i,iz,:] = compute_spectrum(ql[iz,:,:], dx)
        
        k1d,spec_wqt[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, qt[iz,:,:])
        k1d,spec_wthl[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, thl[iz,:,:])
        k1d,spec_wthlv[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, thlv[iz,:,:])
        
        # k1d,spec_wqt[i,iz,:] = compute_spectrum(wf[iz,:,:]*qt[iz,:,:], dx,sqrt=True)
        # k1d,spec_wthl[i,iz,:] = compute_spectrum(wf[iz,:,:]*thl[iz,:,:], dx,sqrt=True)
        # k1d,spec_wthlv[i,iz,:] = compute_spectrum(wf[iz,:,:]*thlv[iz,:,:], dx,sqrt=True)
        
    gc.collect()

#%% Average over time
itav = 4 # number of time steps to average over -> len(plttime) MUST BE MULTIPLE OF THIS

spec_qt_mn = block_reduce(spec_qt,(itav,1,1),func=np.mean)
spec_thl_mn = block_reduce(spec_thl,(itav,1,1),func=np.mean)
spec_thlv_mn = block_reduce(spec_thlv,(itav,1,1),func=np.mean)
spec_w_mn = block_reduce(spec_w,(itav,1,1),func=np.mean)
spec_ql_mn = block_reduce(spec_ql,(itav,1,1),func=np.mean)

spec_wqt_mn = block_reduce(spec_wqt,(itav,1,1),func=np.mean)
spec_wthl_mn = block_reduce(spec_wthl,(itav,1,1),func=np.mean)
spec_wthlv_mn = block_reduce(spec_wthlv,(itav,1,1),func=np.mean)

plttime_mn = plttime[::itav]

#%% Plot
izpl = 0

# Variances
plot_spectrum(k1d, spec_qt_mn, r"$k\widehat{q}_t'^2$", plttime_mn)
plot_spectrum(k1d, spec_thl_mn, r"$k\widehat{\theta}_l'^2$", plttime_mn)
plot_spectrum(k1d, spec_thlv_mn, r"$k\widehat{\theta}_{lv}'^2$", plttime_mn)
plot_spectrum(k1d, spec_w_mn, r"$k\widehat{w}'^2$", plttime_mn)
plot_spectrum(k1d, spec_ql_mn, r"$k\widehat{q_l}'^2$", plttime_mn)

# Fluxes
plot_spectrum(k1d, spec_wqt_mn, r"$k\widehat{wq}_t'$", plttime_mn)
plot_spectrum(k1d, spec_wthl_mn, r"$k\widehat{w\theta}_l'$", plttime_mn)
plot_spectrum(k1d, spec_wthlv_mn, r"$k\widehat{w\theta}_{lv}'$", plttime_mn)

#%% Plot in same spectrum FIXME is not yet implemented

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_wthlv_mn[-1,izpl,:],label=r"$w\theta_{lv}$")
ax.loglog(k1d,spec_wthlv_t_mn[-1,izpl,:],label=r"$w\theta_{lv}'$")
ax.loglog(k1d,spec_wthlv_r_mn[-1,izpl,:],label=r"$w'''\theta_{lv}'''$")
# ax.set_ylim((1e-4,1e2))
ax.set_ylabel(r"$k\widehat{w\theta}_{lv}'$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()

#%% wthlv scale decomposition (run cells above first)
klp = 4

spec_wthlv_l = np.zeros((len(plttime),len(zflim),N2))
spec_wthlv_c = np.zeros((len(plttime),len(zflim),N2))
spec_wthlv_r = np.zeros((len(plttime),len(zflim),N2))

# Mask for low-[ass filtering
circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

for i in range(len(plttime)):
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    
    thl_av = np.mean(thl,axis=(1,2))
    thlv = thl + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt
    thlv = thlv - np.mean(thlv,axis=(1,2))
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    del wh

    # Low-pass filter (and identify high-pass filtered remainder)    
    wff = lowPass(wf, circ_mask)
    wfp = wf - wff
    del wf
                
    thlvf = lowPass(thlv, circ_mask)
    thlvp = thlv - thlvf
    del thlv
    
    gc.collect()
    
    for iz in range(len(zflim)):
        # k1d,spec_wthlv_l[i,iz,:] = compute_spectrum(wff[iz,:,:]*thlvf[iz,:,:], dx)
        # k1d,spec_wthlv_c[i,iz,:] = compute_spectrum(wff[iz,:,:]*thlvp[iz,:,:]+
        #                                             wfp[iz,:,:]*thlvf[iz,:,:], dx)
        # k1d,spec_wthlv_r[i,iz,:] = compute_spectrum(wfp[iz,:,:]*thlvp[iz,:,:], dx)

        k1d,spec_wthlv_l[i,iz,:] = compute_spectrum(wff[iz,:,:], dx,thlvf[iz,:,:])
        k1d,spec_wthlv_c[i,iz,:] = compute_spectrum(wff[iz,:,:], dx, thlvp[iz,:,:])
        _,spec_wthlv_c2 = compute_spectrum(wfp[iz,:,:], dx, thlvf[iz,:,:])
        spec_wthlv_c[i,iz,:] += spec_wthlv_c2
        k1d,spec_wthlv_r[i,iz,:] = compute_spectrum(wfp[iz,:,:], dx, thlvp[iz,:,:])


spec_wthlv_l_mn = block_reduce(spec_wthlv_l,(itav,1,1),func=np.mean)
spec_wthlv_c_mn = block_reduce(spec_wthlv_c,(itav,1,1),func=np.mean)
spec_wthlv_r_mn = block_reduce(spec_wthlv_r,(itav,1,1),func=np.mean)
sumtest = spec_wthlv_l_mn + spec_wthlv_c_mn + spec_wthlv_r_mn

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_wthlv_mn[-1,izpl,:],label=r"$w'\theta_{lv}'$")
ax.loglog(k1d,spec_wthlv_l_mn[-1,izpl,:],label=r"$\widetilde{w'}\widetilde{\theta_{lv}'}$")
ax.loglog(k1d,spec_wthlv_c_mn[-1,izpl,:],label=r"$\widetilde{w'}\theta_{lv}'''+w'''\widetilde{\theta_{lv}'}$")
ax.loglog(k1d,spec_wthlv_r_mn[-1,izpl,:],label=r"$w'''\theta_{lv}'''$")
# ax.loglog(k1d,sumtest[-1,izpl,:],label=r"sum")
# ax.set_ylim((1e-4,1e2))
ax.set_ylabel(r"$k\widehat{w\theta}_{lv}'$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()
