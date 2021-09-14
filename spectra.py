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
itmin = 59#23
itmax = 63
di    = 1
izmin = 39
izmax = 40

klp = 4

lp = '/scratch-shared/janssens/bomex200_e12'
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
spec_wqt = np.zeros((len(plttime),len(zflim),N2))
spec_wthl = np.zeros((len(plttime),len(zflim),N2))
spec_wthlv = np.zeros((len(plttime),len(zflim),N2))

for i in range(len(plttime)):
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    # qlp = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    # u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    # v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    
    thl_av = np.mean(thl,axis=(1,2))
    thlv = thl + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    del wh
    
    for iz in range(len(zflim)):
        k1d,spec_qt[i,iz,:] = compute_spectrum(qt[iz,:,:], dx)
        k1d,spec_thl[i,iz,:] = compute_spectrum(thl[iz,:,:], dx)
        k1d,spec_thlv[i,iz,:] = compute_spectrum(thlv[iz,:,:], dx)
        
        k1d,spec_wqt[i,iz,:] = compute_spectrum(wf[iz,:,:]*qt[iz,:,:], dx)
        k1d,spec_wthl[i,iz,:] = compute_spectrum(wf[iz,:,:]*thl[iz,:,:], dx)
        k1d,spec_wthlv[i,iz,:] = compute_spectrum(wf[iz,:,:]*thlv[iz,:,:], dx)        
    
    gc.collect()

#%% Average over time
itav = 4 # number of time steps to average over -> MUST BE MULTIPLE OF len(plttime)

spec_qt_mn = block_reduce(spec_qt,(itav,1,1),func=np.mean)
spec_thl_mn = block_reduce(spec_thl,(itav,1,1),func=np.mean)
spec_thlv_mn = block_reduce(spec_thlv,(itav,1,1),func=np.mean)

spec_wqt_mn = block_reduce(spec_wqt,(itav,1,1),func=np.mean)
spec_wthl_mn = block_reduce(spec_wthl,(itav,1,1),func=np.mean)
spec_wthlv_mn = block_reduce(spec_wthlv,(itav,1,1),func=np.mean)

plttime_mn = plttime[::itav]

#%% Plot
izpl = 0

plot_spectrum(k1d, spec_qt_mn, r"$\hat{q}_t'^2$", plttime_mn)
plot_spectrum(k1d, spec_thl_mn, r"$\hat{\theta}_l'^2$", plttime_mn)
plot_spectrum(k1d, spec_thlv_mn, r"$\hat{\theta}_{lv}'^2$", plttime_mn)

plot_spectrum(k1d, spec_wqt_mn, r"$\hat{wq}_t'^2$", plttime_mn)
plot_spectrum(k1d, spec_wthl_mn, r"$\hat{w\theta}_l'^2$", plttime_mn)
plot_spectrum(k1d, spec_wthlv_mn, r"$\hat{w\theta}_{lv}'^2$", plttime_mn)
