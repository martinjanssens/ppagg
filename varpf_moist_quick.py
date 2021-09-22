#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 11:07:59 2021

@author: janssens
"""


import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import gc
import sys
sys.path.insert(1, '/home/janssens/scripts/pp3d/')
from functions import *

itmin = 11
itmax = 33
di    = 1
izmin = 0
izmax = 80

klp = 4

lp = '/scratch-shared/janssens/bomex100'
ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])

plttime = np.arange(itmin, itmax, di)
zflim = zf[izmin:izmax]

#%% Dry/moist regions

qtpf_moist_time = np.zeros((plttime.size,izmax-izmin))
thlvpf_moist_time = np.zeros((plttime.size,izmax-izmin))
wf_moist_time = np.zeros((plttime.size,izmax-izmin))

# Mask for low-[ass filtering
circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

for i in range(len(plttime)):
    it1d = np.argmin(np.abs(time1d/3600 - time[plttime[i]]))
    
    # 1D fields
    rhobfi = rhobf[it1d,izmin:izmax]
    
    # 3D fields
    qtp  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thlp =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    # qlp = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    # u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    # v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    # e12 = np.ma.getdata(ds.variables['e12'][plttime[i],izmin:izmax,:,:])

    # Slab averaging
    thl_av = np.mean(thlp,axis=(1,2))
    qt_av = np.mean(qtp,axis=(1,2))   
    # ql_av = np.mean(qlp,axis=(1,2))

    # thlv
    thlvp = thlp + 0.608*thl_av[:,np.newaxis,np.newaxis]*qtp
    thlv_av = np.mean(thlvp,axis=(1,2))
    
    # twp
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qtp[:,:,:],zflim,axis=0)
    
    # subtract mean
    thlvp = thlvp - thlv_av[:,np.newaxis,np.newaxis]
    qtp = qtp - qt_av[:,np.newaxis,np.newaxis]
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    del wh

    # Low-pass filter (and identify high-pass filtered remainder)    
               
    thlvpf = lowPass(thlvp, circ_mask)
    thlvpp = thlvp - thlvpf
    del thlvp

    qtpf = lowPass(qtp, circ_mask)
    qtpp = qtp - qtpf
    del qtp
    
    wff = lowPass(wf, circ_mask)
    wfp = wf - wff
    del wf
    
    # Moist/dry averaging, over the large/small scales
    twp = lowPass(twp, circ_mask)
    mask_moist = np.zeros(twp.shape)
    mask_moist[twp - np.mean(twp) > 0] = 1
    mask_dry = 1 - mask_moist
    
    # Moist, large
    thlvpf_moist_time[i,:] = mean_mask(thlvpf,mask_moist)
    qtpf_moist_time[i,:] = qtpf_moist = mean_mask(qtpf,mask_moist)
    wf_moist_time[i,:] = mean_mask(wff,mask_moist)
    # w_moist_h = mean_mask(whf,mask_moist)
    # qlpf_moist = mean_mask(qlpf,mask_moist)
    
    gc.collect()
    
#%% Plot it

fig,axs = plt.subplots(ncols=3,sharey=True,figsize=(6,5))
for i in range(len(plttime)):
    col = plt.cm.cubehelix(i/len(plttime))
     
    axs[0].plot(qtpf_moist_time[i,:], zflim, color=col,linestyle='-')
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].set_xlabel(r"$\widetilde{q_t'}$")
    axs[0].set_xlim((0,6e-4))
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[1].plot(thlvpf_moist_time[i,:], zflim, color=col,linestyle='-')
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].set_xlabel(r"$\widetilde{\theta_{lv}'}$")
    axs[1].set_xlim((-4e-2,1e-2))
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[2].plot(wf_moist_time[i,:], zflim,color=col,linestyle='-',label='t=%.2f'%time[plttime[i]])
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].set_xlabel(r"$\widetilde{w'}$")
    axs[2].set_xlim((-1e-2,1.7e-2))
    axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[2].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

