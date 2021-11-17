#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copied from https://github.com/fjansson/cloudbotany/blob/cloud_botany/thermo.py

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from functions import lowPass, mean_mask, getRad
from thermofunctions import *

# Python version of the thermodynamic calculations in DALES

# Inputs
itmin = 23
itmax = 24
di    = 1
izmin = 0
izmax = 80
klp = 4
qlmin = 1e-10
lp = '/scratch-shared/janssens/bomex200_e12'

ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')
ds0= nc.Dataset(lp+'/tmser.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])
rhobh = np.ma.getdata(ds1.variables['rhobh'][:])

plttime = np.arange(itmin, itmax, di)
zflim = zf[izmin:izmax]

# Mask for low-[ass filtering
circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

for i in range(len(plttime)):
    it1d = np.argmin(np.abs(time1d/3600 - time[plttime[i]]))
    
    # 1D fields
    rhobfi = rhobf[it1d,izmin:izmax]
    rhobhi = rhobh[it1d,izmin:izmax]
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    # wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    # u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    # v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    
    # Compute T
    presh  = np.ma.getdata(ds1.variables['presh'][it1d,izmin:izmax])
    presf  = (presh[1:]+presh[:-1])*0.5
    exnf   = (presf/1e5)**(rd/cp)
    T = exnf[:,np.newaxis,np.newaxis]*thl[:-1,:,:] + rlv/cp*ql[:-1,:,:]
    
    # Compute qs
    qs = qsatur(T,presf[:,np.newaxis,np.newaxis])

    # Test qs -> Should give the same ql
    qltest = np.maximum(qt[:-1,:,:] - qs, 0.)
    plt.plot(np.mean(qltest,axis=(1,2)),zflim[:-1])
    plt.plot(np.mean(ql,axis=(1,2)),zflim)
    
    # Clausius-Claperyon equation
    qst = rlv*qs/rv/T**2
    
    # Slab-mean, fluctuation
    cm = ql.copy()
    cm[ql>=qlmin] = 1.
    cm[ql<qlmin] = 0.
    cf_av = np.mean(cm,axis=(1,2))
    
    ql_av = np.mean(ql,axis=(1,2))
    qlp = ql - ql_av[:,np.newaxis,np.newaxis]
    
    qs_av = np.mean(qs,axis=(1,2))
    qsp = qs - qs_av[:,np.newaxis,np.newaxis]
    
    thl_av = np.mean(thl,axis=(1,2))
    thlp = thl - thl_av[:,np.newaxis,np.newaxis]
    
    qt_av = np.mean(qt,axis=(1,2))
    qtp = qt - qt_av[:,np.newaxis,np.newaxis]
    
    T_av = np.mean(T,axis=(1,2))
    
    # Cloud-conditioned average and fluctuation
    cmnan = cm.copy()
    cmnan[cmnan==0.] = np.nan
    
    # a-verged over c-louds
    thlac = np.nanmean(thl*cmnan,axis=(1,2))
    qlac = np.nanmean(ql*cmnan,axis=(1,2))
    qtac = np.nanmean(qt*cmnan,axis=(1,2))
    qsac = np.nanmean(qs*cmnan[:-1,:,:],axis=(1,2))
    
    # f-luctuation over c-louds
    thlpc = cmnan*thl - thlac[:,np.newaxis,np.newaxis]
    qlpc = cmnan*ql - qlac[:,np.newaxis,np.newaxis]
    qtpc = cmnan*qt - qtac[:,np.newaxis,np.newaxis]
    qspc = cmnan[:-1,:,:]*qs - qsac[:,np.newaxis,np.newaxis]
    
    # Decomposition of qlp in cloudy regions into contributions from qt and thl
    qlpc_thl = -exnf[:,np.newaxis,np.newaxis]*qst/(rlv/cp*qst + 1)*thlpc[:-1,:,:]
    qlpc_qt = 1/(rlv/cp*qst + 1)*qtpc[:-1,:,:]
    
    # And using xpc = xp*cmnan + x_av*(cf_av-1)/cf_av
    qlpc_thl = qlpc_thl = -exnf[:,np.newaxis,np.newaxis]*qst/(rlv/cp*qst + 1)*((thlp*cmnan)[:-1,:,:] + (thl_av-thlac)[:-1,np.newaxis,np.newaxis])
    qlpc_qt = qlpc_qt = 1/(rlv/cp*qst + 1)*((qtp*cmnan)[:-1,:,:] + (qt_av-qtac)[:-1,np.newaxis,np.newaxis])

    # Moist/dry averaging, over the large/small scales AND clouds
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qt[:,:,:],zflim,axis=0)
    twp = lowPass(twp, circ_mask)
    mask_moist = np.zeros(twp.shape)
    mask_moist[twp - np.mean(twp) > 0] = 1
    mask_dry = 1 - mask_moist
    
    mask_moist_cl = cm.copy()
    for k in range(len(zflim)):
        mask_moist_cl[k,mask_moist==0.] = 0.
    mask_moist_cl[mask_moist_cl==0.] = np.nan
    
    cf_moist = np.nansum(mask_moist_cl,axis=(1,2))/np.sum(mask_moist)
    
    ql_moist_cl = np.nanmean(mask_moist_cl*ql,axis=(1,2))
    qt_moist_cl = np.nanmean(mask_moist_cl*qt,axis=(1,2))
    qs_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qs,axis=(1,2))
    
    # Test: qlp_moist = ql_moist_cl*cf_moist - ql_av
    #                 = (qt_moist_cl-qs_moist_cl)*cf_moist - (qt_cl - qs_cl)*cf_av
    plt.figure()
    plt.plot(mean_mask(qlp, mask_moist),zflim)
    plt.plot(ql_moist_cl*cf_moist - ql_av,zflim)
    plt.plot((qt_moist_cl[:-1]-qs_moist_cl)*cf_moist[:-1] - 
             (qtac[:-1]-qsac)*cf_av[:-1],zflim[:-1])
    
    
        
    # qlp_thl = -qst*thlp[:-1,:,:]
    # qlp_qt  = 1   /(rlv/cp*qst + 1)*qtp[:-1,:,:]
    
    # qlpthl_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qlp_thl,axis=(1,2))
    # qlpqt_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qlp_qt,axis=(1,2))
    # qlp_moist_cl = np.nanmean(mask_moist_cl*qlp,axis=(1,2))
    
    
