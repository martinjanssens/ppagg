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
itmin = 63
itmax = 64
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

Ms = np.zeros((len(plttime),zflim.size))
wqt_mf = np.zeros((len(plttime),zflim.size))
wthlv_mf = np.zeros((len(plttime),zflim.size))
wqt_av_time = np.zeros((len(plttime),zflim.size))
wthlv_av_time = np.zeros((len(plttime),zflim.size))
for i in range(len(plttime)):
    it1d = np.argmin(np.abs(time1d/3600 - time[plttime[i]]))
    
    # 1D fields
    rhobfi = rhobf[it1d,izmin:izmax]
    rhobhi = rhobh[it1d,izmin:izmax]
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    buoy =  np.ma.getdata(ds.variables['buoy'][plttime[i],izmin:izmax,:,:])
    # u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    # v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    
    # thlv
    thlv = thl + (rv/rd-1)*qt
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    wh = wh[:-1,:,:]
    
    # Moist/dry definition
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qt[:,:,:],zflim,axis=0)
    twp = lowPass(twp, circ_mask)
    mask_moist = np.zeros(twp.shape)
    mask_moist[twp - np.mean(twp) > 0] = 1
    mask_dry = 1 - mask_moist
    
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
    
    thlv_av = np.mean(thl,axis=(1,2))
    thlvp = thlv - thlv_av[:,np.newaxis,np.newaxis]
    
    qt_av = np.mean(qt,axis=(1,2))
    qtp = qt - qt_av[:,np.newaxis,np.newaxis]
    qtpf = lowPass(qtp,circ_mask)
    
    T_av = np.mean(T,axis=(1,2))
    
    # Cloud-conditioned average and fluctuation
    cmnan = cm.copy()
    cmnan[cmnan==0.] = np.nan
    
    # Core    
    cmcore = cmnan.copy()
    cmcore[buoy<=0] = np.nan
    cf_core = np.nansum(cmcore,axis=(1,2))/xf.size**2
    
    # a-verged over c-louds
    thlca = np.nanmean(thl*cmnan,axis=(1,2))
    thlvca = np.nanmean(thlv*cmnan,axis=(1,2))
    qlca = np.nanmean(ql*cmnan,axis=(1,2))
    qtca = np.nanmean(qt*cmnan,axis=(1,2))
    qsca = np.nanmean(qs*cmnan[:-1,:,:],axis=(1,2))
    
    # Slab-Average
    # a-veraged over c-loud c-ores
    wcca = np.nanmean(wf*cmcore,axis=(1,2))
    qtcca = np.nanmean(qt*cmcore,axis=(1,2))
    thlvcca = np.nanmean(thlv*cmcore,axis=(1,2))
    
    # Mass flux
    Ms[i,:] = wcca*cf_core
    wqt_mf[i,:] = wcca*cf_core*(qtcca-qt_av)
    wthlv_mf[i,:] = wcca*cf_core*(thlvcca - thlv_av)

    wqt_av_time[i,:] = np.mean(wf*qtp,axis=(1,2))
    wthlv_av_time[i,:] = np.mean(wf*thlvp,axis=(1,2))

    # Moist
    cmc_moist = cmcore.copy()
    cmc_moist[:,mask_moist==0.] = np.nan
    cfc_moist = np.nansum(cmc_moist,axis=(1,2))/np.sum(mask_moist)
    wcc_moist = np.nanmean(wf*cmc_moist,axis=(1,2))
    qtcc_moist = np.nanmean(qt*cmc_moist,axis=(1,2))
    thlvcc_moist = np.nanmean(thlv*cmc_moist,axis=(1,2))
    qt_moist = mean_mask(qt,mask_moist)
    thlv_moist = mean_mask(thlv,mask_moist)
    Ms_moist = wcc_moist*cfc_moist
    wqt_mf_moist = Ms_moist*(qtcc_moist - qt_moist)
    wthlv_mf_moist = Ms_moist*(thlvcc_moist - thlv_moist)    
    
    wqt_moist = mean_mask(wf*qtp,mask_moist)
    wthlv_moist = mean_mask(wf*thlvp,mask_moist)
    
    # # p-erturbation from c-loud average
    # thlcp = cmnan*thl - thlca[:,np.newaxis,np.newaxis]
    # thlvcp = cmnan*thlv - thlvca[:,np.newaxis,np.newaxis]
    # qlcp = cmnan*ql - qlca[:,np.newaxis,np.newaxis]
    # qtcp = cmnan*qt - qtca[:,np.newaxis,np.newaxis]
    # qscp = cmnan[:-1,:,:]*qs - qsca[:,np.newaxis,np.newaxis]
    
    # # Decomposition of qlp in cloudy regions into contributions from qt and thl    
    # # Using:
    # # - qlcp = qtcp - qscp
    # # - qscp = qst*Tcp
    # # - thlcp = 1/exnf*(Tcp - rlv/cp*qlcp)
    # # - xcp = xp*cmnan + x_av - xac
    # qlcp_thl = -exnf[:,np.newaxis,np.newaxis]*qst/(rlv/cp*qst + 1)*((thlp*cmnan)[:-1,:,:] + (thl_av-thlca)[:-1,np.newaxis,np.newaxis])
    # qlcp_qt = 1/(rlv/cp*qst + 1)*((qtp*cmnan)[:-1,:,:] + (qt_av-qtca)[:-1,np.newaxis,np.newaxis])

    # #  Perturbations from slab-average in cloudy cells
    # qlpc_mod = qlcp_thl + qlcp_qt - (ql_av + qlca)[:-1,np.newaxis,np.newaxis]
    
    # plt.figure()
    # plt.plot(np.nanmean((cmnan*qlp)**2,axis=(1,2)),zflim)
    # plt.plot(np.nanmean(qlpc_mod**2,axis=(1,2)),zflim[:-1])
    # plt.plot(np.nanmean((qlpc_mod-qlcp_thl)**2,axis=(1,2)),zflim[:-1])
    # plt.plot(np.nanmean((qlpc_mod-qlcp_qt)**2,axis=(1,2)),zflim[:-1])
    # plt.plot(np.nanmean((qlpc_mod+(ql_av + qlca)[:-1,np.newaxis,np.newaxis])**2,axis=(1,2)),zflim[:-1])

    # # Model in terms of thlv
    # qlcp_thlv = -exnf[:,np.newaxis,np.newaxis]*qst/(rlv/cp*qst + 1)*((thlvp*cmnan)[:-1,:,:] + (thlv_av-thlvca)[:-1,np.newaxis,np.newaxis])
    # qlcp_qt1 = (1+(rd/rv-1)*exnf[:,np.newaxis,np.newaxis]*qst)/(rlv/cp*qst + 1)*((qtp*cmnan)[:-1,:,:] + (qt_av-qtca)[:-1,np.newaxis,np.newaxis])
    # qlpc_mod1 = qlcp_thlv + qlcp_qt1 - (ql_av + qlca)[:-1,np.newaxis,np.newaxis]

    # plt.figure()
    # plt.plot(np.nanmean((cmnan*qlp)**2,axis=(1,2)),zflim)
    # plt.plot(np.nanmean(qlpc_mod1**2,axis=(1,2)),zflim[:-1])
    # plt.plot(np.nanmean((qlpc_mod-qlcp_thlv)**2,axis=(1,2)),zflim[:-1])
    
    # for k in range(len(zflim)-1):
    #     qlpc_mod1[k,np.isnan(qlpc_mod1[k,:,:])] = -ql_av[k]
    
    # plt.figure()
    # plt.plot(mean_mask(cm*qlp,mask_moist),zflim)
    # plt.plot(mean_mask(qlpc_mod1,mask_moist),zflim[:-1])
    # # plt.plot(np.nanmean((qlpc_mod-qlcp_thlv)**2,axis=(1,2)),zflim[:-1])
    # ## This is still not working.
    
    # # Mesoscale fluxes  
    
    # wqlp = wf*qlp
    # wqlpc_mod = wf[:-1,:,:]*qlpc_mod1
    # wqlpc_qt1 = wf[:-1,:,:]*qlcp_qt1
    
    # # Test
    # wqlpf = lowPass(wqlp,circ_mask)
    # wqlpf_moist = mean_mask(wqlpf,mask_moist)
        
    # wqlpc_modf = lowPass(wqlpc_mod,circ_mask)
    # wqlpc_modf_moist = mean_mask(wqlpc_modf,mask_moist)
    
    # wqlpc_qt1[np.isnan(wqlpc_qt1)] = 0.
    # wqlpc_qtf = lowPass(wqlpc_qt1,circ_mask)
    # wqlpc_qtf_moist = mean_mask(wqlpc_qtf,mask_moist)
    

    # # Moist/dry averaging, over the large/small scales AND cloud
    
    # mask_moist_cl = cm.copy()
    # for k in range(len(zflim)):
    #     mask_moist_cl[k,mask_moist==0.] = 0.
    # mask_moist_cl[mask_moist_cl==0.] = np.nan
    
    # cf_moist = np.nansum(mask_moist_cl,axis=(1,2))/np.sum(mask_moist)
    
    # ql_moist_cl = np.nanmean(mask_moist_cl*ql,axis=(1,2))
    # qt_moist_cl = np.nanmean(mask_moist_cl*qt,axis=(1,2))
    # qs_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qs,axis=(1,2))
    
    # # Test: qlp_moist = ql_moist_cl*cf_moist - ql_av
    # #                 = (qt_moist_cl-qs_moist_cl)*cf_moist - (qt_cl - qs_cl)*cf_av
    # plt.figure()
    # plt.plot(mean_mask(qlp, mask_moist),zflim)
    # plt.plot(ql_moist_cl*cf_moist - ql_av,zflim)
    # plt.plot((qt_moist_cl[:-1]-qs_moist_cl)*cf_moist[:-1] - 
    #          (qtac[:-1]-qsac)*cf_av[:-1],zflim[:-1])
    
    
        
    # qlp_thl = -qst*thlp[:-1,:,:]
    # qlp_qt  = 1   /(rlv/cp*qst + 1)*qtp[:-1,:,:]
    
    # qlpthl_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qlp_thl,axis=(1,2))
    # qlpqt_moist_cl = np.nanmean(mask_moist_cl[:-1,:,:]*qlp_qt,axis=(1,2))
    # qlp_moist_cl = np.nanmean(mask_moist_cl*qlp,axis=(1,2))
    
    
