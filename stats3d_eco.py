#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 20:24:33 2020

@author: martinjanssens
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import gc
import sys
sys.path.insert(1, '/home/janssens/scripts/pp3d/')
from functions import *
import argparse

parseFlag = False

if parseFlag:
    parser = argparse.ArgumentParser(description="Merge cross-section and field dump DALES output from parallel runs")
    parser.add_argument("--dir", metavar="DIR", type=str, default=".", help="Directory to load/store data from/to")
    parser.add_argument("--itmin", metavar="N", type=int, default=0, help="First time index")
    parser.add_argument("--itmax", metavar="N", type=int, default=-1, help="Last time index")
    parser.add_argument("--dt", metavar="N", type=int, default=1, help="Time sampling interval")
    parser.add_argument("--izmin", metavar="N", type=int, default=0, help="First height index")
    parser.add_argument("--izmax", metavar="N", type=int, default=80, help="Last height index")
    parser.add_argument("--klp", metavar="N", type=int, default=4, help="Cutoff wavenumber for lw-pass filter")
    parser.add_argument("--store", action="store_true", default=False, help="Process only fielddump")

    args = parser.parse_args()

    lp = args.dir
    itmin = args.itmin
    itmax = args.itmax
    di = args.dt
    izmin = args.izmin
    izmax = args.izmax
    klp = args.klp
    store = args.store

else:
    lp = '/scratch-shared/janssens/bomex200_e12'

ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')
ds0= nc.Dataset(lp+'/tmser.001.nc')
ilp = np.loadtxt(lp+'/lscale.inp.001')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])
rhobh = np.ma.getdata(ds1.variables['rhobh'][:])

dx = np.diff(xf)[0]
dy = np.diff(yf)[0] # Assumes uniform horizontal spacing
dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain

delta = (dx*dy*np.diff(zh))**(1./3)

# Larger-scale subsidence
wfls = ilp[:,3]

#%% Dry/moist regions

if not parseFlag:
    itmin = 46#23
    itmax = 47
    di    = 1
    izmin = 0
    izmax = 80
    store = False
    klp = 4

plttime = np.arange(itmin, itmax, di)
zflim = zf[izmin:izmax]

qtpf_moist_time = np.zeros((plttime.size,izmax-izmin))
qtpf_dry_time = np.zeros((plttime.size,izmax-izmin))
qtpf_prod_moist_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_prod_dry_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_prod_moist_wex_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_prod_dry_wex_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_vdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_vdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_hdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_hdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_subs_moist_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_subs_dry_time = np.zeros((plttime.size,izmax-izmin-2))
qtpf_diff_moist_time = np.zeros((plttime.size,izmax-izmin-4))
qtpf_diff_dry_time = np.zeros((plttime.size,izmax-izmin-4))

thlvpf_moist_time = np.zeros((plttime.size,izmax-izmin))
thlvpf_dry_time = np.zeros((plttime.size,izmax-izmin))
thlvpf_prod_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_prod_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_vdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_vdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_hdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_hdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_subs_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_subs_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpf_diff_moist_time = np.zeros((plttime.size,izmax-izmin-4))
thlvpf_diff_dry_time = np.zeros((plttime.size,izmax-izmin-4))

thlvpp_moist_time = np.zeros((plttime.size,izmax-izmin))
thlvpp_dry_time = np.zeros((plttime.size,izmax-izmin))
thlvpp_prod_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_prod_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_vdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_vdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_hdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_hdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_subs_moist_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_subs_dry_time = np.zeros((plttime.size,izmax-izmin-2))
thlvpp_diff_moist_time = np.zeros((plttime.size,izmax-izmin-4))
thlvpp_diff_dry_time = np.zeros((plttime.size,izmax-izmin-4))

wthlvpf_prod_moist_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_prod_dry_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_vdiv_moist_time = np.zeros((plttime.size,izmax-izmin-3))
wthlvpf_vdiv_dry_time = np.zeros((plttime.size,izmax-izmin-3))
wthlvpf_hdiv_moist_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_hdiv_dry_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_buoy_moist_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_buoy_dry_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_subs_moist_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_subs_dry_time = np.zeros((plttime.size,izmax-izmin-2))
wthlvpf_diff_moist_time = np.zeros((plttime.size,izmax-izmin-4))
wthlvpf_diff_dry_time = np.zeros((plttime.size,izmax-izmin-4))

thl_av_time = np.zeros((plttime.size,izmax-izmin))
qt_av_time = np.zeros((plttime.size,izmax-izmin))
thlv_av_time = np.zeros((plttime.size,izmax-izmin))

thlpf_moist_time = np.zeros((plttime.size,izmax-izmin))
thlpf_dry_time = np.zeros((plttime.size,izmax-izmin))
wff_moist_time = np.zeros((plttime.size,izmax-izmin))
wff_dry_time = np.zeros((plttime.size,izmax-izmin))
qlpf_moist_time = np.zeros((plttime.size,izmax-izmin))
qlpf_dry_time = np.zeros((plttime.size,izmax-izmin))

qtpp_moist_time = np.zeros((plttime.size,izmax-izmin))
qtpp_dry_time = np.zeros((plttime.size,izmax-izmin))
thlpp_moist_time = np.zeros((plttime.size,izmax-izmin))
thlpp_dry_time = np.zeros((plttime.size,izmax-izmin))
wfp_moist_time = np.zeros((plttime.size,izmax-izmin))
wfp_dry_time = np.zeros((plttime.size,izmax-izmin))
qlpp_moist_time = np.zeros((plttime.size,izmax-izmin))
qlpp_dry_time = np.zeros((plttime.size,izmax-izmin))

wthlpf_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlpf_dry_time = np.zeros((plttime.size,izmax-izmin))
wqtpf_moist_time = np.zeros((plttime.size,izmax-izmin))
wqtpf_dry_time = np.zeros((plttime.size,izmax-izmin))
wqlpf_moist_time = np.zeros((plttime.size,izmax-izmin))
wqlpf_dry_time = np.zeros((plttime.size,izmax-izmin))
wthlvp_av_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_dry_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_l_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_l_dry_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_c_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_c_dry_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_r_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlvpf_r_dry_time = np.zeros((plttime.size,izmax-izmin))
wthlvpp_moist_time = np.zeros((plttime.size,izmax-izmin))
wthlvpp_dry_time = np.zeros((plttime.size,izmax-izmin))

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
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thlp =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    qlp = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    # e12 = np.ma.getdata(ds.variables['e12'][plttime[i],izmin:izmax,:,:])

    # Slab averaging
    thl_av = np.mean(thlp,axis=(1,2))
    qt_av = np.mean(qt,axis=(1,2))   
    ql_av = np.mean(qlp,axis=(1,2))
    thlv_av = thl_av*(1+0.608*qt_av)
    
    thl_av_time[i,:] = thl_av
    qt_av_time[i,:] = qt_av
    thlv_av_time[i,:] = thlv_av
    
    # -> need presf for thv_av, taken from nearest 1d data time and half-level
    presh  = np.ma.getdata(ds1.variables['presh'][it1d,izmin:izmax])
    presf  = (presh[1:]+presh[:-1])*0.5
    exnf   = (presf/1e5)**(Rd/cp)
    thv_av = (thl_av[:-1] + (Lv/cp)*ql_av[:-1]/exnf)*(1.+(Rv/Rd-1)*qt_av[:-1] -Rv/Rd*ql_av[:-1])
    
    # Eddy diffusivities
    # dthvdz = compute_dthvdz(thlp, qt, qlp, exnf, dzh)
    # _,ekhp = compute_ek(e12, dthvdz, thv_av, delta[izmin:izmax-1])
    # ekh_av = np.mean(ekhp,axis=(1,2))
    # ekhp = ekhp - ekh_av[:,np.newaxis,np.newaxis]
    # del e12
    # del dthvdz
    
    # Define thlv
    thlvp = thlp + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt
    
    # Perturbations
    qtp = qt - qt_av[:,np.newaxis,np.newaxis]
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qt[:,:,:],zflim,axis=0)
    del qt
    
    gc.collect()
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    wh = wh[:-1,:,:]

    thlp = thlp - thl_av[:,np.newaxis,np.newaxis]
    thlvp = thlvp - thlv_av[:,np.newaxis,np.newaxis]
    qlp = qlp - ql_av[:,np.newaxis,np.newaxis]
    
    # Slab average resolved fluxes
    wqt_av = np.mean(wf*qtp,axis=(1,2))
    wthl_av = np.mean(wf*thlp,axis=(1,2))
    # wql_av = np.mean(wf*qlp,axis=(1,2))
 
    # Low-pass filter (and identify high-pass filtered remainder)    
    qtpf = lowPass(qtp, circ_mask)
    qtpp = qtp - qtpf
    del qtp
    
    whf = lowPass(wh, circ_mask)
    whp = wh - whf
    del wh
    
    wff = lowPass(wf, circ_mask)
    wfp = wf - wff
    del wf
            
    thlpf = lowPass(thlp, circ_mask)
    thlpp = thlp - thlpf
    del thlp
    
    thlvpf = lowPass(thlvp, circ_mask)
    thlvpp = thlvp - thlvpf
    del thlvp
    
    qlpf = lowPass(qlp, circ_mask)
    qlpp = qlp - qlpf
    del qlp

    gc.collect()
    

    # Moist/dry averaging, over the large/small scales
    twp = lowPass(twp, circ_mask)
    mask_moist = np.zeros(twp.shape)
    mask_moist[twp - np.mean(twp) > 0] = 1
    mask_dry = 1 - mask_moist
    
    # Moist, large
    thlpf_moist = mean_mask(thlpf,mask_moist)
    qtpf_moist = mean_mask(qtpf,mask_moist)
    wf_moist = mean_mask(wff,mask_moist)
    # w_moist_h = mean_mask(whf,mask_moist)
    qlpf_moist = mean_mask(qlpf,mask_moist)
    
    # Dry, large
    thlpf_dry = mean_mask(thlpf,mask_dry)
    qtpf_dry = mean_mask(qtpf,mask_dry)
    wf_dry = mean_mask(wff,mask_dry)
    qlpf_dry = mean_mask(qlpf,mask_dry)
    
    # Moist, small
    thlpp_moist = mean_mask(thlpp,mask_moist)
    qtpp_moist = mean_mask(qtpp,mask_moist)
    wp_moist = mean_mask(wfp,mask_moist)
    # w_moist_h = mean_mask(whf,mask_moist)
    qlpp_moist = mean_mask(qlpp,mask_moist)
    
    # Dry, small
    thlpp_dry = mean_mask(thlpp,mask_dry)
    qtpp_dry = mean_mask(qtpp,mask_dry)
    wp_dry = mean_mask(wfp,mask_dry)
    qlpp_dry = mean_mask(qlpp,mask_dry)

    # Store per variable
    thlpf_moist_time[i,:] = thlpf_moist
    thlpf_dry_time[i,:] = thlpf_dry
    thlpp_moist_time[i,:] = thlpp_moist
    thlpp_dry_time[i,:] = thlpp_dry

    qtpf_moist_time[i,:] = qtpf_moist
    qtpf_dry_time[i,:] = qtpf_dry
    qtpp_moist_time[i,:] = qtpp_moist
    qtpp_dry_time[i,:] = qtpp_dry

    thlvpf_moist_time[i,:] = thlpf_moist + 0.608*thl_av*qtpf_moist
    thlvpf_dry_time[i,:] = thlpf_dry + 0.608*thl_av*qtpf_dry
    thlvpp_moist_time[i,:] = thlpp_moist + 0.608*thl_av*qtpp_moist
    thlvpp_dry_time[i,:] = thlpp_dry + 0.608*thl_av*qtpp_dry

    wff_moist_time[i,:] = wf_moist
    wff_dry_time[i,:] = wf_dry
    wfp_moist_time[i,:] = wp_moist
    wfp_dry_time[i,:] = wp_dry
    
    qlpf_moist_time[i,:] = qlpf_moist
    qlpf_dry_time[i,:] = qlpf_dry
    qlpp_moist_time[i,:] = qlpp_moist
    qlpp_dry_time[i,:] = qlpp_dry
    
    ## Fluxes 
    # FIXME no sgs here yet!!
    wthlpf = lowPass((wff+wfp)*(thlpf+thlpp), circ_mask)
    wqtpf = lowPass((wff+wfp)*(qtpf+qtpp), circ_mask)
    wqlpf = lowPass((wff+wfp)*(qlpf+qlpp), circ_mask)

    wthlpf_moist = mean_mask(wthlpf, mask_moist)
    wthlpf_dry = mean_mask(wthlpf, mask_dry)
    
    wqtpf_moist = mean_mask(wqtpf, mask_moist)
    wqtpf_dry = mean_mask(wqtpf, mask_dry)
    
    wqlpf_moist = mean_mask(wqlpf, mask_moist)
    wqlpf_dry = mean_mask(wqlpf, mask_dry)

    wthlvp = (wff+wfp)*((thlpf+thlpp) + 0.608*thl_av[:,np.newaxis,np.newaxis]*(qtpf+qtpp))
    wthlvp_av = np.mean(wthlvp,axis=(1,2))
    wthlvpp = wthlvp - (wthlpf + 0.608*thl_av[:,np.newaxis,np.newaxis]*wqtpf) + wthlvp_av[:,np.newaxis,np.newaxis]

    wthlvpf_moist = wthlpf_moist + 0.608*thl_av*wqtpf_moist
    wthlvpf_dry = wthlpf_dry + 0.608*thl_av*wqtpf_dry

    wthlvpp_moist = mean_mask(wthlvpp, mask_moist)
    wthlvpp_dry = mean_mask(wthlvpp, mask_dry)
    
    # Scale decompose wthlvf contributions FIXME need to make Galilean invariant
    wthlvpf_l, wthlvpf_c, wthlvpf_r = scaleDecomposeFlux(wff , wfp, thlvpf, thlvpp, circ_mask)
    
    wthlvpf_l_moist_time[i,:] = mean_mask(wthlvpf_l, mask_moist)
    wthlvpf_l_dry_time[i,:] = mean_mask(wthlvpf_l, mask_dry)
    wthlvpf_c_moist_time[i,:] = mean_mask(wthlvpf_c, mask_moist)
    wthlvpf_c_dry_time[i,:] = mean_mask(wthlvpf_c, mask_dry)
    wthlvpf_r_moist_time[i,:] = mean_mask(wthlvpf_r, mask_moist)
    wthlvpf_r_dry_time[i,:] = mean_mask(wthlvpf_r, mask_dry)
    
    wthlpf_moist_time[i,:] = wthlpf_moist
    wthlpf_dry_time[i,:] = wthlpf_dry
    
    wqtpf_moist_time[i,:] = wqtpf_moist
    wqtpf_dry_time[i,:] = wqtpf_dry
    
    wqlpf_moist_time[i,:] = wqlpf_moist
    wqlpf_dry_time[i,:] = wqlpf_dry
    
    wthlvp_av_time[i,:] = wthlvp_av
    wthlvpf_moist_time[i,:] = wthlvpf_moist
    wthlvpf_dry_time[i,:] = wthlvpf_dry
    wthlvpp_moist_time[i,:] = wthlvpp_moist
    wthlvpp_dry_time[i,:] = wthlvpp_dry
    
    
    del wthlpf
    del wqtpf
    del wqlpf
    del wthlvp
    del wthlvpp
    gc.collect()
        
    ## BUDGET TERMS
    
    # Gradient production
    
    # Mean gradients
    Gamma_thl = (thl_av[1:] - thl_av[:-1])/dzh
    Gamma_qt = (qt_av[1:] - qt_av[:-1])/dzh
    # Gamma_ql = (ql_av[1:] - ql_av[:-1])/dzh
    Gamma_thlv = (thlv_av[1:] - thlv_av[:-1])/dzh
    
    Gamma_thl_f = (Gamma_thl[1:] + Gamma_thl[:-1])*0.5
    Gamma_qt_f = (Gamma_qt[1:] + Gamma_qt[:-1])*0.5
    # Gamma_ql_f = (Gamma_ql[1:] + Gamma_ql[:-1])*0.5
    Gamma_thlv_f = (Gamma_thlv[1:] + Gamma_thlv[:-1])*0.5
    
    # Large-scale thlv production
    thlvpf_prod_moist = wf_moist[1:-1]*Gamma_thlv_f
    thlvpf_prod_dry = wf_dry[1:-1]*Gamma_thlv_f
    
    thlvpf_prod_moist_time[i,:] = thlvpf_prod_moist
    thlvpf_prod_dry_time[i,:] = thlvpf_prod_dry
    
    # Small-scale
    thlvpp_prod_moist = wp_moist[1:-1]*Gamma_thlv_f
    thlvpp_prod_dry = wp_dry[1:-1]*Gamma_thlv_f
    
    thlvpp_prod_moist_time[i,:] = thlvpp_prod_moist
    thlvpp_prod_dry_time[i,:] = thlvpp_prod_dry
    
    # Moisture production term with actual w'
    qtpf_prod_wex_moist = wf_moist[1:-1]*Gamma_qt_f
    qtpf_prod_wex_dry = wf_dry[1:-1]*Gamma_qt_f

    qtpf_prod_moist_wex_time[i,:] = qtpf_prod_wex_moist
    qtpf_prod_dry_wex_time[i,:] = qtpf_prod_wex_dry
    
    # wthlv
    wthlvpf_prod = lowPass(wfp**2,circ_mask)[1:-1]*Gamma_thlv_f[:,np.newaxis,np.newaxis]
    wthlvp_prod_av = np.mean(wthlvpf_prod,axis=(1,2))
    wthlvpf_prod_moist_time[i,:] = mean_mask(wthlvpf_prod, mask_moist) - wthlvp_prod_av
    wthlvpf_prod_dry_time[i,:] = mean_mask(wthlvpf_prod, mask_dry) - wthlvp_prod_av
    
    del wthlvpf_prod
    gc.collect()
    
    # Reynolds vertical flux divergence anomaly (with second order scheme)
    
    # thlv
    div_wthlv_r = ddzwx_2nd(whp, thlpp+0.608*thl_av[:,np.newaxis,np.newaxis]*qtpp,
                            dzh, rhobf=rhobfi)
    div_wthlv_av = np.mean(ddzwx_2nd(whf+whp, 
                                     thlpf+thlpp+0.608*thl_av[:,np.newaxis,np.newaxis]*(qtpf+qtpp),
                                     dzh, rhobf=rhobfi),axis=(1,2))
    div_wthlv_rf = lowPass(div_wthlv_r, circ_mask)
    div_wthlv_rp = div_wthlv_r - div_wthlv_rf # Since div_wthlv_rf still includes the mean flux, this is already the anomalous flux

    # Moist/dry and large/small scale
    div_wthlv_rf_moist = mean_mask(div_wthlv_rf, mask_moist)
    div_wthlv_rf_dry = mean_mask(div_wthlv_rf, mask_dry)
    div_wthlv_rp_moist = mean_mask(div_wthlv_rp, mask_moist)
    div_wthlv_rp_dry = mean_mask(div_wthlv_rp, mask_dry)
    
    thlvpf_vdiv_moist_time[i,:] = div_wthlv_rf_moist - div_wthlv_av
    thlvpf_vdiv_dry_time[i,:] = div_wthlv_rf_dry - div_wthlv_av
    thlvpp_vdiv_moist_time[i,:] = div_wthlv_rp_moist
    thlvpp_vdiv_dry_time[i,:] = div_wthlv_rp_dry
    
    # wthlv
    wdiv_wthlvf_r = lowPass(wfp[1:-1,:,:]*div_wthlv_r,circ_mask)
    # wdiv_wthlv_av = lowPass(wfp[1:-1,:,:]*div_wthlv_av[:,np.newaxis,np.newaxis],circ_mask) # <-- basically zero
    wdiv_wthlv_av = np.mean(wfp[1:-1,:,:]*div_wthlv_r,axis=(1,2))

    div_ww_r = ddzww_2nd(whp, dzh, rhobf=rhobfi, rhobh=rhobhi) # At half levels
    div_ww_r = (div_ww_r[1:,:,:] + div_ww_r[:-1,:,:])*0.5 # At full levels (lose highest one)
    
    thlvpdiv_wwf_r = lowPass(thlvpp[1:-2,:,:]*div_ww_r,circ_mask)
    thlvpdiv_ww_av = np.mean(thlvpp[1:-2,:,:]*div_ww_r,axis=(1,2))
    
    wthlvpf_vdiv_moist_time[i,:] = (mean_mask(wdiv_wthlvf_r[:-1,:,:]+thlvpdiv_wwf_r, mask_moist) - 
                               (wdiv_wthlv_av[:-1]+thlvpdiv_ww_av))
    wthlvpf_vdiv_dry_time[i,:] = (mean_mask(wdiv_wthlvf_r[:-1,:,:]+thlvpdiv_wwf_r, mask_dry) - 
                               (wdiv_wthlv_av[:-1]+thlvpdiv_ww_av))
    
    # qt
    div_wqt_r = ddzwx_2nd(whp, qtpp, dzh, rhobf=rhobfi)
    div_wqt_av = np.mean(ddzwx_2nd(whf+whp, qtpf+qtpp, dzh, rhobf=rhobfi),axis=(1,2))
    div_wqt_rf = lowPass(div_wqt_r, circ_mask)
    
    div_wqt_rf_moist = mean_mask(div_wqt_rf,mask_moist)
    div_wqt_rf_dry = mean_mask(div_wqt_rf,mask_dry)
    
    qtpf_vdiv_moist_time[i,:] = div_wqt_rf_moist - div_wqt_av
    qtpf_vdiv_dry_time[i,:] = div_wqt_rf_dry - div_wqt_av

    del div_wthlv_rf
    del div_wthlv_rp    
    del wdiv_wthlvf_r
    del div_ww_r
    del thlvpdiv_wwf_r
    del div_wqt_r
    del div_wqt_rf
    gc.collect()

    # Moisture instability term model (WTG for thlv and qtpf model for flux anomaly div)
    w_mod = lowPass(np.abs(whp),circ_mask)
    div_wthlvfa_mod = ddzwx_2nd(w_mod, -0.608*thl_av[:,np.newaxis,np.newaxis]*qtpf, dzh, rhobf=rhobfi)
    div_wthlvfa_mod_moist = mean_mask(div_wthlvfa_mod, mask_moist)
    div_wthlvfa_mod_dry = mean_mask(div_wthlvfa_mod, mask_dry)
    del w_mod
    del div_wthlvfa_mod
    gc.collect()

    qtpf_prod_moist = div_wthlvfa_mod_moist*Gamma_qt_f/Gamma_thlv_f
    qtpf_prod_dry = div_wthlvfa_mod_dry*Gamma_qt_f/Gamma_thlv_f
    
    # Model that just relies on the WTG
    # qtpf_prod_moist = (div_wthlv_r_moist - div_wthlv_av)*Gamma_qt_f/Gamma_thlv_f
    # qtpf_prod_dry = (div_wthlv_r_dry - div_wthlv_av)*Gamma_qt_f/Gamma_thlv_f
    
    qtpf_prod_moist_time[i,:] = qtpf_prod_moist
    qtpf_prod_dry_time[i,:] = qtpf_prod_dry

    # Horizontal advection

    # Horizontal thlv advection
    div_uhthlvp = ddxhuha_2nd(u, v, thlpf+thlpp+0.608*thl_av[:,np.newaxis,np.newaxis]*(qtpf+qtpp), dx, dy)
    div_uhthlvpf = lowPass(div_uhthlvp, circ_mask)
    
    # moist/dry and large/small scale
    div_uhthlvpf_moist = mean_mask(div_uhthlvpf, mask_moist)
    div_uhthlvpf_dry = mean_mask(div_uhthlvpf, mask_dry)
    div_uhthlvpp_moist = mean_mask(div_uhthlvp - div_uhthlvpf, mask_moist)
    div_uhthlvpp_dry = mean_mask(div_uhthlvp - div_uhthlvpf, mask_dry)

    thlvpf_hdiv_moist_time[i,:] = div_uhthlvpf_moist[1:-1]
    thlvpf_hdiv_dry_time[i,:] = div_uhthlvpf_dry[1:-1]
    thlvpp_hdiv_moist_time[i,:] = div_uhthlvpp_moist[1:-1]
    thlvpp_hdiv_dry_time[i,:] = div_uhthlvpp_dry[1:-1]
    
    # TODO Add wthlvpf anomaly here
    
    # Horizontal moisture advection
    # intra-scale contribution largest, but entire term kept for now
    div_uhqtp = lowPass(ddxhuha_2nd(u, v, qtpf+qtpp, dx, dy), circ_mask)
    div_uhqtp_moist = mean_mask(div_uhqtp,mask_moist)
    div_uhqtp_dry = mean_mask(div_uhqtp,mask_dry)

    qtpf_hdiv_moist_time[i,:] = div_uhqtp_moist[1:-1]
    qtpf_hdiv_dry_time[i,:] = div_uhqtp_dry[1:-1]
    
    del div_uhthlvp
    del div_uhthlvpf
    del div_uhqtp
    del u
    del v
    gc.collect()

    # Subsidence warming
    wsubdthlvpdz = wsubdxdz(wfls[izmin:izmax],thlpf+thlpp+0.608*thl_av[:,np.newaxis,np.newaxis]*(qtpf+qtpp), dzh)
    wsubdthlvpdzf = lowPass(wsubdthlvpdz, circ_mask)
    
    # moist/dry and large/small
    wsubdthlvpdzf_moist = mean_mask(wsubdthlvpdzf, mask_moist)
    wsubdthlvpdzf_dry = mean_mask(wsubdthlvpdzf, mask_dry)
    wsubdthlvpdzp_moist = mean_mask(wsubdthlvpdz - wsubdthlvpdzf, mask_moist)
    wsubdthlvpdzp_dry = mean_mask(wsubdthlvpdz - wsubdthlvpdzf, mask_dry)
    
    thlvpf_subs_moist_time[i,:] = wsubdthlvpdzf_moist[1:]
    thlvpf_subs_dry_time[i,:] = wsubdthlvpdzf_dry[1:]
    thlvpp_subs_moist_time[i,:] = wsubdthlvpdzp_moist[1:]
    thlvpp_subs_dry_time[i,:] = wsubdthlvpdzp_dry[1:]
    
    # Subsidence drying
    wsubdqtpdzf = lowPass(wsubdxdz(wfls[izmin:izmax], qtpf+qtpp, dzh),circ_mask)
    wsubdqtpdzf_moist = mean_mask(wsubdqtpdzf,mask_moist)
    wsubdqtpdzf_dry = mean_mask(wsubdqtpdzf,mask_dry)
    
    qtpf_subs_moist_time[i,:] = wsubdqtpdzf_moist[1:]
    qtpf_subs_dry_time[i,:] = wsubdqtpdzf_dry[1:]
    
    del wsubdthlvpdz
    del wsubdthlvpdzf
    del wsubdqtpdzf
    gc.collect()
    
    # Buoyancy effect on wthlvpf
    wthlvp_buoy = lowPass(thlvpp**2*grav/thlv_av[:,np.newaxis,np.newaxis], circ_mask)
    

    # SFS diffusion
    # Heat
    # diff_thlvp = (diffeka(ekhp+ekh_av[:,np.newaxis,np.newaxis], 
    #                       thlpf+thlpp+0.608*thl_av[:,np.newaxis,np.newaxis]*(qtpf+qtpp),
    #                       dx, dy, zf, rhobfi, rhobhi) +
    #               diffzeka(ekhp, thl_av[:,np.newaxis,np.newaxis]*(1+0.608*qt_av[:,np.newaxis,np.newaxis]),
    #                        dzh, rhobfi, rhobhi))
    # diff_thlvpf = lowPass(diff_thlvp, circ_mask)
    
    # # moist/dry and large/small
    # diff_thlvpf_moist = mean_mask(diff_thlvpf, mask_moist)
    # diff_thlvpf_dry = mean_mask(diff_thlvpf, mask_dry)
    # diff_thlvpp_moist = mean_mask(diff_thlvp - diff_thlvpf, mask_moist)
    # diff_thlvpp_dry = mean_mask(diff_thlvp - diff_thlvpf, mask_dry)

    # thlvpf_diff_moist_time[i,:] = diff_thlvpf_moist
    # thlvpf_diff_dry_time[i,:] = diff_thlvpf_dry
    # thlvpp_diff_moist_time[i,:] = diff_thlvpp_moist
    # thlvpp_diff_dry_time[i,:] = diff_thlvpp_dry
    
    # # Moisture
    # diff_qtpf = lowPass(diffeka(ekhp+ekh_av[:,np.newaxis,np.newaxis], qtpf+qtpp, dx, dy, zf, rhobfi, rhobhi)+
    #                     diffzeka(ekhp, qt_av[:,np.newaxis,np.newaxis], dzh, rhobfi, rhobhi),
    #                     circ_mask)
    # diff_qtpf_moist = mean_mask(diff_qtpf,mask_moist)
    # diff_qtpf_dry = mean_mask(diff_qtpf,mask_dry)

    # qtpf_diff_moist_time[i,:] = diff_qtpf_moist
    # qtpf_diff_dry_time[i,:] = diff_qtpf_dry
if store:
    np.save(lp+'/time.npy',time[plttime])
    np.save(lp+'/plttime.npy',plttime)
    np.save(lp+'/zf.npy',zflim)
    
    np.save(lp+'/qtpf_moist_time.npy',qtpf_moist_time)
    np.save(lp+'/qtpf_dry_time.npy',qtpf_dry_time)
    np.save(lp+'/qtpf_prod_moist_time.npy',qtpf_prod_moist_time)
    np.save(lp+'/qtpf_prod_dry_time.npy',qtpf_prod_dry_time)
    np.save(lp+'/qtpf_prod_moist_wex_time.npy',qtpf_prod_moist_wex_time)
    np.save(lp+'/qtpf_prod_dry_wex_time.npy',qtpf_prod_dry_wex_time)
    np.save(lp+'/qtpf_vdiv_moist_time.npy',qtpf_vdiv_moist_time)
    np.save(lp+'/qtpf_vdiv_dry_time.npy',qtpf_vdiv_dry_time)
    np.save(lp+'/qtpf_hdiv_moist_time.npy',qtpf_hdiv_moist_time)
    np.save(lp+'/qtpf_hdiv_dry_time.npy',qtpf_hdiv_dry_time)
    np.save(lp+'/qtpf_subs_moist_time.npy',qtpf_subs_moist_time)
    np.save(lp+'/qtpf_subs_dry_time.npy',qtpf_subs_dry_time)
    np.save(lp+'/qtpf_diff_moist_time.npy',qtpf_diff_moist_time)
    np.save(lp+'/qtpf_diff_dry_time.npy',qtpf_diff_dry_time)
    
    np.save(lp+'/thlvpf_moist_time.npy',thlvpf_moist_time)
    np.save(lp+'/thlvpf_dry_time.npy',thlvpf_dry_time)
    np.save(lp+'/thlvpf_prod_moist_time.npy',thlvpf_prod_moist_time)
    np.save(lp+'/thlvpf_prod_dry_time.npy',thlvpf_prod_dry_time)
    np.save(lp+'/thlvpf_vdiv_moist_time.npy',thlvpf_vdiv_moist_time)
    np.save(lp+'/thlvpf_vdiv_dry_time.npy',thlvpf_vdiv_dry_time)
    np.save(lp+'/thlvpf_hdiv_moist_time.npy',thlvpf_hdiv_moist_time)
    np.save(lp+'/thlvpf_hdiv_dry_time.npy',thlvpf_hdiv_dry_time)
    np.save(lp+'/thlvpf_subs_moist_time.npy',thlvpf_subs_moist_time)
    np.save(lp+'/thlvpf_subs_dry_time.npy',thlvpf_subs_dry_time)
    np.save(lp+'/thlvpf_diff_moist_time.npy',thlvpf_diff_moist_time)
    np.save(lp+'/thlvpf_diff_dry_time.npy',thlvpf_diff_dry_time)
    
    np.save(lp+'/thlvpp_moist_time.npy',thlvpp_moist_time)
    np.save(lp+'/thlvpp_dry_time.npy',thlvpp_dry_time)
    np.save(lp+'/thlvpp_prod_moist_time.npy',thlvpp_prod_moist_time)
    np.save(lp+'/thlvpp_prod_dry_time.npy',thlvpp_prod_dry_time)
    np.save(lp+'/thlvpp_vdiv_moist_time.npy',thlvpp_vdiv_moist_time)
    np.save(lp+'/thlvpp_vdiv_dry_time.npy',thlvpp_vdiv_dry_time)
    np.save(lp+'/thlvpp_hdiv_moist_time.npy',thlvpp_hdiv_moist_time)
    np.save(lp+'/thlvpp_hdiv_dry_time.npy',thlvpp_hdiv_dry_time)
    np.save(lp+'/thlvpp_subs_moist_time.npy',thlvpp_subs_moist_time)
    np.save(lp+'/thlvpp_subs_dry_time.npy',thlvpp_subs_dry_time)
    np.save(lp+'/thlvpp_diff_moist_time.npy',thlvpp_diff_moist_time)
    np.save(lp+'/thlvpp_diff_dry_time.npy',thlvpp_diff_dry_time)
    
    np.save(lp+'/thl_av_time.npy',thl_av_time)
    np.save(lp+'/thlv_av_time.npy',thlv_av_time)
    np.save(lp+'/qt_av_time.npy',qt_av_time)
    
    np.save(lp+'/thlpf_moist_time.npy',thlpf_moist_time)
    np.save(lp+'/thlpf_dry_time.npy',thlpf_dry_time)
    np.save(lp+'/wff_moist_time.npy',wff_moist_time)
    np.save(lp+'/wff_dry_time.npy',wff_dry_time)
    np.save(lp+'/qlpf_moist_time.npy',qlpf_moist_time) 
    np.save(lp+'/qlpf_dry_time.npy',qlpf_dry_time)
    
    np.save(lp+'/thlpp_moist_time.npy',thlpp_moist_time)
    np.save(lp+'/thlpp_dry_time.npy',thlpp_dry_time)
    np.save(lp+'/wfp_moist_time.npy',wfp_moist_time)
    np.save(lp+'/wfp_dry_time.npy',wfp_dry_time)
    np.save(lp+'/qlpp_moist_time.npy',qlpp_moist_time) 
    np.save(lp+'/qlpp_dry_time.npy',qlpp_dry_time)
    
    np.save(lp+'/wthlpf_moist_time.npy',wthlpf_moist_time)
    np.save(lp+'/wthlpf_dry_time.npy',wthlpf_dry_time)
    
    np.save(lp+'/wqtpf_moist_time.npy',wqtpf_moist_time)
    np.save(lp+'/wqtpf_dry_time.npy',wqtpf_dry_time)
    
    np.save(lp+'/wqlpf_moist_time.npy',wqlpf_moist_time)
    np.save(lp+'/wqlpf_dry_time.npy',wqlpf_dry_time)
    
    np.save(lp+'/wthlvp_av_time.npy',wthlvp_av_time)
    np.save(lp+'/wthlvpf_moist_time.npy',wthlvpf_moist_time)
    np.save(lp+'/wthlvpf_dry_time.npy',wthlvpf_dry_time)
    np.save(lp+'/wthlvpf_l_moist_time.npy',wthlvpf_l_moist_time)
    np.save(lp+'/wthlvpf_l_dry_time.npy',wthlvpf_l_dry_time)
    np.save(lp+'/wthlvpf_c_moist_time.npy',wthlvpf_c_moist_time)
    np.save(lp+'/wthlvpf_c_dry_time.npy',wthlvpf_c_dry_time)
    np.save(lp+'/wthlvpf_r_moist_time.npy',wthlvpf_r_moist_time)
    np.save(lp+'/wthlvpf_r_dry_time.npy',wthlvpf_r_dry_time)
    np.save(lp+'/wthlvpp_moist_time.npy',wthlvpp_moist_time)
    np.save(lp+'/wthlvpp_dry_time.npy',wthlvpp_dry_time)
    

