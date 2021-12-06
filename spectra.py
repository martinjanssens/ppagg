#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:23:26 2021

@author: janssens
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import gc
from functions import *
import argparse

parseFlag = True

if parseFlag:
    parser = argparse.ArgumentParser(description="Compute spectra of selected variables from DALES 3D fields")
    parser.add_argument("--dir", metavar="DIR", type=str, default=".", help="Directory to load/store data from/to")
    parser.add_argument("--itmin", metavar="N", type=int, default=0, help="First time index")
    parser.add_argument("--itmax", metavar="N", type=int, default=-1, help="Last time index")
    parser.add_argument("--dt", metavar="N", type=int, default=1, help="Time sampling interval")
    parser.add_argument("--izmin", metavar="N", type=int, default=0, help="First height index")
    parser.add_argument("--izmax", metavar="N", type=int, default=80, help="Last height index")
    parser.add_argument("--dz", metavar="N", type=int, default=1, help="Height isampling interval")
    parser.add_argument("--klp", metavar="N", type=int, default=4, help="Cutoff wavenumber for lw-pass filter")
    parser.add_argument("--store", action="store_true", default=False, help="Saves the output if given")

    args = parser.parse_args()

    lp = args.dir
    itmin = args.itmin
    itmax = args.itmax
    dti = args.dt
    izmin = args.izmin
    izmax = args.izmax
    dzi = args.dz
    klp = args.klp
    store = args.store
else:
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

#%% Compute spectra at a given height
if not parseFlag:
    # Run specifics
    itmin = 95#23
    itmax = 96
    dti   = 1
    izmin = 0
    izmax = 80
    dzi   = 4
    klp = 4
    store=True

plttime = np.arange(itmin, itmax, dti)
zflim = zf[izmin:izmax:dzi]
N = xf.size; N2 = N//2

spec_qt = np.zeros((len(plttime),len(zflim),N2))
spec_thl = np.zeros((len(plttime),len(zflim),N2))
spec_thlv = np.zeros((len(plttime),len(zflim),N2))
spec_w = np.zeros((len(plttime),len(zflim),N2))
spec_ql = np.zeros((len(plttime),len(zflim),N2))
spec_wqt = np.zeros((len(plttime),len(zflim),N2))
spec_wthl = np.zeros((len(plttime),len(zflim),N2))
spec_wthlv = np.zeros((len(plttime),len(zflim),N2))
spec_wql = np.zeros((len(plttime),len(zflim),N2))

c = 0
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

    thl = thl - thl_av[:,np.newaxis,np.newaxis]
    thlv = thlv - np.mean(thlv,axis=(1,2))[:,np.newaxis,np.newaxis]
    qt = qt - np.mean(qt,axis=(1,2))[:,np.newaxis,np.newaxis]
    ql = ql - np.mean(qt,axis=(1,2))[:,np.newaxis,np.newaxis]
    
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    del wh
    
    for iz in range(len(zflim)):
        if c%10 == 0:
            print('Computing spectra at time step', i+1, '/', len(plttime),
                  ', height', iz+1,'/',len(zflim), 'total', c+1,'/',len(zflim)*len(plttime))
        k1d,spec_qt[i,iz,:] = compute_spectrum(qt[iz,:,:], dx)
        k1d,spec_thl[i,iz,:] = compute_spectrum(thl[iz,:,:], dx)
        k1d,spec_thlv[i,iz,:] = compute_spectrum(thlv[iz,:,:], dx)
        k1d,spec_w[i,iz,:] = compute_spectrum(wf[iz,:,:], dx)
        k1d,spec_ql[i,iz,:] = compute_spectrum(ql[iz,:,:], dx)
        
        k1d,spec_wqt[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, qt[iz,:,:])
        k1d,spec_wthl[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, thl[iz,:,:])
        k1d,spec_wthlv[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, thlv[iz,:,:])
        k1d,spec_wql[i,iz,:] = compute_spectrum(wf[iz,:,:], dx, ql[iz,:,:])
        
        # These are not co-spectra, even if they are probably more accurate
        # k1d,spec_wqt[i,iz,:] = compute_spectrum(wf[iz,:,:]*qt[iz,:,:], dx,sqrt=True)
        # k1d,spec_wthl[i,iz,:] = compute_spectrum(wf[iz,:,:]*thl[iz,:,:], dx,sqrt=True)
        # k1d,spec_wthlv[i,iz,:] = compute_spectrum(wf[iz,:,:]*thlv[iz,:,:], dx,sqrt=True)
        
        c += 1
        
    gc.collect()
if store:
    np.save(lp+'/time_spec.npy',time)
    np.save(lp+'/plttime_spec.npy',plttime)
    np.save(lp+'/zf_spec.npy',zflim)
    np.save(lp+'/k1d.npy',k1d)
    
    np.save(lp+'/spec_qt.npy',spec_qt)
    np.save(lp+'/spec_thl.npy',spec_thl)
    np.save(lp+'/spec_thlv.npy',spec_thlv)
    np.save(lp+'/spec_w.npy',spec_w)
    np.save(lp+'/spec_ql.npy',spec_ql)
    np.save(lp+'/spec_wqt.npy',spec_wqt)
    np.save(lp+'/spec_wthl.npy',spec_wthl)
    np.save(lp+'/spec_wthlv.npy',spec_wthlv)
    np.save(lp+'/spec_wql.npy',spec_wql)
    
