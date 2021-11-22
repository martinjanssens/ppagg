#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 19:28:54 2021

@author: janssens
"""


import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
sys.path.insert(1, '/home/janssens/scripts/pp3d/')
from functions import *
from thermofunctions import qsatur, rlv, rd, rv, cp

lp = '/scratch-shared/janssens/bomex200_e12'
sp = lp+'/figs'
itmin = 57
itmax = 58
di    = 1
izmin = 0
izmax = 80
klp = 4

ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])
rhobh = np.ma.getdata(ds1.variables['rhobh'][:])

dx = np.diff(xf)[0]
dy = np.diff(xf)[0] # Assumes uniform horizontal spacing
dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain

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
    wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    
    presh  = np.ma.getdata(ds1.variables['presh'][it1d,izmin:izmax])
    presf  = (presh[1:]+presh[:-1])*0.5
    exnf   = (presf/1e5)**(Rd/cp)
    
    T = exnf[:,np.newaxis,np.newaxis]*thl[:-1,:,:] + rlv/cp*ql[:-1,:,:]
    
    # Compute qs
    qs = qsatur(T,presf[:,np.newaxis,np.newaxis])
    
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qt[:,:,:],zflim,axis=0)
    twpp=twp-np.mean(twp)
    twppf = lowPass(twpp,circ_mask)
    
    # Sort by twp
    itwpsort = np.argsort(twppf.flatten())
    ql_sort = ql.reshape(zflim.size,-1)[:,itwpsort]
    
    fig = plt.figure(figsize=(6,4)); ax = plt.gca()
    sc = ax.imshow(np.flipud(ql_sort),
                   extent=[twppf.min(),twppf.max(),zflim[0],zflim[-1]],
                   aspect='auto',cmap='YlGnBu',
                   vmin=0,
                   vmax=5e-5)
    ax.set_xlabel(r"$TWP_m'$ [kg/m$^2$]")
    ax.set_ylabel(r'z [m]')
    pos = ax.get_position()
    cbax = fig.add_axes([0.95, pos.ymin, 0.01, pos.height])
    cb1 = fig.colorbar(sc, cax=cbax)
    cb1.ax.set_ylabel(r"$q_l$ [kg/kg]", rotation=270, labelpad=15)
    plt.savefig(sp+'/ql_twppf.pdf',bbox_inches='tight',dpi=300)