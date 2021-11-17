#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:56:45 2021

@author: janssens
"""


import numpy as np
import pandas as pd
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import seaborn as sns
from functions import getRad, lowPass

lp = '/scratch-shared/janssens/bomex200aswitch/a2'
sp = lp+'/figs'

klp = 4
qlc = 1e-7

ds = nc.Dataset(lp+'/cape2d.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)

extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000

circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

#%% Plot twp, twppf and clouds at a single time step

tPlot = 24

it = np.argmin(abs(tPlot-time))

# Albedo
qli = np.ma.getdata(ds.variables['lwp'][it,:,:])
# tau = 0.19*qli**(5./6)*7e7**(1/3)
# alb = tau/(6.8+tau)
alb = qli.copy()
alb[qli<qlc] = 0
alb[qli>=qlc] = 1

# TWP fluctuation
twpp = np.ma.getdata(ds.variables['twp'][it,:,:])
twpp -= np.mean(twpp)

# Low-pass filtered TWP fluctuation
twppf = lowPass(twpp, circ_mask)

fig,axs = plt.subplots(ncols=2,figsize=(8,4),sharey=True)
sc = axs[0].imshow(twpp,extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
axs[0].set_xlabel('x [km]')
axs[0].set_ylabel('y [km]')

axs[1].imshow(twppf,extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
axs[1].contour(twppf,levels=[0],extent=extent,origin='upper')
axs[1].set_xlabel('x [km]')

# axs[2].contour(twppf,levels=[0],extent=extent,origin='upper',colors='white')
# axs[2].imshow(alb,extent=extent,vmin=0,vmax=1,cmap='gist_gray')
# axs[2].set_xlabel('x [km]')

pos1 = axs[1].get_position()
cbax = fig.add_axes([.92, pos1.ymin, 0.01, pos1.height])
# cbax = fig.add_axes([1, 0.1, 0.02, 0.85])
# cbax = fig.add_axes([-0.06, 0.1, 0.02, 0.85])
cb = fig.colorbar(sc, cax=cbax)
cb.ax.set_ylabel(r"Total Water Path fluctuation [kg/kg/m$^2$]", rotation=270, labelpad=15) #-65
# plt.tight_layout()
plt.savefig(sp+'/twpfluct.pdf',bbox_inches='tight',dpi=300)

#%% Plot the time evolution of twpp

tPlot = np.arange(6,18,3)

fig,axs = plt.subplots(ncols=len(tPlot),nrows=2,figsize=(3.5*len(tPlot),7),
                       sharex=True,sharey=True,squeeze=False)

for j in range(len(tPlot)):
    it = np.argmin(abs(tPlot[j]-time))
    twpp = np.ma.getdata(ds.variables['twp'][it,:,:])
    twpp -= np.mean(twpp)
    twppf = lowPass(twpp, circ_mask)
    
    cm = np.ma.getdata(ds.variables['cldtop'][it,:,:])
    # cm = np.zeros(qli.shape)
    # cm[qli<qlc] = 0
    # cm[qli>=qlc] = 1

    sc1 = axs[0,j].imshow(twpp, extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
    
    sc2 = axs[1,j].imshow(cm  , extent=extent,vmin=0 ,vmax=2000,cmap='gray')
    
    if j > 1:
        axs[0,j].contour(twppf,levels=[0],extent=extent,origin='upper',
                         linewidths=0.1,colors='black')
        axs[1,j].contour(twppf,levels=[0],extent=extent,origin='upper',
                         linewidths=0.1,colors='white')
    
    axs[1,j].set_xlabel('x [km]')
    axs[0,j].set_title('t = %.1f hr'%tPlot[j])
    if j == 0:
        axs[0,j].set_ylabel('y [km]')
        axs[1,j].set_ylabel('y [km]')

    if j == len(tPlot)-1:
        pos1 = axs[0,j].get_position()
        cbax1 = fig.add_axes([.92, pos1.ymin, 0.01, pos1.height])
        cb1 = fig.colorbar(sc1, cax=cbax1)
        cb1.ax.set_ylabel(r"$TWP'$ [kg/kg/m$^2$]", rotation=270, labelpad=15)
        
        pos2 = axs[1,j].get_position()
        cbax2 = fig.add_axes([.92, pos2.ymin, 0.01, pos2.height])
        cb2 = fig.colorbar(sc2, cax=cbax2)
        cb2.ax.set_ylabel(r"Cloud-top height [m]", rotation=270, labelpad=15)

plt.savefig(sp+'/twp_cld_evo.pdf', bbox_inches='tight',dpi=300)
    
    

