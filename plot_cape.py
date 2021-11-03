#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:56:45 2021

@author: janssens
"""


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import netCDF4 as nc
import seaborn as sns
from functions import getRad, lowPass

lp = '/scratch-shared/janssens/bomex200aswitch/a2'
sp = lp+'/figs'

tPlot = 24
klp = 4

ds = nc.Dataset(lp+'/cape2d.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)

it = np.argmin(abs(tPlot-time))
extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000

# Albedo
qli = np.ma.getdata(ds.variables['lwp'][it,:,:])
tau = 0.19*qli**(5./6)*7e7**(1/3)
alb = tau/(6.8+tau)
alb[qli<1e-7] = 0
alb[qli>=1e-7] = 1

# TWP fluctuation
twpp = np.ma.getdata(ds.variables['twp'][it,:,:])
twpp -= np.mean(twpp)

# Low-pass filtered TWP fluctuation
circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1
twppf = lowPass(twpp, circ_mask)


fig,axs = plt.subplots(ncols=3,figsize=(9,3),sharey=True)
sc = axs[0].imshow(twpp,extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
axs[0].set_xlabel('x [km]')
axs[0].set_ylabel('y [km]')

axs[1].imshow(twppf,extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
axs[1].contour(twppf,levels=[0],extent=extent,origin='upper')
axs[1].set_xlabel('x [km]')

axs[2].contour(twppf,levels=[0],extent=extent,origin='upper',colors='white')
axs[2].imshow(alb,extent=extent,vmin=0,vmax=1,cmap='gist_gray')
axs[2].set_xlabel('x [km]')

# cbax = fig.add_axes([1, 0.1, 0.02, 0.85])
cbax = fig.add_axes([-0.06, 0.1, 0.02, 0.85])
cb = fig.colorbar(sc, cax=cbax)
cb.ax.set_ylabel(r"Total Water Path fluctuation [kg/kg/m$^2$]", rotation=90, labelpad=-65)
plt.tight_layout()
plt.savefig(sp+'/twpfluct.pdf',bbox_inches='tight')


