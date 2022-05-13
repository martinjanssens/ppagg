#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:56:45 2021

@author: janssens
"""


import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib.animation as animation
from functions import getRad, lowPass, mean_mask

lp = '/scratch-shared/janssens/eurec4a_old/eurec4a_mean_ssthet'
sp = lp+'/figs'

klp = 4
qlc = 1e-7

ds = nc.Dataset(lp+'/fielddump.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell edges (h in mhh)
zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)

extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000

circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

#%% Make movie over time
tpltmin = 0.0
tpltmax = 24.0

var = 'thl'
zplt = 20.

itpltmin = np.where(time>=tpltmin)[0][0]
itpltmax = np.where(time<tpltmax)[0][-1]+1
iz = np.argmin(np.abs(zplt-zf))

tPlot = time[itpltmin:itpltmax]

def animate(i):
    ax.collections = []
    twpp = np.ma.getdata(ds.variables[var][itpltmin+i,iz,:,:])
    twpp -= np.mean(twpp)
    # twpp = lowPass(twpp, circ_mask)
    
    sc1 = ax.imshow(twpp, extent=extent,vmin=-0.3,vmax=0.3,cmap='RdYlBu_r')
    # ax.contour(twppf,levels=[0],extent=extent,origin='upper',linewidths=1,colors='black')
    ax.set_xlabel('x [km]')
    ax.set_xlabel('y [km]')
    ax.set_title('t = %.1f hr'%tPlot[i])

    pos1 = ax.get_position()
    cbax1 = fig.add_axes([.92, pos1.ymin, 0.01, pos1.height])
    cb1 = fig.colorbar(sc1, cax=cbax1)
    cb1.ax.set_ylabel(r"Lowest model level $\theta_l'$ [K]", rotation=270, labelpad=15)

fig, ax = plt.subplots(figsize=(6, 5))

ani = animation.FuncAnimation(fig, animate, interval=100, frames=len(tPlot))

ani.save(sp+"/twppf.mp4")
