#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:56:45 2021

@author: janssens
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib.animation as animation
from functions import getRad, lowPass, mean_mask

lp = '/scratch-shared/janssens/eurec4a_mean_200km'
sp = lp+'/figs'

klp = 4
qlc = 1e-7
izmin = 0
izmax = 80

ds = nc.Dataset(lp+'/cape2d.001.nc')
ds3 = nc.Dataset(lp+'/fielddump.001.nc')
ds1 = nc.Dataset(lp+'/profiles.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)
#zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)

time3 = np.ma.getdata(ds3.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds3.variables['zt'][:]) # Cell edges (h in mhh)
zflim = zf[izmin:izmax]
dzh = np.zeros(zf.shape)
dzh[1:] = zf[1:] - zf[:-1] # First value is difference mid 1st cell and mid 1st cell below ground
dzh[0] = 2*zf[1]
dzhlim = dzh[izmin:izmax]

time1 = np.ma.getdata(ds1.variables['time'][:]) / 3600

extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000

circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

#%% Make movie over time
tpltmin = 0.0
tpltmax = 24.0

# var = 'thl'
# zplt = 20.

itpltmin = np.where(time>=tpltmin)[0][0]
itpltmax = np.where(time<tpltmax)[0][-1]+1
# iz = np.argmin(np.abs(zplt-zf))

tPlot = time[itpltmin:itpltmax]

it3p = 0
def animate(i):
    global it3p
    global div_sc_ci
    if i == 0:
        it3p=0
    ti = tPlot[i]
    it3d = np.argmin(np.abs(ti-time3))
    print('time 2d:',ti)
    print('time 3d:',time3[it3d])
    print('time 3d prev:',time3[it3p])
    if it3d > it3p or it3d == 0:
        print('Making 3D plot')
        it1d = np.argmin(np.abs(time3[it3d] - time1))
        
        # Determine mesoscale divergence
        wf = np.ma.getdata(ds3.variables['w'][it3d,izmin:izmax+1,:,:])
        wf = (wf[1:,:,:] + wf[:-1,:,:])*0.5
        wff = lowPass(wf, circ_mask)
        div_w = -(wff[1:,:,:] - wff[:-1,:,:])/dzhlim[1:,np.newaxis,np.newaxis]
        div_w = (div_w[1:,:,:] + div_w[:-1,:,:])*0.5

        # Determine cloud base and inversion base
        ql_av = np.ma.getdata(ds1.variables['ql'][it1d,izmin:izmax])
        try:
            zcb = zflim[ql_av>1e-8][0]
        except:
            zcb = 500
        izcb = np.argmin(np.abs(zcb - zflim[1:-1]))
        
        thl_av = np.ma.getdata(ds1.variables['thl'][it1d,izmin:izmax])
        izinv = np.argmax((thl_av[1:] - thl_av[:-1])/dzhlim[1:])
                
        # Determine subcloud-layer/cloud-layer divergence
        div_sc = np.mean(div_w[:izcb,:,:],axis=0)
        div_cl = np.mean(div_w[izcb:izinv,:,:],axis=0)
        div_cl_sc = div_sc/div_cl

        # Subcloud layer divergence in columns with converging/diverging structure
        div_sc_ci = div_sc
        div_sc_ci[div_cl_sc > 0] = np.nan

    it3p = it3d

    #axs.collections = []
    twpp = np.ma.getdata(ds.variables['twp'][itpltmin+i,:,:])
    twpp -= np.mean(twpp)
    # twpp = lowPass(twpp, circ_mask)
    
    cldtop = np.ma.getdata(ds.variables['cldtop'][itpltmin+i,:,:])

    sc1 = axs[0].imshow(twpp, extent=extent,vmin=-1.25,vmax=1.25,cmap='RdYlBu')
    # ax.contour(twppf,levels=[0],extent=extent,origin='upper',linewidths=1,colors='black')
    axs[0].set_xlabel('x [km]')
    axs[0].set_ylabel('y [km]')
    # axs[0].set_title('t = %.1f hr'%tPlot[i])

    pos1 = axs[0].get_position()
    cbax1 = fig.add_axes([.355, pos1.ymin, 0.01, pos1.height])
    cb1 = fig.colorbar(sc1, cax=cbax1)
    cb1.ax.set_ylabel(r"Total water path anomaly", rotation=270, labelpad=15)
    
    sc2 = axs[1].imshow(cldtop, extent=extent,vmin=0,vmax=2000,cmap='Greys_r')
    # ax.contour(twppf,levels=[0],extent=extent,origin='upper',linewidths=1,colors='black')
    axs[1].set_xlabel('x [km]')
    axs[1].set_title('t = %.1f hr'%tPlot[i])

    pos1 = axs[1].get_position()
    cbax2 = fig.add_axes([.625, pos1.ymin, 0.01, pos1.height])
    cb2 = fig.colorbar(sc2, cax=cbax2)
    cb2.ax.set_ylabel(r"Cloud-top height [m]", rotation=270, labelpad=15)

    cmap = mpl.cm.get_cmap("RdBu_r").copy()
    cmap.set_bad('white')
    sc3 = axs[2].imshow(div_sc_ci,cmap=cmap,vmin=-5e-5,vmax=5e-5,extent=extent)

#ax.contour(twppf,levels=[0],extent=extent,origin='upper',
#                 linewidths=1,colors='black')

    pos3 = axs[2].get_position()
    cbax3 = fig.add_axes([.9, pos3.ymin, 0.01, pos3.height])
    cb1 = fig.colorbar(sc3, cax=cbax3)
    cb1.ax.set_ylabel(r"$\mathcal{D}$ s$^{-1}$", rotation=270, labelpad=15)

    axs[0].set_xlabel('x [km]')

fig, axs = plt.subplots(figsize=(18, 5),ncols=3, sharey=True)

ani = animation.FuncAnimation(fig, animate, interval=100, frames=len(tPlot))

ani.save(sp+"/twp_cldtop.mp4")
