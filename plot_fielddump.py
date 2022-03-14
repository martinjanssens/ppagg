#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 19:28:54 2021

@author: janssens
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
import netCDF4 as nc
import sys
sys.path.insert(1, '/home/janssens/scripts/pp3d/')
from functions import *
from thermofunctions import qsatur, rlv, rd, rv, cp

lp = '/scratch-shared/janssens/eurec4a_test'
sp = lp+'/figs'
itmin = 63
itmax = 64
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
    u = np.ma.getdata(ds.variables['u'][plttime[i],izmin:izmax,:,:])
    v = np.ma.getdata(ds.variables['v'][plttime[i],izmin:izmax,:,:])
    wf = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][plttime[i],izmin:izmax,:,:])
    if 'p' in ds.variables.keys():
        p = np.ma.getdata(ds.variables['p'][plttime[i],izmin:izmax,:,:])
    
    # w  h-> f-levels
    wf = (wf[1:,:,:] + wf[:-1,:,:])*0.5
    
    # Tempereature and saturation specific humidity
    presh  = np.ma.getdata(ds1.variables['presh'][it1d,izmin:izmax])
    presf  = (presh[1:]+presh[:-1])*0.5
    exnf   = (presf/1e5)**(Rd/cp)
    
    T = exnf[:,np.newaxis,np.newaxis]*thl[:-1,:,:] + rlv/cp*ql[:-1,:,:]
    qs = qsatur(T,presf[:,np.newaxis,np.newaxis])
    
    # Scale decompose
    twp = np.trapz(rhobfi[:,np.newaxis,np.newaxis]*qt[:,:,:],zflim,axis=0)
    twpp=twp-np.mean(twp)
    twppf = lowPass(twpp,circ_mask)
    
    wff = lowPass(wf, circ_mask)
    wfp = wf - wff
    
    up = u - np.mean(u,axis=(1,2))[:,np.newaxis,np.newaxis]
    upf = lowPass(up, circ_mask)
    
    qt_av = np.mean(qt,axis=(1,2))
    qtp = qt - qt_av[:,np.newaxis,np.newaxis]
    qtpf = lowPass(qtp, circ_mask)
    
    thl_av = np.mean(thl,axis=(1,2))
    
    thlv = thl + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt
    thlv_av = np.mean(thlv,axis=(1,2))
    thlvp = thlv - thlv_av[:,np.newaxis,np.newaxis]
    
    ql_av = np.mean(ql,axis=(1,2))
    qlp = ql - ql_av[:,np.newaxis,np.newaxis]
    
    wql = wf*ql
    wqlf = lowPass(wql,circ_mask)
    wqlp = wf*qlp
    wqlpf = lowPass(wqlp,circ_mask)
    wthlvp = wf*thlvp
    wthlvpf = lowPass(wthlvp,circ_mask)
    
    ddxhuh = (np.roll(u,-1,axis=2) - u)/dx + (np.roll(v,-1,axis=1) - v)/dy
    
    # Sort by twp
    itwpsort = np.argsort(twppf.flatten())
    ql_sort = ql.reshape(zflim.size,-1)[:,itwpsort]
    wf_sort = wf.reshape(zflim.size,-1)[:,itwpsort]
    
    fig,axs = plt.subplots(nrows=2,sharex=True,figsize=(6,8))
    sc0 = axs[0].imshow(np.flipud(ql_sort),
                        extent=[twppf.min(),twppf.max(),zflim[0],zflim[-1]],
                        aspect='auto',cmap='YlGnBu',
                        vmin=0,
                        vmax=5e-5)
    axs[0].set_ylabel(r'z [m]')
    pos0 = axs[0].get_position()
    cbax0 = fig.add_axes([0.95, pos0.ymin, 0.01, pos0.height])
    cb0 = fig.colorbar(sc0, cax=cbax0)
    cb0.ax.set_ylabel(r"$q_l$ [kg/kg]", rotation=270, labelpad=15)
    
    sc1 = axs[1].imshow(np.flipud(wf_sort),
                       extent=[twppf.min(),twppf.max(),zflim[0],zflim[-1]],
                       aspect='auto',cmap='RdYlBu_r',
                       vmin=-0.15,
                       vmax=0.15)
    axs[1].set_xlabel(r"$TWP_m'$ [kg/m$^2$]")
    axs[1].set_ylabel(r'z [m]')
    pos1 = axs[1].get_position()
    cbax1 = fig.add_axes([0.95, pos1.ymin, 0.01, pos1.height])
    cb1 = fig.colorbar(sc1, cax=cbax1)
    cb1.ax.set_ylabel(r"$w$ [m/s]", rotation=270, labelpad=15)
    
    plt.savefig(sp+'/ql_w_twppf.pdf',bbox_inches='tight',dpi=300)
    plt.show()
    
    # Flux origins
    extent=[xf.min()/1000,xf.max()/1000,zflim[0],zflim[-1]-1]
    ix = 380
    iz = 37
    fig,axs = plt.subplots(figsize=(10,4),nrows=2,sharex=True,gridspec_kw={'height_ratios':[1,3]})
    sc0 = axs[1].imshow(np.flipud(qt[:-1,ix,:]-qs[:,ix,:]),
                            extent=extent,
                            aspect='auto',cmap='YlGnBu',
                            vmin=-0.01,vmax=0.003)
    pos = axs[1].get_position()
    cbax0=fig.add_axes([.92, pos.ymin, 0.007, pos.height])
    cb0 = fig.colorbar(sc0, cax=cbax0)
    cb0.ax.set_ylabel(r"$q_t- q_s$ [kg/kg]", rotation=270, labelpad=15)

    sc1 = axs[1].contour(qtpf[:-1,ix,:],extent=extent,origin='lower',linewidths=0.75,cmap='Blues',levels=np.linspace(0.0002,0.003,6))
    cbax1=fig.add_axes([1.02, pos.ymin, 0.007, pos.height])
    norm = colors.Normalize(vmin=sc1.cvalues.min(), vmax=sc1.cvalues.max())
    sm = plt.cm.ScalarMappable(norm=norm, cmap=sc1.cmap)
    sm.set_array([])
    cb1 = fig.colorbar(sm, cax=cbax1)
    cb1.ax.set_ylabel(r"$q_{t_m}'$ [kg/kg]", rotation=270, labelpad=15)
    cb1.locator = ticker.MaxNLocator(nbins=6)
    cb1.update_ticks()
    axs[1].contour(ql[:-1,ix,:],levels=[1e-7],extent=extent,origin='lower',linewidths=0.5,colors='black')
    axs[1].axhline(zflim[iz],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
    axs[1].set_xlabel('x [km]')
    axs[1].set_ylabel('z [m]')

    axs[0].plot(xf/1000,(wthlvp[iz,ix,:]-np.mean(wthlvp[iz,:,:])),c='maroon',label=r"$F_{\theta_{lv}'}$")
    axs[0].plot(xf/1000,(wthlvpf[iz,ix,:]-np.mean(wthlvp[iz,:,:])),c='maroon',linestyle='dashed',label=r"$F_{{\theta_{lv}'}_m}$")
    axs[0].plot(xf/1000,-7*thl_av[iz]*(wqlp[iz,ix,:]-np.mean(wqlp[iz,:,:])),c='seagreen',label=r"$-7\overline{\theta_l}F_{q_l'}$")
    axs[0].plot(xf/1000,-7*thl_av[iz]*(wqlpf[iz,ix,:]-np.mean(wqlp[iz,:,:])),c='seagreen',linestyle='dashed',label=r"$-7\overline{\theta_l}F_{{q_l'}_m}$")
    axs[0].set_ylim((-1,0))
    axs[0].set_ylabel(r"$F_{\theta_{lv}}$, [Km/s]")
    axs[0].legend(bbox_to_anchor=(1,1),loc='upper left',ncol=2)

    plt.savefig(sp+'/structure.pdf',bbox_inches='tight',dpi=300)
    plt.show()
    
    # Conceptual figure
    extent=[xf.min()/1000,xf.max()/1000,zflim[0],zflim[-1]]
    ix = 372
    iz0 = 15
    iz1 = 37
    # fig,axs = plt.subplots(figsize=(10,4),nrows=1)
    fig = plt.figure(figsize=(10,4)); axs = plt.gca()
    sc0 = axs.imshow(np.flipud(qt[:,ix,:]),
                            extent=extent,
                            aspect='auto',cmap='RdYlBu',
                            vmin=0.01,vmax=0.018)
    pos = axs.get_position()
    cbax0=fig.add_axes([.92, pos.ymin, 0.007, pos.height])
    cb0 = fig.colorbar(sc0, cax=cbax0)
    cb0.ax.set_ylabel(r"$q_t$ [kg/kg]", rotation=270, labelpad=15)
    
    # sc1 = axs.contour(qtpf[:-1,ix,:],extent=extent,origin='lower',linewidths=0.75,cmap='Blues',levels=np.linspace(0.0002,0.003,6))
    sc1 = axs.contour(wqlpf[:-1,ix,:],extent=extent,origin='lower',linewidths=1,cmap='Blues',levels=np.linspace(0.0001,0.0002,6))
    cbax1=fig.add_axes([1.02, pos.ymin, 0.007, pos.height])
    norm = colors.Normalize(vmin=sc1.cvalues.min(), vmax=sc1.cvalues.max())
    sm = plt.cm.ScalarMappable(norm=norm, cmap=sc1.cmap)
    sm.set_array([])
    cb1 = fig.colorbar(sm, cax=cbax1)
    cb1.ax.set_ylabel(r"${w'q_l'}_m$ [m/s kg/kg]", rotation=270, labelpad=15)
    cb1.locator = ticker.MaxNLocator(nbins=6)
    cb1.update_ticks()
    axs.contour(ql[:,ix,:],levels=[1e-7],extent=extent,origin='lower',linewidths=1,colors='black')
    axs.contour(qtpf[:,ix,:],levels=[0.0004],extent=extent,origin='lower',linewidths=1,colors='black',linestyles='dashed')
    xfstr = np.linspace(xf[0],xf[-1],len(xf))/1000
    [X,Y] = np.meshgrid(xfstr,zflim)
    speed = np.sqrt((upf[:,ix,:]/1000)**2+wff[:,ix,:]**2)
    lws = np.maximum(.75,3.*speed/np.max(speed))
    st = axs.streamplot(X,Y,upf[:,ix,:]/1000,wff[:,ix,:],linewidth=lws, density=1,color='grey')
    axs.axhline(zflim[iz0],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
    axs.axhline(zflim[iz1],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
    
    axs.set_xlabel('x [km]')
    axs.set_ylabel('z [m]')
    axs.set_xlim(extent[0:2])
    plt.savefig(sp+'/concept.pdf',bbox_inches='tight',dpi=300)
    plt.show()
