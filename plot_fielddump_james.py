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

lp = '/scratch-shared/janssens/bomex200_e12'
sp = lp+'/figs'
it = 63
izmin = 0
izmax = 80
klp = 4
ix = 372
iz = 37
iz0 = 15
iz1 = 37
fq = 1000 # kg/kg => g/kg
fwql = 7*300 # kg/kg m/s => K m/s
generate_slice = False

if generate_slice:
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
    
    # Vertical differences
    dzf = np.zeros(zh.shape)
    dzf[:-1] = zh[1:] - zh[:-1] # First value is difference top 1st cell and surface
    dzf[-1] = dzf[-2]
    
    dzh = np.zeros(zf.shape)
    dzh[1:] = zf[1:] - zf[:-1] # First value is difference mid 1st cell and mid 1st cell below ground
    dzh[0] = 2*zf[1]
    
    delta = (dx*dy*np.diff(zh))**(1./3)
    
    dzflim = dzf[izmin:izmax]
    dzhlim = dzh[izmin:izmax]
    
    zflim = zf[izmin:izmax]
    zhlim = zh[izmin:izmax]
    
    # Mask for low-[ass filtering
    circ_mask = np.zeros((xf.size,xf.size))
    rad = getRad(circ_mask)
    circ_mask[rad<=klp] = 1
    
    it1d = np.argmin(np.abs(time1d/3600 - time[it]))
    
    # 1D fields
    rhobfi = rhobf[it1d,izmin:izmax]
    rhobhi = rhobh[it1d,izmin:izmax]
    
    # 3D fields
    qt  = np.ma.getdata(ds.variables['qt'][it,izmin:izmax,:,:])
    u = np.ma.getdata(ds.variables['u'][it,izmin:izmax,:,:])
    v = np.ma.getdata(ds.variables['v'][it,izmin:izmax,:,:])
    wh = np.ma.getdata(ds.variables['w'][it,izmin:izmax+1,:,:])
    thl =  np.ma.getdata(ds.variables['thl'][it,izmin:izmax,:,:])
    ql = np.ma.getdata(ds.variables['ql'][it,izmin:izmax,:,:])
    if 'p' in ds.variables.keys():
        p = np.ma.getdata(ds.variables['p'][it,izmin:izmax,:,:])
    
    # w  h-> f-levels
    wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
    wh = wh[:-1,:,:]
    
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
    
    div_wql = ddzwx_2nd(wh, qlp, dzflim, dzhlim, rhobf=rhobfi)
    div_wql_av = np.mean(div_wql,axis=(1,2))
    div_wqlf = lowPass(div_wql, circ_mask)

    # Store necesasry fields to make figures
    qtplt = qt[:,ix,:]
    qlplt = ql[:,ix,:]
    qtpfplt = qtpf[:,ix,:]
    upfplt = upf[:,ix,:]
    wffplt = wff[:,ix,:]
    qsplt = qs[:,ix,:]
    wqlfplt = wqlf[:,ix,:]
    wthlvpplt = wthlvp[iz,ix,:]-np.mean(wthlvp[iz,:,:])
    wthlvpfplt = wthlvpf[iz,ix,:]-np.mean(wthlvp[iz,:,:])
    wqlpplt = 7*thl_av[iz]*(wqlp[iz,ix,:]-np.mean(wqlp[iz,:,:]))
    wqlpfplt = 7*thl_av[iz]*(wqlpf[iz,ix,:]-np.mean(wqlp[iz,:,:]))
    
    np.save(lp+'/fig_data/qtplf.npy',qtplt)
    np.save(lp+'/fig_data/qlplt.npy',qlplt)
    np.save(lp+'/fig_data/qtpfplt.npy',qtpfplt)
    np.save(lp+'/fig_data/upfplt.npy',upfplt)
    np.save(lp+'/fig_data/wffplt.npy',wffplt)
    np.save(lp+'/fig_data/zflim.npy',zflim)
    np.save(lp+'/fig_data/xf.npy',xf)
    np.save(lp+'/fig_data/qsplt.npy',qsplt)
    np.save(lp+'/fig_data/wqlfplt.npy',wqlfplt)
    np.save(lp+'/fig_data/wthlvpplt.npy',wthlvpplt)
    np.save(lp+'/fig_data/wthlvpfplt.npy',wthlvpfplt)
    np.save(lp+'/fig_data/wqlpplt.npy',wqlpplt)
    np.save(lp+'/fig_data/wqlpfplt.npy',wqlpfplt)
    

if not generate_slice:
    qtplt = np.load(lp+'/fig_data/qtplf.npy')
    qlplt = np.load(lp+'/fig_data/qlplt.npy')
    qtpfplt = np.load(lp+'/fig_data/qtpfplt.npy')
    upfplt = np.load(lp+'/fig_data/upfplt.npy')
    wffplt = np.load(lp+'/fig_data/wffplt.npy')
    zflim = np.load(lp+'/fig_data/zflim.npy')
    xf = np.load(lp+'/fig_data/xf.npy')
    qsplt = np.load(lp+'/fig_data/qsplt.npy')
    wqlfplt = np.load(lp+'/fig_data/wqlfplt.npy')
    wthlvpplt = np.load(lp+'/fig_data/wthlvpplt.npy')
    wthlvpfplt = np.load(lp+'/fig_data/wthlvpfplt.npy')
    wqlpplt = np.load(lp+'/fig_data/wqlpplt.npy')
    wqlpfplt = np.load(lp+'/fig_data/wqlpfplt.npy')


#%% First inset
extent=[20,70,20,2500]
fig = plt.figure(figsize=(5,4)); axs = plt.gca()
sc0 = axs.imshow(np.flipud(qtplt[:62,100:350])*fq,
                        extent=extent,
                        aspect='auto',cmap='RdYlBu',
                        vmin=0.01*fq,vmax=0.018*fq)
pos = axs.get_position()
cbax0=fig.add_axes([.92, pos.ymin, 0.01, pos.height])
cb0 = fig.colorbar(sc0, cax=cbax0)
cb0.ax.set_ylabel(r"Total moisture [g/kg]", rotation=270, labelpad=15)

axs.contour(qlplt[:62,100:350],levels=[1e-7],extent=extent,origin='lower',linewidths=1,colors='black')

axs.contour(qtpfplt[:62,100:350],levels=[0.0004],extent=extent,origin='lower',linewidths=1,colors='black')# ,linestyles='dashed')
xfstr = np.linspace(20,102300,512)/1000
[X,Y] = np.meshgrid(xfstr,np.arange(20,2500,40))
speed = np.sqrt((upfplt[:62,:]/1000)**2+wffplt[:62,:]**2)
lws = np.maximum(.75,3.*speed/np.max(speed))
st = axs.streamplot(X,Y,upfplt[:62,:]/1000,wffplt[:62,:],linewidth=lws, density=1,color='grey')
axs.axhline(zflim[iz0],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
axs.axhline(zflim[iz1],linestyle='--',color='k',alpha=0.8,linewidth=0.5)

axs.set_xlabel('x [km]')
axs.set_ylabel('z [m]')
axs.set_xlim(extent[0:2])
plt.savefig(sp+'/concept_1.svg',bbox_inches='tight',dpi=300)
plt.show()

#%% Second inset
extent=[30,50,zflim[10],2500]
fig = plt.figure(figsize=(5,4)); axs = plt.gca()
sc0 = axs.imshow(np.flipud(qtplt[10:62,150:275]/qsplt[10:62,150:275]),
                        extent=extent,
                        aspect='auto',cmap='YlGnBu',
                        vmin=0.8,vmax=1.1)
pos = axs.get_position()
cbax0=fig.add_axes([.92, pos.ymin, 0.01, pos.height])
cb0 = fig.colorbar(sc0, cax=cbax0)
cb0.ax.set_ylabel(r"Relative humidity [g/kg]", rotation=270, labelpad=15,fontsize=16)
cb0.ax.tick_params(labelsize=14)

sc1 = axs.contour(wqlfplt[10:62,150:275]*fwql,extent=extent,origin='lower',linewidths=1.5,cmap='Purples',levels=np.linspace(0.00005*fwql,0.0002*fwql,8))
cbax1=fig.add_axes([1.09, pos.ymin, 0.01, pos.height])
norm = colors.Normalize(vmin=sc1.cvalues.min(), vmax=sc1.cvalues.max())
sm = plt.cm.ScalarMappable(norm=norm, cmap=sc1.cmap)
sm.set_array([])
cb1 = fig.colorbar(sm, cax=cbax1)
cb1.ax.set_ylabel(r"Liquid water flux [Km/s]", rotation=270, labelpad=15, fontsize=16)
cb1.ax.tick_params(labelsize=14)
cb1.locator = ticker.MaxNLocator(nbins=6)
cb1.update_ticks()
axs.contour(qlplt[10:62,150:275],levels=[1e-7],extent=extent,origin='lower',linewidths=1,colors='black')

axs.axhline(zflim[iz0],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
axs.axhline(zflim[iz1],linestyle='--',color='k',alpha=0.8,linewidth=0.5)

axs.set_xlabel('x [km]',fontsize=16)
axs.set_ylabel('z [m]',fontsize=16)
axs.tick_params(axis='both', labelsize=14)

axs.set_xlim(extent[0:2])
plt.savefig(sp+'/concept_2.svg',bbox_inches='tight',dpi=300)
plt.show()

#%% Flux origins
extent=[xf[0]/1000,xf[-1]/1000,zflim[0],zflim[-1]-1]
fig,axs = plt.subplots(figsize=(10,4),nrows=2,sharex=True,gridspec_kw={'height_ratios':[1,3]})
sc0 = axs[1].imshow(np.flipud(qtplt[:-1,:]-qsplt),
                        extent=extent,
                        aspect='auto',cmap='YlGnBu',
                        vmin=-0.01,vmax=0.003)
pos = axs[1].get_position()
cbax0=fig.add_axes([.92, pos.ymin, 0.007, pos.height])
cb0 = fig.colorbar(sc0, cax=cbax0)
cb0.ax.set_ylabel(r"$q_t- q_s$ [kg/kg]", rotation=270, labelpad=15)

sc1 = axs[1].contour(qtpfplt[:-1,:],extent=extent,origin='lower',linewidths=0.75,cmap='Blues',levels=np.linspace(0.0002,0.003,6))
cbax1=fig.add_axes([1.02, pos.ymin, 0.007, pos.height])
norm = colors.Normalize(vmin=sc1.cvalues.min(), vmax=sc1.cvalues.max())
sm = plt.cm.ScalarMappable(norm=norm, cmap=sc1.cmap)
sm.set_array([])
cb1 = fig.colorbar(sm, cax=cbax1)
cb1.ax.set_ylabel(r"$q_{t_m}'$ [kg/kg]", rotation=270, labelpad=15)
cb1.locator = ticker.MaxNLocator(nbins=6)
cb1.update_ticks()
axs[1].contour(qlplt[:-1,:],levels=[1e-7],extent=extent,origin='lower',linewidths=0.5,colors='black')
axs[1].axhline(zflim[iz],linestyle='--',color='k',alpha=0.8,linewidth=0.5)
axs[1].set_xlabel('x [km]')
axs[1].set_ylabel('z [m]')

axs[0].plot(xf/1000,wthlvpplt,c='maroon',label=r"$F_{\theta_{lv}'}$")
axs[0].plot(xf/1000,wthlvpfplt,c='maroon',linestyle='dashed',label=r"$F_{{\theta_{lv}'}_m}$")
axs[0].plot(xf/1000,-wqlpplt,c='seagreen',label=r"$-7\overline{\theta_l}F_{q_l'}$")
axs[0].plot(xf/1000,-wqlpfplt,c='seagreen',linestyle='dashed',label=r"$-7\overline{\theta_l}F_{{q_l'}_m}$")
axs[0].set_ylim((-1,0))
axs[0].set_ylabel(r"$F_{\theta_{lv}}$, [Km/s]")
axs[0].legend(bbox_to_anchor=(1,1),loc='upper left',ncol=2)

plt.savefig(sp+'/structure.pdf',bbox_inches='tight',dpi=300)
plt.show()

