#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 14:56:45 2021

@author: janssens
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import netCDF4 as nc
import matplotlib.animation as animation
from functions import getRad, lowPass, mean_mask

# lp = '/scratch-shared/janssens/bomex200aswitch/a2'
# lp = '/scratch-shared/janssens/eurec4a_old/eurec4a_mean_ssthet'
lp = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/BOMEXStability/bomex200_e12'
sp = lp+'/figs'

klp = 4
qlc = 1e-7

ds = nc.Dataset(lp+'/cape2d.001.nc')
ds1 = nc.Dataset(lp+'/profiles.001.nc')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)

extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000

circ_mask = np.zeros((xf.size,xf.size))
rad = getRad(circ_mask)
circ_mask[rad<=klp] = 1

# Calculate column-averaged density
zf = ds1['zt'][:].data
rhob = ds1['rhobf'][0,:].data
rho0 = np.trapz(rhob,zf)

#%% Plot twp, twppf and clouds at a single time step


tPlot = 24
fq=1e3

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
sc = axs[0].imshow(twpp*fq/rho0,extent=extent,vmin=-2*fq/rho0,vmax=2*fq/rho0,cmap='RdYlBu')
axs[0].set_xlabel('x [km]')
axs[0].set_ylabel('y [km]')

axs[1].imshow(twppf*fq/rho0,extent=extent,vmin=-2*fq/rho0,vmax=2*fq/rho0,cmap='RdYlBu')
axs[1].contour(twppf*fq/rho0,levels=[0],extent=extent,origin='upper')
axs[1].set_xlabel('x [km]')

# axs[2].contour(twppf,levels=[0],extent=extent,origin='upper',colors='white')
# axs[2].imshow(alb,extent=extent,vmin=0,vmax=1,cmap='gist_gray')
# axs[2].set_xlabel('x [km]')

pos1 = axs[1].get_position()
cbax = fig.add_axes([.92, pos1.ymin, 0.01, pos1.height])
# cbax = fig.add_axes([1, 0.1, 0.02, 0.85])
# cbax = fig.add_axes([-0.06, 0.1, 0.02, 0.85])
cb = fig.colorbar(sc, cax=cbax)
cb.ax.set_ylabel(r" $\langle q_t'\rangle$ [g/kg]", rotation=270, labelpad=15) #-65
# plt.tight_layout()
plt.savefig(sp+'/twpfluct.pdf',bbox_inches='tight',dpi=300)

#%% Plot the time evolution of twpp

tPlot = np.arange(6,16,3)

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
    cf = np.count_nonzero(cm) / cm.size

    buoycb = np.ma.getdata(ds.variables['buoycb'][it,:,:])

    sc1 = axs[0,j].imshow(twpp*fq/rho0, extent=extent,vmin=-2*fq/rho0,vmax=2*fq/rho0,cmap='RdYlBu')
    
    sc2 = axs[1,j].imshow(cm  , extent=extent,vmin=0 ,vmax=2000,cmap='gray')
    
    # sc3 = axs[2,j].imshow(buoycb  , extent=extent,vmin=-1 ,vmax=1,cmap='RdYlBu_r')
    
    if j > 1:
        axs[0,j].contour(twppf*fq/rho0,levels=[0],extent=extent,origin='upper',
                         linewidths=1,colors='black')
        axs[1,j].contour(twppf*fq/rho0,levels=[0],extent=extent,origin='upper',
                         linewidths=1,colors='white')
    
        # axs[2,j].contour(twppf,levels=[0],extent=extent,origin='upper',
                         # linewidths=1,colors='black')

    axs[1,j].set_xlabel('x [km]')
    # axs[2,j].set_xlabel('x [km]')
    axs[0,j].set_title('t = %.1f hr'%tPlot[j])
    if j == 0:
        axs[0,j].set_ylabel('y [km]')
        axs[1,j].set_ylabel('y [km]')

    if j == len(tPlot)-1:
        pos1 = axs[0,j].get_position()
        cbax1 = fig.add_axes([.91, pos1.ymin, 0.006, pos1.height])
        cb1 = fig.colorbar(sc1, cax=cbax1)
        cb1.ax.set_ylabel(r"$\langle q_t'\rangle$ [g/kg]", rotation=270, labelpad=15)
        
        pos2 = axs[1,j].get_position()
        cbax2 = fig.add_axes([.91, pos2.ymin, 0.006, pos2.height])
        cb2 = fig.colorbar(sc2, cax=cbax2)
        cb2.ax.set_ylabel(r"Cloud-top height [m]", rotation=270, labelpad=15)

        # pos3 = axs[2,j].get_position()
        # cbax3 = fig.add_axes([.92, pos2.ymin, 0.01, pos3.height])
        # cb3 = fig.colorbar(sc3, cax=cbax3)
        # cb3.ax.set_ylabel(r"Cloud-base buoyancy [K]", rotation=270, labelpad=15)
    axs[1,j].annotate('cloud fraction: %.2f'%cf, (0.05,1.05),xycoords='axes fraction')
plt.savefig(sp+'/twp_cld_evo.pdf', bbox_inches='tight',dpi=300)

#%% Make cth distribution
tPlot = np.arange(6,19,3)
itmin = np.argmin(abs(tPlot[0]-time))
itmax = np.argmin(abs(tPlot[-1]-time))

#fig,axs = plt.subplots(ncols=1,nrows=1,figsize=(5,5),squeeze=False)

cth = np.ma.getdata(ds.variables['cldtop'][itmin:itmax,:,:])
ax = sns.histplot(y=cth[cth>0],bins=np.arange(500,2040,40),element="poly", fill=False,stat='density',color='k')


#%% Make movie over time
tpltmin = 0.0
tpltmax = 24.0

itpltmin = np.where(time>=tpltmin)[0][0]
itpltmax = np.where(time<tpltmax)[0][-1]+1

tPlot = time[itpltmin:itpltmax]

def animate(i):
    ax.collections = []
    twpp = np.ma.getdata(ds.variables['twp'][itpltmin+i,:,:])
    twpp -= np.mean(twpp)
    twppf = lowPass(twpp, circ_mask)
    
    sc1 = ax.imshow(twpp, extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
    ax.contour(twppf,levels=[0],extent=extent,origin='upper',linewidths=1,colors='black')
    ax.set_xlabel('x [km]')
    ax.set_xlabel('y [km]')
    ax.set_title('t = %.1f hr'%tPlot[i])

fig, ax = plt.subplots(figsize=(5, 5))

ani = animation.FuncAnimation(fig, animate, interval=100, frames=len(tPlot))

ani.save(sp+"/twppf.mp4")

#%% Plot time evolution of fraction of twppf_moist
tPlot = np.arange(6,24,0.25)
alpha=0.5
lw=2
col_moist = plt.cm.RdYlBu(0.99)

twppf_moist_frac = np.zeros(tPlot.shape)
twppf_moist = np.zeros(tPlot.shape)
cf = np.zeros(tPlot.shape)

for j in range(len(tPlot)):
    it = np.argmin(abs(tPlot[j]-time))
    twpp = np.ma.getdata(ds.variables['twp'][it,:,:])
    twpp -= np.mean(twpp)
    twppf = lowPass(twpp, circ_mask)
    mask_moist = np.zeros(twppf.shape)
    mask_moist[twppf - np.mean(twppf) > 0] = 1
    twppf_moist_frac[j] = np.count_nonzero(mask_moist) / mask_moist.size
    twppf_moist[j] = mean_mask(twppf, mask_moist)
    
    cm = np.ma.getdata(ds.variables['cldtop'][it,:,:])
    cf[j] = np.count_nonzero(cm) / cm.size

fig = plt.figure(figsize=(5,10/3)); axs = plt.gca()
# axs.plot(tPlot,twppf_moist,c=col_moist,lw=lw,alpha=alpha,label=labs[i])
axs.plot(tPlot,twppf_moist_frac,c=col_moist,lw=lw,alpha=alpha,label='Fraction of moist, mesoscale columns')

ax2 = axs.twinx()
ax2.plot(tPlot,cf,c='k',lw=lw,alpha=alpha,label='Cloud fraction')
plt.show()
#%% Plot time evolution of twp, for a number of simulations

# time 50m res:np.arange(600,115200,600)

lps = ['/scratch-shared/janssens/bomex200aswitch/a2',
       '/scratch-shared/janssens/bomex100',
        '/scratch-shared/janssens/bomex50',
        '/scratch-shared/janssens/bomex200a5',
       ]
sp = lps[-1]+'/figs'

labs = [r'$\Delta x = 200m$',
         r'$\Delta x = 100m$',
          r'$\Delta x = 50m$',
         r'$O(5)$ advection']
ls = ['-','--',':','-.']

tmin = 6.
tmax = [36., 
        24.,
        36.,
        36.,
        36.]

klp = 4
qlc = 1e-7

alpha=0.5
lw=2
col_moist = plt.cm.RdYlBu(0.99)
col_dry = plt.cm.RdYlBu(0)

# fig,axs = plt.subplots(ncols=len(tPlot),nrows=len(lps),figsize=(3.5*len(tPlot),3.5*len(lps)),
#                        sharex=True,sharey=True,squeeze=False)

f1 = plt.figure(figsize=(5,10/3)); axs1 = plt.gca()

for i in range(len(lps)):
    ds = nc.Dataset(lps[i]+'/cape2d.001.nc')
    
    time  = np.ma.getdata(ds.variables['time'][:]) / 3600
    # xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
    # xh    = np.ma.getdata(ds.variables['xm'][:]) # Cell edges (h in mhh)
    # yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)
    # yh    = np.ma.getdata(ds.variables['ym'][:]) # Cell edges (h in mhh)
    
    itpltmin = np.where(time>=tmin)[0][0]
    itpltmax = np.where(time<tmax[i])[0][-1]+1
    plttime_var = np.arange(itpltmin,itpltmax,1)

    
    # extent = np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000
    sz = ds.dimensions['xt'].size
    circ_mask = np.zeros((sz,sz))
    rad = getRad(circ_mask)
    circ_mask[rad<=klp] = 1
    
    # for j in range(len(tPlot)):
    #     it = np.argmin(abs(tPlot[j]-time))
    #     twpp = np.ma.getdata(ds.variables['twp'][it,:,:])
    #     twpp -= np.mean(twpp)
    #     twppf = lowPass(twpp, circ_mask)
        
    #     sc1 = axs[i,j].imshow(twpp, extent=extent,vmin=-2,vmax=2,cmap='RdYlBu')
        
    #     if j > 1:
    #         axs[i,j].contour(twppf,levels=[0],extent=extent,origin='upper',
    #                          linewidths=0.5,colors='black')
    #     if i == 0:
    #         axs[0,j].set_title('t = %.1f hr'%tPlot[j])
    #     if i == len(lps)-1:
    #         axs[1,j].set_xlabel('x [km]')
    #     if j == 0:
    #         axs[i,j].set_ylabel('y [km]')
    
    #     if j == len(tPlot)-1:
    #         pos1 = axs[1,j].get_position()
    #         cbax1 = fig.add_axes([.92, pos1.ymin, 0.01, pos1.height])
    #         cb1 = fig.colorbar(sc1, cax=cbax1)
    #         cb1.ax.set_ylabel(r"$TWP'$ [kg/kg/m$^2$]", rotation=270, labelpad=15)
    
    twppf_moist = np.zeros(len(plttime_var))
    twppf_dry = np.zeros(len(plttime_var))
    for j in range(len(plttime_var)):
        twpp = np.ma.getdata(ds.variables['twp'][plttime_var[j],:,:])
        twpp -= np.mean(twpp)
        twppf = lowPass(twpp, circ_mask)
        mask_moist = np.zeros(twppf.shape)
        mask_moist[twppf > 0] = 1
        mask_dry = 1 - mask_moist
    
        twppf_moist[j] = mean_mask(twppf, mask_moist)
        twppf_dry[j] = mean_mask(twppf, mask_dry)
    
    axs1.plot(time[plttime_var],twppf_moist,c=col_moist,linestyle=ls[i],lw=lw,alpha=alpha,label=labs[i])
    axs1.plot(time[plttime_var],twppf_dry,c=col_dry,linestyle=ls[i],lw=lw,alpha=alpha)
axs1.set_xlabel('Time [hr]')
axs1.set_ylabel(r"$TWP_m'$ [kg/m$^2$]")
axs1.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig(sp+'/twp_evo_num.pdf',bbox_inches='tight')
