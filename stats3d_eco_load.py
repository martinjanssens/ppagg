#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 11:18:54 2021

@author: janssens
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.optimize import curve_fit
from skimage.measure import block_reduce


# Run specifics (now have to match the run you've done to create the files to be loaded)
itmin = 11#23
itmax = 40
di    = 1
izmin = 0
izmax = 80

klp = 4

lp = '/scratch-shared/janssens/bomex100'
ds = nc.Dataset(lp+'/fielddump.001.nc')
ds1= nc.Dataset(lp+'/profiles.001.nc')
ds0= nc.Dataset(lp+'/tmser.001.nc')
ilp = np.loadtxt(lp+'/lscale.inp.001')

time  = np.ma.getdata(ds.variables['time'][:]) / 3600
zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])

dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain

# Larger-scale subsidence
wfls = ilp[:,3]


plttime = np.arange(itmin, itmax, di)
zflim = zf[izmin:izmax]

qtpf_moist_time = np.load(lp+'/qtpf_moist_time.npy')
qtpf_dry_time = np.load(lp+'/qtpf_dry_time.npy')
qtpf_prod_moist_time = np.load(lp+'/qtpf_prod_moist_time.npy')
qtpf_prod_dry_time = np.load(lp+'/qtpf_prod_dry_time.npy')
qtpf_prod_moist_wex_time = np.load(lp+'/qtpf_prod_moist_wex_time.npy')
qtpf_prod_dry_wex_time = np.load(lp+'/qtpf_prod_dry_wex_time.npy')
qtpf_vdiv_moist_time = np.load(lp+'/qtpf_vdiv_moist_time.npy')
qtpf_vdiv_dry_time = np.load(lp+'/qtpf_vdiv_dry_time.npy')
qtpf_hdiv_moist_time = np.load(lp+'/qtpf_hdiv_moist_time.npy')
qtpf_hdiv_dry_time = np.load(lp+'/qtpf_hdiv_dry_time.npy')
qtpf_subs_moist_time = np.load(lp+'/qtpf_subs_moist_time.npy')
qtpf_subs_dry_time = np.load(lp+'/qtpf_subs_dry_time.npy')

thlvpf_moist_time = np.load(lp+'/thlvpf_moist_time.npy')
thlvpf_dry_time = np.load(lp+'/thlvpf_dry_time.npy')
thlvpf_prod_moist_time = np.load(lp+'/thlvpf_prod_moist_time.npy')
thlvpf_prod_dry_time = np.load(lp+'/thlvpf_prod_dry_time.npy')
thlvpf_vdiv_moist_time = np.load(lp+'/thlvpf_vdiv_moist_time.npy')
thlvpf_vdiv_dry_time = np.load(lp+'/thlvpf_vdiv_dry_time.npy')
thlvpf_hdiv_moist_time = np.load(lp+'/thlvpf_hdiv_moist_time.npy')
thlvpf_hdiv_dry_time = np.load(lp+'/thlvpf_hdiv_dry_time.npy')
thlvpf_subs_moist_time = np.load(lp+'/thlvpf_subs_moist_time.npy')
thlvpf_subs_dry_time = np.load(lp+'/thlvpf_subs_dry_time.npy')

thlpf_moist_time = np.load(lp+'/thlpf_moist_time.npy')
thlpf_dry_time = np.load(lp+'/thlpf_dry_time.npy')

wff_moist_time = np.load(lp+'/wff_moist_time.npy')
wff_dry_time = np.load(lp+'/wff_dry_time.npy')

qlpf_moist_time = np.load(lp+'/qlpf_moist_time.npy') 
qlpf_dry_time = np.load(lp+'/qlpf_dry_time.npy')

#%% Plotprofiles of  mesoscale-filtered variables in time

fig,axs = plt.subplots(ncols=5,sharey=True,figsize=(10,5))
for i in range(len(plttime)):
    col = plt.cm.cubehelix(i/len(plttime))
     
    axs[0].plot(qtpf_moist_time[i,:], zflim, color=col,linestyle='-')
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].set_xlabel(r"$q_t'$")
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[1].plot(thlvpf_moist_time[i,:], zflim, color=col,linestyle='-')
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].set_xlabel(r"$\theta_{lv}'$")
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[2].plot(wff_moist_time[i,:], zflim,color=col,linestyle='-')
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].set_xlabel(r"$w'$")
    axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[3].plot(thlpf_moist_time[i,:], zflim,color=col,linestyle='-')
    axs[3].axvline(0,color='gray',linestyle='dotted')
    axs[3].set_xlabel(r"$\theta_l'$")
    axs[3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[4].plot(qlpf_moist_time[i,:], zflim, label='t=%.2f'%time[plttime[i]],color=col,linestyle='-')
    axs[4].axvline(0,color='gray',linestyle='dotted')
    axs[4].set_xlabel(r"$q_l'$")
    axs[4].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[4].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

#%%

# Average budget contributions over time dimension
itav = 16 # four time steps (15 min for bomex -> 1 hour average)

qtpfmn_prod_moist_wex = block_reduce(qtpf_prod_moist_wex_time,(itav,1),func=np.mean)
qtpfmn_vdiv_moist = block_reduce(qtpf_vdiv_moist_time,(itav,1),func=np.mean)
qtpfmn_hdiv_moist = block_reduce(qtpf_hdiv_moist_time,(itav,1),func=np.mean)
qtpfmn_subs_moist = block_reduce(qtpf_subs_moist_time,(itav,1),func=np.mean)
# qtpfmn_diff_moist = block_reduce(qtpf_diff_moist_time,(itav,1),func=np.mean)
qtpfmn_budg_moist = (-qtpfmn_prod_moist_wex[:,1:-1] - qtpfmn_vdiv_moist[:,1:-1]
                     -qtpfmn_hdiv_moist[:,1:-1] - qtpfmn_subs_moist[:,1:-1])
                     # +qtpfmn_diff_moist)
qtpfmn_prod_dry_wex = block_reduce(qtpf_prod_dry_wex_time,(itav,1),func=np.mean)
qtpfmn_vdiv_dry = block_reduce(qtpf_vdiv_dry_time,(itav,1),func=np.mean)
qtpfmn_hdiv_dry = block_reduce(qtpf_hdiv_dry_time,(itav,1),func=np.mean)
qtpfmn_subs_dry = block_reduce(qtpf_subs_dry_time,(itav,1),func=np.mean)
# qtpfmn_diff_dry = block_reduce(qtpf_diff_dry_time,(itav,1),func=np.mean)
qtpfmn_budg_dry = (-qtpfmn_prod_dry_wex[:,1:-1] - qtpfmn_vdiv_dry[:,1:-1]
                     -qtpfmn_hdiv_dry[:,1:-1] - qtpfmn_subs_dry[:,1:-1])
                     # +qtpfmn_diff_dry)

thlvpfmn_prod_moist = block_reduce(thlvpf_prod_moist_time,(itav,1),func=np.mean)
thlvpfmn_vdiv_moist = block_reduce(thlvpf_vdiv_moist_time,(itav,1),func=np.mean)
thlvpfmn_hdiv_moist = block_reduce(thlvpf_hdiv_moist_time,(itav,1),func=np.mean)
thlvpfmn_subs_moist = block_reduce(thlvpf_subs_moist_time,(itav,1),func=np.mean)
# thlvpfmn_diff_moist = block_reduce(thlvpf_diff_moist_time,(itav,1),func=np.mean)
thlvpfmn_budg_moist = (-thlvpfmn_prod_moist[:,1:-1] - thlvpfmn_vdiv_moist[:,1:-1]
                     -thlvpfmn_hdiv_moist[:,1:-1] - thlvpfmn_subs_moist[:,1:-1])
                     # +thlvpfmn_diff_moist)
thlvpfmn_prod_dry = block_reduce(thlvpf_prod_dry_time,(itav,1),func=np.mean)
thlvpfmn_vdiv_dry = block_reduce(thlvpf_vdiv_dry_time,(itav,1),func=np.mean)
thlvpfmn_hdiv_dry = block_reduce(thlvpf_hdiv_dry_time,(itav,1),func=np.mean)
thlvpfmn_subs_dry = block_reduce(thlvpf_subs_dry_time,(itav,1),func=np.mean)
# thlvpfmn_diff_dry = block_reduce(thlvpf_diff_dry_time,(itav,1),func=np.mean)
thlvpfmn_budg_dry = (-thlvpfmn_prod_dry[:,1:-1] - thlvpfmn_vdiv_dry[:,1:-1]
                     -thlvpfmn_hdiv_dry[:,1:-1] - thlvpfmn_subs_dry[:,1:-1])
                     # +thlvpfmn_diff_dry)

# Budget terms
terms = [r"$\frac{\partial\langle\tilde{q_t'}\rangle}{\partial t}$",
         r"$-\tilde{w'}\frac{\partial \overline{q_t}}{\partial z}$",
         r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'''q_t'''}-\overline{w'q_t'}\right)\right)$",
         r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'q_t'}\right)$",
         r"$-\overline{w_{LS}}\frac{\partial \tilde{q_t'}}{\partial z}$",
         r"$\widetilde{\frac{\partial}{\partial x_j}\left(K_h\frac{\partial q_t'}{\partial x_j}\right)}+\widetilde{\frac{\partial}{\partial x_j}\left(K_h'\frac{\partial \overline{q_t}}{\partial x_j}\right)}$"
         ]

# Budget at kth avergaed time step
k=-1
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(qtpfmn_budg_moist[k,:], zflim[2:-2],c='midnightblue')
axs[0].plot(-qtpfmn_prod_moist_wex[k,:], zflim[1:-1],c='darkseagreen')
axs[0].plot(-qtpfmn_vdiv_moist[k,:], zflim[1:-1],c='maroon')
axs[0].plot(-qtpfmn_hdiv_moist[k,:], zflim[1:-1],c='peachpuff')
axs[0].plot(-qtpfmn_subs_moist[k,:], zflim[1:-1],c='olive')
# axs[0].plot(qtpfmn_diff_moist[k,:], zflim[2:-2],c='skyblue')
axs[0].set_xlabel(r"Contribution to $\frac{\partial\tilde{q_t'}}{\partial t}$")

axs[1].plot(qtpfmn_budg_dry[k,:], zflim[2:-2],c='midnightblue',label=terms[0])
axs[1].plot(-qtpfmn_prod_dry_wex[k,:], zflim[1:-1],c='darkseagreen',label=terms[1])
axs[1].plot(-qtpfmn_vdiv_dry[k,:], zflim[1:-1],c='maroon',label=terms[2])
axs[1].plot(-qtpfmn_hdiv_dry[k,:], zflim[1:-1],c='peachpuff',label=terms[3])
axs[1].plot(-qtpfmn_subs_dry[k,:], zflim[1:-1],c='olive',label=terms[4])
# axs[1].plot(qtpfmn_diff_dry[k,:], zflim[2:-2],c='skyblue',label=terms[5])
axs[1].set_xlabel(r"Contribution to $\frac{\partial\tilde{q_t'}}{\partial t}$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(thlvpfmn_budg_moist[k,:], zflim[2:-2],c='midnightblue')
axs[0].plot(-thlvpfmn_prod_moist[k,:], zflim[1:-1],c='darkseagreen')
axs[0].plot(-thlvpfmn_vdiv_moist[k,:], zflim[1:-1],c='maroon')
axs[0].plot(-thlvpfmn_hdiv_moist[k,:], zflim[1:-1],c='peachpuff')
axs[0].plot(-thlvpfmn_subs_moist[k,:], zflim[1:-1],c='olive')
# axs[0].plot(thlvpfmn_diff_moist[k,:], zflim[2:-2],c='skyblue')
axs[0].set_xlabel(r"Contribution to $\frac{\partial\tilde{\theta_{lv}'}}{\partial t}$")

axs[1].plot(thlvpfmn_budg_dry[k,:], zflim[2:-2],c='midnightblue',label='Tendency')
axs[1].plot(-thlvpfmn_prod_dry[k,:], zflim[1:-1],c='darkseagreen',label='Gradient production')
axs[1].plot(-thlvpfmn_vdiv_dry[k,:], zflim[1:-1],c='maroon',label='Anomalous vertical flux divergence')
axs[1].plot(-thlvpfmn_hdiv_dry[k,:], zflim[1:-1],c='peachpuff',label='Horizontal divergence')
axs[1].plot(-thlvpfmn_subs_dry[k,:], zflim[1:-1],c='olive',label='Subsidence')
# axs[1].plot(thlvpfmn_diff_dry[k,:], zflim[2:-2],c='skyblue',label='SFS diffusion')
axs[1].set_xlabel(r"Contribution to $\frac{\partial\tilde{\theta_{lv}'}}{\partial t}$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

#%% Vertically integrated statistics
def vint(field,rhob,z,tmin=0,tmax=-1):
    
    if len(field.shape) == 3:
        var = np.trapz(rhob[:,np.newaxis,np.newaxis]*field[:,:,:],z,axis=0)
    elif len(field.shape) == 4:
        var = np.trapz(rhob[np.newaxis,:,np.newaxis,np.newaxis]*
                       field[tmin:tmax,:,:,:],z,axis=1)
    elif len(field.shape) == 2:
        var = np.trapz(rhob[np.newaxis,:]*field[tmin:tmax,:],z,axis=1)
    elif len(field.shape) == 1:
        var = np.trapz(rhob*field,z)
    return var

it1d = np.argmin(np.abs(time1d/3600 - time[plttime[-1]]))
    
# 1D fields
rhobfi = rhobf[it1d,izmin:izmax]

izminint = 8 # 500m for BOMEX, above the big gradient jump
qtpfi_moist = vint(qtpf_moist_time,rhobfi,zflim,tmax=len(plttime))
qtpfi_dry = vint(qtpf_dry_time,rhobfi,zflim,tmax=len(plttime))

# Moistening gradient production per simplified WTG budget
qtpfi_prod_moist = vint(qtpf_prod_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
qtpfi_prod_dry = vint(qtpf_prod_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

qtpfi_prod_wex_moist = vint(qtpf_prod_moist_wex_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
qtpfi_prod_wex_dry = vint(qtpf_prod_dry_wex_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Moistening through anomalous vertical small-scale fluxes
# FIXME offset zf in integration by 1 from field
qtpfi_vdiv_moist = vint(qtpf_vdiv_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
qtpfi_vdiv_dry = vint(qtpf_vdiv_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Moistening through horizontal advection
qtpfi_hdiv_moist = vint(qtpf_hdiv_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
qtpfi_hdiv_dry = vint(qtpf_hdiv_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Moistening through subsidence
qtpfi_subs_moist = vint(qtpf_subs_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
qtpfi_subs_dry = vint(qtpf_subs_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Moistening through SFS diffusion
# qtpfi_diff_moist = vint(qtpf_diff_moist_time,rhobfi,zflim,tmax=len(plttime))
# qtpfi_diff_dry = vint(qtpf_diff_dry_time,rhobfi,zflim,tmax=len(plttime))

# Fit the moisture fluctuation's evolution
[exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime]*3600, 
                                       qtpfi_moist,
                                       p0=[1,1e-5])

[exp_dry,fac_dry], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime]*3600, 
                                       qtpfi_dry,
                                       p0=[-1,-1e-5])

# And differentiate to estimate its tendency
qtpfi_tend_moist = fac_moist*exp_moist*(time[plttime]*3600)**(exp_moist-1)
qtpfi_tend_dry = fac_dry*exp_dry*(time[plttime]*3600)**(exp_dry-1)

# Estimate residual
qtpfi_resid_moist = qtpfi_tend_moist + qtpfi_prod_wex_moist + qtpfi_vdiv_moist + qtpfi_hdiv_moist + qtpfi_subs_moist #- qtpfi_diff_moist
qtpfi_resid_dry = qtpfi_tend_dry + qtpfi_prod_wex_dry + qtpfi_vdiv_dry + qtpfi_hdiv_dry + qtpfi_subs_dry #- qtpfi_diff_dry

# Temporal plot
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
axs[0].plot(time[plttime],qtpfi_tend_moist,c='midnightblue')
# axs[0].plot(time[plttime],qtpfi_prod_moist,c='darkseagreen')
axs[0].plot(time[plttime],-qtpfi_prod_wex_moist,c='darkseagreen')
axs[0].plot(time[plttime],-qtpfi_vdiv_moist,c='maroon')
axs[0].plot(time[plttime],-qtpfi_hdiv_moist,c='peachpuff')
axs[0].plot(time[plttime],-qtpfi_subs_moist,c='olive')
# axs[0].plot(time[plttime],qtpfi_diff_moist,c='skyblue')
axs[0].plot(time[plttime],qtpfi_resid_moist,c='slategray')
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist region')

axs[1].plot(time[plttime],qtpfi_tend_dry,c='midnightblue',label=terms[0])
# axs[1].plot(time[plttime],qtpfi_prod_dry,c='darkseagreen',label=r"$F_{\langle\tilde{q_t'}\rangle}$")
axs[1].plot(time[plttime],-qtpfi_prod_wex_dry,c='darkseagreen',label=terms[1])
axs[1].plot(time[plttime],-qtpfi_vdiv_dry,c='maroon',label=terms[2])
axs[1].plot(time[plttime],-qtpfi_hdiv_dry,c='peachpuff',label=terms[3])
axs[1].plot(time[plttime],-qtpfi_subs_dry,c='olive',label=terms[4])
# axs[1].plot(time[plttime],qtpfi_diff_dry,c='skyblue',label=terms[5])
axs[1].plot(time[plttime],qtpfi_resid_dry,c='slategray',label=r"Residual")
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry region')

axs[0].set_ylabel('Large-scale moistening rate [kg/kg/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

### And now for thlv
thlvpfi_moist = vint(thlvpf_moist_time,rhobfi,zflim,tmax=len(plttime))
thlvpfi_dry = vint(thlvpf_dry_time,rhobfi,zflim,tmax=len(plttime))

# thlv gradient production
thlvpfi_prod_moist = vint(thlvpf_prod_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
thlvpfi_prod_dry = vint(thlvpf_prod_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# heating through anomalous vertical fluxes
# FIXME offset zf in integration by 1 from field
thlvpfi_vdiv_moist = vint(thlvpf_vdiv_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
thlvpfi_vdiv_dry = vint(thlvpf_vdiv_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Heating through horizontal advection
thlvpfi_hdiv_moist = vint(thlvpf_hdiv_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
thlvpfi_hdiv_dry = vint(thlvpf_hdiv_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Heating through subsidence
thlvpfi_subs_moist = vint(thlvpf_subs_moist_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))
thlvpfi_subs_dry = vint(thlvpf_subs_dry_time,rhobfi[1:-1],zflim[1:-1],tmax=len(plttime))

# Moistening through SFS diffusion
# thlvpfi_diff_moist = vint(thlvpf_diff_moist_time,rhobfi,zflim,tmax=len(plttime))
# thlvpfi_diff_dry = vint(thlvpf_diff_dry_time,rhobfi,zflim,tmax=len(plttime))


# Fit the haeting fluctuation's evolution
[exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime]*3600, 
                                       thlvpfi_moist,
                                       p0=[1,0])

[exp_dry,fac_dry], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime]*3600, 
                                       thlvpfi_dry,
                                       p0=[-3,0])

# And differentiate to estimate its tendency
thlvpfi_tend_moist = fac_moist*exp_moist*(time[plttime]*3600)**(exp_moist-1)
thlvpfi_tend_dry = fac_dry*exp_dry*(time[plttime]*3600)**(exp_dry-1)

# Estimate residual
thlvpfi_resid_moist = thlvpfi_tend_moist + thlvpfi_prod_moist + thlvpfi_vdiv_moist + thlvpfi_hdiv_moist + thlvpfi_subs_moist
thlvpfi_resid_dry = thlvpfi_tend_dry + thlvpfi_prod_dry + thlvpfi_vdiv_dry + thlvpfi_hdiv_dry + thlvpfi_subs_dry

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
axs[0].plot(time[plttime],thlvpfi_tend_moist,c='midnightblue')
# axs[0].plot(time[plttime],qtpfi_prod_moist,c='darkseagreen')
axs[0].plot(time[plttime],-thlvpfi_prod_moist,c='maroon')
axs[0].plot(time[plttime],-thlvpfi_vdiv_moist,c='peachpuff')
axs[0].plot(time[plttime],-thlvpfi_hdiv_moist,c='olive')
axs[0].plot(time[plttime],-thlvpfi_subs_moist,c='skyblue')
axs[0].plot(time[plttime],thlvpfi_resid_moist,c='slategray')
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist region')

axs[1].plot(time[plttime],qtpfi_tend_dry,c='midnightblue',label=r"$\frac{\partial\langle\tilde{\theta_{lv}'}\rangle}{\partial t}$")
# axs[1].plot(time[plttime],qtpfi_prod_dry,c='darkseagreen',label=r"$F_{\langle\tilde{q_t'}\rangle}$")
axs[1].plot(time[plttime],-thlvpfi_prod_dry,c='maroon',label=r"$-\tilde{w'}\frac{\partial \overline{\theta_{lv}}}{\partial z}$")
axs[1].plot(time[plttime],-thlvpfi_vdiv_dry,c='peachpuff',label=r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'\theta_{lv}'}-\overline{w'\theta_{lv}'}\right)\right)$")
axs[1].plot(time[plttime],-thlvpfi_hdiv_dry,c='olive',label=r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'\theta_{lv}'}\right)$")
axs[1].plot(time[plttime],-thlvpfi_subs_dry,c='skyblue',label=r"$-\overline{w_{LS}}\frac{\partial \tilde{\theta_{lv}'}}{\partial z}$")
axs[1].plot(time[plttime],thlvpfi_resid_dry,c='slategray',label=r"Residual")
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry region')

axs[0].set_ylabel('Large-scale heating rate [K/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))