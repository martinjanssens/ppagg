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

# Run specifics
lp = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/BOMEXStability/bomex200_e12/ppagg_ep'
ds1= nc.Dataset(lp+'/../profiles.001.nc')
ilp = np.loadtxt(lp+'/../lscale.inp.001')

time = np.load(lp+'/time.npy')

time1d = np.ma.getdata(ds1.variables['time'][:])
rhobf = np.ma.getdata(ds1.variables['rhobf'][:])

# Larger-scale processes
zf = ilp[:,0]
wfls = ilp[:,3]
qtavp_ls = ilp[:,6]
thlavp_ls = ilp[:,7]

plttime = np.load(lp+'/plttime.npy')
zflim = np.load(lp+'/zf.npy')

dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain

izmin = np.where(zflim[0] == zf)[0][0]
izmax = np.where(zflim[-1] == zf)[0][0]+1

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
qtpf_diff_moist_time = np.load(lp+'/qtpf_diff_moist_time.npy')
qtpf_diff_dry_time = np.load(lp+'/qtpf_diff_dry_time.npy')

thlvpf_moist_time = np.load(lp+'/thlvpf_moist_time.npy')
thlvpf_dry_time = np.load(lp+'/thlvpf_dry_time.npy')
thlvpf_prod_moist_time = np.load(lp+'/thlvpf_prod_moist_time.npy')
thlvpf_prod_dry_time = np.load(lp+'/thlvpf_prod_dry_time.npy')
thlvpf_vdiv_moist_time = np.load(lp+'/thlvpf_vdiv_moist_time.npy')
thlvpf_vdiv_dry_time = np.load(lp+'/thlvpf_vdiv_dry_time.npy')
thlvpf_hdiv_moist_time = np.load(lp+'/thlvpf_hdiv_moist_time.npy')
thlvpf_hdiv_dry_time = np.load(lp+'/thlvpf_hdiv_dry_time.npy')
thlvpf_subs_moist_time = np.load(lp+'/thlvpf_subs_moist_time.npy')
thlvpf_subs_dry_time = np.load(lp+'/thlvpf_subs_moist_time.npy')
thlvpf_diff_moist_time = np.load(lp+'/thlvpf_diff_dry_time.npy')
thlvpf_diff_dry_time = np.load(lp+'/thlvpf_diff_dry_time.npy')

thlvpp_moist_time = np.load(lp+'/thlvpp_moist_time.npy')
thlvpp_dry_time = np.load(lp+'/thlvpp_dry_time.npy')
thlvpp_prod_moist_time = np.load(lp+'/thlvpp_prod_moist_time.npy')
thlvpp_prod_dry_time = np.load(lp+'/thlvpp_prod_dry_time.npy')
thlvpp_vdiv_moist_time = np.load(lp+'/thlvpp_vdiv_moist_time.npy')
thlvpp_vdiv_dry_time = np.load(lp+'/thlvpp_vdiv_dry_time.npy')
thlvpp_hdiv_moist_time = np.load(lp+'/thlvpp_hdiv_moist_time.npy')
thlvpp_hdiv_dry_time = np.load(lp+'/thlvpp_hdiv_dry_time.npy')
thlvpp_subs_moist_time = np.load(lp+'/thlvpp_subs_moist_time.npy')
thlvpp_subs_dry_time = np.load(lp+'/thlvpp_subs_dry_time.npy')
thlvpp_diff_moist_time = np.load(lp+'/thlvpp_diff_moist_time.npy')
thlvpp_diff_dry_time = np.load(lp+'/thlvpp_diff_dry_time.npy')

wthlvpf_prod_moist_time = np.load(lp+'/wthlvpf_prod_moist_time.npy')
wthlvpf_prod_dry_time =  np.load(lp+'/wthlvpf_prod_dry_time.npy')
wthlvpf_vdiv_moist_time =  np.load(lp+'/wthlvpf_vdiv_moist_time.npy')
wthlvpf_vdiv_dry_time = np.load(lp+'/wthlvpf_vdiv_dry_time.npy')
wthlvpf_hdiv_moist_time = np.load(lp+'/wthlvpf_hdiv_moist_time.npy')
wthlvpf_hdiv_dry_time = np.load(lp+'/wthlvpf_hdiv_dry_time.npy')
wthlvpf_buoy_moist_time = np.load(lp+'/wthlvpf_buoy_moist_time.npy')
wthlvpf_buoy_dry_time = np.load(lp+'/wthlvpf_buoy_dry_time.npy')
wthlvpf_pres_moist_time = np.load(lp+'/wthlvpf_pres_moist_time.npy')
wthlvpf_pres_dry_time = np.load(lp+'/wthlvpf_pres_dry_time.npy')
wthlvpf_subs_moist_time = np.load(lp+'/wthlvpf_subs_moist_time.npy')
wthlvpf_subs_dry_time = np.load(lp+'/wthlvpf_subs_dry_time.npy')
wthlvpf_diff_moist_time = np.load(lp+'/wthlvpf_diff_moist_time.npy')
wthlvpf_diff_dry_time = np.load(lp+'/wthlvpf_diff_dry_time.npy')

thl_av_time = np.load(lp+'/thl_av_time.npy')
thlv_av_time = np.load(lp+'/thlv_av_time.npy')
qt_av_time = np.load(lp+'/qt_av_time.npy')

thlpf_moist_time = np.load(lp+'/thlpf_moist_time.npy')
thlpf_dry_time = np.load(lp+'/thlpf_dry_time.npy')
wff_moist_time = np.load(lp+'/wff_moist_time.npy')
wff_dry_time = np.load(lp+'/wff_dry_time.npy')
qlpf_moist_time = np.load(lp+'/qlpf_moist_time.npy') 
qlpf_dry_time = np.load(lp+'/qlpf_dry_time.npy')

thlpp_moist_time = np.load(lp+'/thlpp_moist_time.npy')
thlpp_dry_time = np.load(lp+'/thlpp_dry_time.npy')
wfp_moist_time = np.load(lp+'/wfp_moist_time.npy')
wfp_dry_time = np.load(lp+'/wfp_dry_time.npy')
qlpp_moist_time = np.load(lp+'/qlpp_moist_time.npy') 
qlpp_dry_time = np.load(lp+'/qlpp_dry_time.npy')

wthlpf_moist_time = np.load(lp+'/wthlpf_moist_time.npy')
wthlpf_dry_time = np.load(lp+'/wthlpf_dry_time.npy')

wqtpf_moist_time = np.load(lp+'/wqtpf_moist_time.npy')
wqtpf_dry_time = np.load(lp+'/wqtpf_dry_time.npy')

wqlpf_moist_time = np.load(lp+'/wqlpf_moist_time.npy')
wqlpf_dry_time = np.load(lp+'/wqlpf_dry_time.npy')

wthlvp_av_time = np.load(lp+'/wthlvp_av_time.npy')
wthlvpf_moist_time = np.load(lp+'/wthlvpf_moist_time.npy')
wthlvpf_dry_time = np.load(lp+'/wthlvpf_dry_time.npy')
wthlvpf_l_moist_time = np.load(lp+'/wthlvpf_l_moist_time.npy')
wthlvpf_l_dry_time = np.load(lp+'/wthlvpf_l_dry_time.npy')
wthlvpf_c_moist_time = np.load(lp+'/wthlvpf_c_moist_time.npy')
wthlvpf_c_dry_time = np.load(lp+'/wthlvpf_c_dry_time.npy')
wthlvpf_r_moist_time = np.load(lp+'/wthlvpf_r_moist_time.npy')
wthlvpf_r_dry_time = np.load(lp+'/wthlvpf_r_dry_time.npy')
wthlvpp_moist_time = np.load(lp+'/wthlvpp_moist_time.npy')
wthlvpp_dry_time = np.load(lp+'/wthlvpp_dry_time.npy')

thvpf_moist_time = thlvpf_moist_time + 7*thl_av_time*qlpf_moist_time
thvpf_dry_time = thlvpf_dry_time + 7*thl_av_time*qlpf_dry_time

# Reconstruct mean flux divergence (approximately, i.e. ignoring rho)
wthl_av = ds1['wthlt'][:,izmin:izmax]
wqt_av = ds1['wqtt'][:,izmin:izmax]
thl_av_1d = ds1['thl'][:,izmin:izmax]
wthlv_av = wthl_av + 0.608*thl_av_1d*wqt_av

ddz_wqt_av_time = -((wqt_av[:,1:] - wqt_av[:,:-1])/dzh)
ddz_wthlv_av_time = -((wthlv_av[:,1:] - wthlv_av[:,:-1])/dzh)

ddz_wqt_av_time = (ddz_wqt_av_time[:,1:] + ddz_wqt_av_time[:,:-1])*0.5
ddz_wthlv_av_time = (ddz_wthlv_av_time[:,1:] + ddz_wthlv_av_time[:,:-1])*0.5

# Slopes of mean profiles
Gamma_thlv = thlvpf_prod_moist_time/wff_moist_time[:,1:-1]
Gamma_qt = qtpf_prod_moist_wex_time/wff_moist_time[:,1:-1]

# Tendencies of variables of interest
def tderive(var,time):
    return ((var[1:,1:-1] - var[:-1,1:-1])
           /(time[1:,np.newaxis] - time[:-1,np.newaxis])/3600)

qtpf_tend_moist_time = np.zeros(qtpf_prod_moist_time.shape)
qtpf_tend_dry_time = np.zeros(qtpf_prod_moist_time.shape)
qtpf_tend_moist_time[1:,:] = tderive(qtpf_moist_time, time)
qtpf_tend_dry_time[1:,:] = tderive(qtpf_dry_time, time)

thlvpf_tend_moist_time = np.zeros(thlvpf_prod_moist_time.shape)
thlvpf_tend_dry_time = np.zeros(thlvpf_prod_moist_time.shape)
thlvpf_tend_moist_time[1:,:] = tderive(thlvpf_moist_time, time)
thlvpf_tend_dry_time[1:,:] = tderive(thlvpf_dry_time, time)

#%% Plotprofiles of  mesoscale-filtered variables in time
tpltmin = 6.
tpltmax = 23.
dit = 1.0 # Rounds to closest multiple of dt in time

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig,axs = plt.subplots(ncols=6,sharey=True,figsize=(12,5))
for i in range(len(plttime_var)):
    col = plt.cm.cubehelix(i/len(plttime_var))
     
    axs[0].plot(qtpf_moist_time[plttime_var[i],:], zflim, color=col,linestyle='-')
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].set_xlabel(r"$\widetilde{q_t'}$")
    axs[0].set_xlim((0,6e-4))
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[1].plot(thlvpf_moist_time[plttime_var[i],:], zflim, color=col,linestyle='-')
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].set_xlabel(r"$\widetilde{\theta_{lv}'}$")
    axs[1].set_xlim((-4e-2,1e-2))
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[2].plot(wff_moist_time[plttime_var[i],:], zflim,color=col,linestyle='-')
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].set_xlabel(r"$\widetilde{w'}$")
    axs[2].set_xlim((-1e-2,1.7e-2))
    axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[3].plot(thlpf_moist_time[plttime_var[i],:], zflim,color=col,linestyle='-')
    axs[3].axvline(0,color='gray',linestyle='dotted')
    axs[3].set_xlabel(r"$\widetilde{\theta_l'}$")
    axs[3].set_xlim((-1.2e-1,0))
    axs[3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[4].plot(thvpf_moist_time[plttime_var[i],:], zflim,color=col,linestyle='-')
    axs[4].axvline(0,color='gray',linestyle='dotted')
    axs[4].set_xlabel(r"$\widetilde{\theta_v'}$")
    axs[4].set_xlim((-2.5e-2,2.5e-2))
    axs[4].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[5].plot(qlpf_moist_time[plttime_var[i],:], zflim, label='t=%.2f'%time[plttime_var[i]],color=col,linestyle='-')
    axs[5].axvline(0,color='gray',linestyle='dotted')
    axs[5].set_xlabel(r"$\widetilde{q_l'}$")
    axs[5].set_xlim((0,9e-6))
    axs[5].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[5].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

#%% Plot profiles of small-scale-filtered variables in time
tpltmin = 12.
tpltmax = 13.
dit = 1.0 # Rounds to closest multiple of dt in time

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig,axs = plt.subplots(ncols=4,sharey=True,figsize=(8,5))
for i in range(len(plttime_var)):
    col = plt.cm.cubehelix(i/len(plttime_var))

    axs[0].plot(thlvpp_moist_time[plttime_var[i],:], zflim, color=col,linestyle='-')
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].set_xlabel(r"$\theta_{lv}'''$")
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[1].plot(wfp_moist_time[plttime_var[i],:], zflim,color=col,linestyle='-')
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].set_xlabel(r"$w'''$")
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[2].plot(thlpf_moist_time[plttime_var[i],:], zflim,color=col,linestyle='-')
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].set_xlabel(r"$\theta_l'''$")
    axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[3].plot(qlpp_moist_time[plttime_var[i],:], zflim, label='t=%.2f'%time[plttime_var[i]],color=col,linestyle='-')
    axs[3].axvline(0,color='gray',linestyle='dotted')
    axs[3].set_xlabel(r"$q_l'''$")
    axs[3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[3].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

#%%
# Average budget contributions over time dimension
tpltmin = 6.
tpltmax = 16.

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

qtpfmn_tend_moist = np.mean(qtpf_tend_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_moist_wex = np.mean(qtpf_prod_moist_wex_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_vdiv_moist = np.mean(qtpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_hdiv_moist = np.mean(qtpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_subs_moist = np.mean(qtpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_diff_moist = np.mean(qtpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_budg_moist = (-qtpfmn_prod_moist_wex[1:-1] - qtpfmn_vdiv_moist[1:-1]
                     -qtpfmn_hdiv_moist[1:-1] - qtpfmn_subs_moist[1:-1]
                     +qtpfmn_diff_moist)
qtpfmn_resi_moist = qtpfmn_tend_moist[1:-1] - qtpfmn_budg_moist

qtpfmn_tend_dry = np.mean(qtpf_tend_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_wex = np.mean(qtpf_prod_dry_wex_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_vdiv_dry = np.mean(qtpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_hdiv_dry = np.mean(qtpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_subs_dry = np.mean(qtpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_diff_dry = np.mean(qtpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_budg_dry = (-qtpfmn_prod_dry_wex[1:-1] - qtpfmn_vdiv_dry[1:-1]
                     -qtpfmn_hdiv_dry[1:-1] - qtpfmn_subs_dry[1:-1]
                     +qtpfmn_diff_dry)
qtpfmn_resi_dry = qtpfmn_tend_dry[1:-1] - qtpfmn_budg_dry

thlvpfmn_tend_moist = np.mean(thlvpf_tend_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_prod_moist = np.mean(thlvpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_vdiv_moist = np.mean(thlvpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_hdiv_moist = np.mean(thlvpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_subs_moist = np.mean(thlvpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_diff_moist = np.mean(thlvpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_budg_moist = (-thlvpfmn_prod_moist[1:-1] - thlvpfmn_vdiv_moist[1:-1]
                       -thlvpfmn_hdiv_moist[1:-1] - thlvpfmn_subs_moist[1:-1]
                       +thlvpfmn_diff_moist)
thlvpfmn_resi_moist = thlvpfmn_tend_moist[1:-1] - thlvpfmn_budg_moist

thlvpfmn_tend_dry = np.mean(thlvpf_tend_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_prod_dry = np.mean(thlvpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_vdiv_dry = np.mean(thlvpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_hdiv_dry = np.mean(thlvpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_subs_dry = np.mean(thlvpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_diff_dry = np.mean(thlvpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_budg_dry = (-thlvpfmn_prod_dry[1:-1] - thlvpfmn_vdiv_dry[1:-1]
                     -thlvpfmn_hdiv_dry[1:-1] - thlvpfmn_subs_dry[1:-1]
                     +thlvpfmn_diff_dry)
thlvpfmn_resi_dry = thlvpfmn_tend_dry[1:-1] - thlvpfmn_budg_dry

thlvppmn_prod_moist = np.mean(thlvpp_prod_moist_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_vdiv_moist = np.mean(thlvpp_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_hdiv_moist = np.mean(thlvpp_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_subs_moist = np.mean(thlvpp_subs_moist_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_diff_moist = np.mean(thlvpp_diff_moist_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_budg_moist = (-thlvppmn_prod_moist[1:-1] - thlvppmn_vdiv_moist[1:-1]
                       -thlvppmn_hdiv_moist[1:-1] - thlvppmn_subs_moist[1:-1]
                       +thlvppmn_diff_moist)
thlvppmn_prod_dry = np.mean(thlvpp_prod_dry_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_vdiv_dry = np.mean(thlvpp_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_hdiv_dry = np.mean(thlvpp_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_subs_dry = np.mean(thlvpp_subs_dry_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_diff_dry = np.mean(thlvpp_diff_dry_time[itpltmin:itpltmax,:],axis=0)
thlvppmn_budg_dry = (-thlvppmn_prod_dry[1:-1] - thlvppmn_vdiv_dry[1:-1]
                     -thlvppmn_hdiv_dry[1:-1] - thlvppmn_subs_dry[1:-1]
                     +thlvppmn_diff_dry)

wthlvpfmn_prod_moist = np.mean(wthlvpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_vdiv_moist = np.mean(wthlvpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_hdiv_moist = np.mean(wthlvpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_buoy_moist = np.mean(wthlvpf_buoy_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_pres_moist = np.mean(wthlvpf_pres_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_subs_moist = np.mean(wthlvpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_diff_moist = np.mean(wthlvpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_budg_moist = (-wthlvpfmn_prod_moist[1:-2] - wthlvpfmn_vdiv_moist[1:-1]
                        -wthlvpfmn_hdiv_moist[1:-2] - wthlvpfmn_subs_moist[1:-2]
                        +wthlvpfmn_buoy_moist[1:-2] - wthlvpfmn_pres_moist[1:-2]
                        +wthlvpfmn_diff_moist)

wthlvpfmn_prod_dry = np.mean(wthlvpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_vdiv_dry = np.mean(wthlvpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_hdiv_dry = np.mean(wthlvpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_buoy_dry = np.mean(wthlvpf_buoy_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_pres_dry = np.mean(wthlvpf_pres_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_subs_dry = np.mean(wthlvpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_diff_dry = np.mean(wthlvpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_budg_dry = (-wthlvpfmn_prod_dry[1:-2] - wthlvpfmn_vdiv_dry[1:-1]
                      -wthlvpfmn_hdiv_dry[1:-2] - wthlvpfmn_subs_dry[1:-2]
                      +wthlvpfmn_buoy_dry[1:-2] - wthlvpfmn_pres_dry[1:-2]
                      +wthlvpfmn_diff_dry)

# Budget terms
terms = [r"$\frac{\partial\langle\tilde{q_t'}\rangle}{\partial t}$",
         r"$-\tilde{w'}\frac{\partial \overline{q_t}}{\partial z}$",
         r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'''q_t'''}-\overline{w'q_t'}\right)\right)$",
         r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'q_t'}\right)$",
         r"$-\overline{w_{LS}}\frac{\partial \tilde{q_t'}}{\partial z}$",
         r"$\widetilde{\frac{\partial}{\partial x_j}\left(K_h\frac{\partial q_t'}{\partial x_j}\right)}+\widetilde{\frac{\partial}{\partial x_j}\left(K_h'\frac{\partial \overline{q_t}}{\partial x_j}\right)}$"
         ]

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(qtpfmn_tend_moist, zflim[1:-1],c='midnightblue')
axs[0].plot(-qtpfmn_prod_moist_wex, zflim[1:-1],c='darkseagreen')
axs[0].plot(-qtpfmn_vdiv_moist, zflim[1:-1],c='maroon')
axs[0].plot(-qtpfmn_hdiv_moist, zflim[1:-1],c='peachpuff')
axs[0].plot(-qtpfmn_subs_moist, zflim[1:-1],c='olive')
axs[0].plot(qtpfmn_diff_moist, zflim[2:-2],c='skyblue')
axs[0].plot(qtpfmn_resi_moist, zflim[2:-2],c='gray')
axs[0].set_xlabel(r"Contribution to $\frac{\partial\tilde{q_t'}}{\partial t}$")

axs[1].plot(qtpfmn_tend_dry, zflim[1:-1],c='midnightblue',label=terms[0])
axs[1].plot(-qtpfmn_prod_dry_wex, zflim[1:-1],c='darkseagreen',label=terms[1])
axs[1].plot(-qtpfmn_vdiv_dry, zflim[1:-1],c='maroon',label=terms[2])
axs[1].plot(-qtpfmn_hdiv_dry, zflim[1:-1],c='peachpuff',label=terms[3])
axs[1].plot(-qtpfmn_subs_dry, zflim[1:-1],c='olive',label=terms[4])
axs[1].plot(qtpfmn_diff_dry, zflim[2:-2],c='skyblue',label=terms[5])
axs[1].plot(qtpfmn_resi_dry, zflim[2:-2],c='gray',label='Residual')
axs[1].set_xlabel(r"Contribution to $\frac{\partial\tilde{q_t'}}{\partial t}$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))


fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(thlvpfmn_tend_moist, zflim[1:-1],c='midnightblue')
axs[0].plot(-thlvpfmn_prod_moist, zflim[1:-1],c='darkseagreen')
axs[0].plot(-thlvpfmn_vdiv_moist, zflim[1:-1],c='maroon')
axs[0].plot(-thlvpfmn_hdiv_moist, zflim[1:-1],c='peachpuff')
axs[0].plot(-thlvpfmn_subs_moist, zflim[1:-1],c='olive')
axs[0].plot( thlvpfmn_diff_moist, zflim[2:-2],c='skyblue')
axs[0].plot( thlvpfmn_resi_moist, zflim[2:-2],c='gray')
axs[0].set_xlabel(r"Contribution to $\frac{\partial\tilde{\theta_{lv}'}}{\partial t}$")

axs[1].plot(thlvpfmn_tend_dry, zflim[1:-1],c='midnightblue',label='Tendency')
axs[1].plot(-thlvpfmn_prod_dry, zflim[1:-1],c='darkseagreen',label='Gradient production')
axs[1].plot(-thlvpfmn_vdiv_dry, zflim[1:-1],c='maroon',label='Anomalous vertical flux divergence')
axs[1].plot(-thlvpfmn_hdiv_dry, zflim[1:-1],c='peachpuff',label='Horizontal divergence')
axs[1].plot(-thlvpfmn_subs_dry, zflim[1:-1],c='olive',label='Subsidence')
axs[1].plot (thlvpfmn_diff_dry, zflim[2:-2],c='skyblue',label='SFS diffusion')
axs[1].plot( thlvpfmn_resi_dry, zflim[2:-2],c='gray',label='Residual')
axs[1].set_xlabel(r"Contribution to $\frac{\partial\tilde{\theta_{lv}'}}{\partial t}$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))


fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(thlvppmn_budg_moist, zflim[2:-2],c='midnightblue')
axs[0].plot(-thlvppmn_prod_moist, zflim[1:-1],c='darkseagreen')
axs[0].plot(-thlvppmn_vdiv_moist, zflim[1:-1],c='maroon')
axs[0].plot(-thlvppmn_hdiv_moist, zflim[1:-1],c='peachpuff')
axs[0].plot(-thlvppmn_subs_moist, zflim[1:-1],c='olive')
# axs[0].plot(thlvppmn_diff_moist[k,:], zflim[2:-2],c='skyblue')
axs[0].set_xlabel(r"Contribution to $\frac{\partial\theta_{lv}'''}{\partial t}$")

axs[1].plot(thlvppmn_budg_dry, zflim[2:-2],c='midnightblue',label='Tendency')
axs[1].plot(-thlvppmn_prod_dry, zflim[1:-1],c='darkseagreen',label='Gradient production')
axs[1].plot(-thlvppmn_vdiv_dry, zflim[1:-1],c='maroon',label='Anomalous vertical flux divergence')
axs[1].plot(-thlvppmn_hdiv_dry, zflim[1:-1],c='peachpuff',label='Horizontal divergence')
axs[1].plot(-thlvppmn_subs_dry, zflim[1:-1],c='olive',label='Subsidence')
# axs[1].plot(thlvppmn_diff_dry[k,:], zflim[2:-2],c='skyblue',label='SFS diffusion')
axs[1].set_xlabel(r"Contribution to $\frac{\partial\theta_{lv}'''}{\partial t}$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(wthlvpfmn_budg_moist, zflim[2:-3],c='midnightblue')
axs[0].plot(-wthlvpfmn_prod_moist, zflim[1:-1],c='darkseagreen')
axs[0].plot(-wthlvpfmn_vdiv_moist, zflim[1:-2],c='maroon')
axs[0].plot(-wthlvpfmn_hdiv_moist, zflim[1:-1],c='peachpuff')
axs[0].plot( wthlvpfmn_buoy_moist, zflim[1:-1],c='mediumvioletred')
axs[0].plot(-wthlvpfmn_pres_moist, zflim[1:-1],c='slategray')
axs[0].plot(-wthlvpfmn_subs_moist, zflim[1:-1],c='olive')
axs[0].plot( wthlvpfmn_diff_moist, zflim[2:-3],c='skyblue')
axs[0].set_xlabel(r"Contribution to $\frac{\partial}{\partial t}\left( \widetilde{w'''\theta_{lv}'''}-\overline{w'''\theta_{lv}'''} \right)$")

axs[1].plot(wthlvpfmn_budg_dry, zflim[2:-3],c='midnightblue',label='Tendency')
axs[1].plot(-wthlvpfmn_prod_dry, zflim[1:-1],c='darkseagreen',label='Gradient production')
axs[1].plot(-wthlvpfmn_vdiv_dry, zflim[1:-2],c='maroon',label='Anomalous vertical flux divergence')
axs[1].plot(-wthlvpfmn_hdiv_dry, zflim[1:-1],c='peachpuff',label='Horizontal divergence')
axs[1].plot( wthlvpfmn_buoy_dry, zflim[1:-1],c='mediumvioletred',label='Buoyancy')
axs[1].plot(-wthlvpfmn_pres_dry, zflim[1:-1],c='slategray',label='Pressure gradient')
axs[1].plot(-wthlvpfmn_subs_dry, zflim[1:-1],c='olive',label='Subsidence')
axs[1].plot( wthlvpfmn_diff_dry, zflim[2:-3],c='skyblue',label='SFS diffusion')
axs[1].set_xlabel(r"Contribution to $\frac{\partial}{\partial t}\left(\widetilde{w'''\theta_{lv}'''}-\overline{w'''\theta_{lv}'''}\right)$")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))
plt.show()

# WTG-based model of moisture variance production

# Vertical velocities

# Exact
wffmn_moist = np.mean(wff_moist_time[itpltmin:itpltmax],axis=0)
wffmn_dry = np.mean(wff_dry_time[itpltmin:itpltmax],axis=0)

# w model with actual thlvpf_vdiv
wff_moist_wtg = -thlvpf_vdiv_moist_time/Gamma_thlv
wff_dry_wtg = -thlvpf_vdiv_dry_time/Gamma_thlv

wffmn_moist_wtg = np.mean(wff_moist_wtg[itpltmin:itpltmax,:],axis=0)
wffmn_dry_wtg = np.mean(wff_dry_wtg[itpltmin:itpltmax,:],axis=0)

# w model with simple reliance on qtpf
wff_moist_mod = -qtpf_prod_moist_time/Gamma_qt
wff_dry_mod = -qtpf_prod_dry_time/Gamma_qt

wffmn_moist_mod = np.mean(wff_moist_mod[itpltmin:itpltmax,:],axis=0)
wffmn_dry_mod = np.mean(wff_dry_mod[itpltmin:itpltmax,:],axis=0)

# Moisture variance production
qtpfmn_prod_moist = np.mean(-qtpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry = np.mean(-qtpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)

qtpf_prod_moist_wtg_time = wff_moist_wtg*Gamma_qt
qtpf_prod_dry_wtg_time = wff_dry_wtg*Gamma_qt

qtpfmn_prod_moist_wtg = np.mean(qtpf_prod_moist_wtg_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_wtg = np.mean(qtpf_prod_dry_wtg_time[itpltmin:itpltmax,:],axis=0)

# w plot
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(wffmn_moist, zflim, c='midnightblue')
axs[0].plot(wffmn_moist_wtg, zflim[1:-1], c='darkseagreen')
# axs[0].plot(wffmn_moist_mod, zflim[1:-1], c='maroon')
axs[0].set_xlabel(r"$\widetilde{w'}$ [m/s], moist region")

axs[1].plot(wffmn_dry, zflim, c='midnightblue', label=r"Ground truth $\widetilde{w'}'$")
axs[1].plot(wffmn_dry_wtg, zflim[1:-1], c='darkseagreen', label=r"WTG model")
# axs[1].plot(wffmn_dry_mod, zflim[1:-1], c='maroon', label=r"Full loop model")
axs[1].set_xlabel(r"$\widetilde{w'}$ [m/s], dry region")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

# Moisture variance production plot
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(qtpfmn_prod_moist_wex, zflim[1:-1], c='darkseagreen',linestyle='-')
axs[0].plot(qtpfmn_prod_moist_wtg, zflim[1:-1], c='darkseagreen',linestyle='--')
# axs[0].plot(qtpfmn_prod_moist, zflim[1:-1], c='darkseagreen',linestyle='dotted')
axs[0].set_xlabel(r"$\widetilde{w'}\frac{\partial\overline{q_t}}{\partial z}$ [kg/kg/s], moist region")

axs[1].plot(qtpfmn_prod_dry_wex, zflim[1:-1], c='darkseagreen',linestyle='-', label=r"Ground truth")
axs[1].plot(qtpfmn_prod_dry_wtg, zflim[1:-1], c='darkseagreen',linestyle='--', label=r"WTG model")
# axs[1].plot(qtpfmn_prod_dry, zflim[1:-1], c='darkseagreen',linestyle='dotted', label=r"Full loop model")
axs[1].set_xlabel(r"$\widetilde{w'}\frac{\partial\overline{q_t}}{\partial z}$ [kg/kg/s], dry region")

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))



#%% Vertically integrated statistics
tpltmin = 6.
tpltmax = 18.
dit = 0.5 # Rounds to closest multiple of dt in time

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

def vint(field,rhob,z,plttime=plttime):
    
    if len(field.shape) == 3:
        var = np.trapz(rhob[:,np.newaxis,np.newaxis]*field[:,:,:],z,axis=0)
    elif len(field.shape) == 4:
        var = np.trapz(rhob[np.newaxis,:,np.newaxis,np.newaxis]*
                       field[plttime,:,:,:],z,axis=1)
    elif len(field.shape) == 2:
        var = np.trapz(rhob[np.newaxis,:]*field[plttime,:],z,axis=1)
    elif len(field.shape) == 1:
        var = np.trapz(rhob*field,z)
    return var
   
# 1D fields
rhobfi = rhobf[0,izmin:izmax] # Won't really change much through time, so ok to take 0 value

qtpfi_moist = vint(qtpf_moist_time,rhobfi,zflim,plttime_var)
qtpfi_dry = vint(qtpf_dry_time,rhobfi,zflim,plttime_var)

# Moistening gradient production per simplified WTG budget
qtpfi_prod_moist = vint(qtpf_prod_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_prod_dry = vint(qtpf_prod_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

qtpfi_prod_wex_moist = vint(qtpf_prod_moist_wex_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_prod_wex_dry = vint(qtpf_prod_dry_wex_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Moistening through anomalous vertical small-scale fluxes
# FIXME offset zf in integration by 1 from field
qtpfi_vdiv_moist = vint(qtpf_vdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_vdiv_dry = vint(qtpf_vdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Moistening through horizontal advection
qtpfi_hdiv_moist = vint(qtpf_hdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_hdiv_dry = vint(qtpf_hdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Moistening through subsidence
qtpfi_subs_moist = vint(qtpf_subs_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_subs_dry = vint(qtpf_subs_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Moistening through SFS diffusion
# qtpfi_diff_moist = vint(qtpf_diff_moist_time,rhobfi,zflim,plttime_var)
# qtpfi_diff_dry = vint(qtpf_diff_dry_time,rhobfi,zflim,plttime_var)

# Fit the moisture fluctuation's evolution over the intended time
[exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime_var]*3600, 
                                       qtpfi_moist,
                                       p0=[1,1e-5])

[exp_dry,fac_dry], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime_var]*3600, 
                                       qtpfi_dry,
                                       p0=[-1,-1e-5])

# And differentiate to estimate its tendency
qtpfi_tend_moist = fac_moist*exp_moist*(time[plttime_var]*3600)**(exp_moist-1)
qtpfi_tend_dry = fac_dry*exp_dry*(time[plttime_var]*3600)**(exp_dry-1)

# Estimate residual
qtpfi_resid_moist = qtpfi_tend_moist + qtpfi_prod_wex_moist + qtpfi_vdiv_moist + qtpfi_hdiv_moist + qtpfi_subs_moist #- qtpfi_diff_moist
qtpfi_resid_dry = qtpfi_tend_dry + qtpfi_prod_wex_dry + qtpfi_vdiv_dry + qtpfi_hdiv_dry + qtpfi_subs_dry #- qtpfi_diff_dry

# Temporal plot
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
axs[0].plot(time[plttime_var],qtpfi_tend_moist,c='midnightblue')
axs[0].plot(time[plttime_var],-qtpfi_prod_wex_moist,c='darkseagreen')
axs[0].plot(time[plttime_var],-qtpfi_vdiv_moist,c='maroon')
axs[0].plot(time[plttime_var],-qtpfi_hdiv_moist,c='peachpuff')
axs[0].plot(time[plttime_var],-qtpfi_subs_moist,c='olive')
# axs[0].plot(time[plttime_var],qtpfi_diff_moist,c='skyblue')
axs[0].plot(time[plttime_var],qtpfi_resid_moist,c='slategray')
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist region')

axs[1].plot(time[plttime_var],qtpfi_tend_dry,c='midnightblue',label=terms[0])
axs[1].plot(time[plttime_var],-qtpfi_prod_wex_dry,c='darkseagreen',label=terms[1])
axs[1].plot(time[plttime_var],-qtpfi_vdiv_dry,c='maroon',label=terms[2])
axs[1].plot(time[plttime_var],-qtpfi_hdiv_dry,c='peachpuff',label=terms[3])
axs[1].plot(time[plttime_var],-qtpfi_subs_dry,c='olive',label=terms[4])
# axs[1].plot(time[plttime_var],qtpfi_diff_dry,c='skyblue',label=terms[5])
axs[1].plot(time[plttime_var],qtpfi_resid_dry,c='slategray',label=r"Residual")
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry region')

axs[0].set_ylabel('Large-scale moistening rate [kg/kg/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

# Model evaluation
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
axs[0].plot(time[plttime_var],qtpfi_tend_moist,c='midnightblue')
axs[0].plot(time[plttime_var],-qtpfi_prod_wex_moist,c='darkseagreen')
axs[0].plot(time[plttime_var],qtpfi_prod_moist,c='darkseagreen',linestyle='dotted')
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist region')

axs[1].plot(time[plttime_var],qtpfi_tend_dry,c='midnightblue',label=terms[0])
axs[1].plot(time[plttime_var],-qtpfi_prod_wex_dry,c='darkseagreen',label=terms[1])
axs[1].plot(time[plttime_var],qtpfi_prod_dry,c='darkseagreen',linestyle='dotted',label='Modelled '+terms[1])
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry region')

axs[0].set_ylabel('Large-scale moistening rate [kg/kg/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

### And now for thlv
thlvpfi_moist = vint(thlvpf_moist_time,rhobfi,zflim,plttime_var)
thlvpfi_dry = vint(thlvpf_dry_time,rhobfi,zflim,plttime_var)

# thlv gradient production
thlvpfi_prod_moist = vint(thlvpf_prod_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
thlvpfi_prod_dry = vint(thlvpf_prod_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# heating through anomalous vertical fluxes
# FIXME offset zf in integration by 1 from field
thlvpfi_vdiv_moist = vint(thlvpf_vdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
thlvpfi_vdiv_dry = vint(thlvpf_vdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Heating through horizontal advection
thlvpfi_hdiv_moist = vint(thlvpf_hdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
thlvpfi_hdiv_dry = vint(thlvpf_hdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Heating through subsidence
thlvpfi_subs_moist = vint(thlvpf_subs_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
thlvpfi_subs_dry = vint(thlvpf_subs_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# Heating through SFS diffusion
# thlvpfi_diff_moist = vint(thlvpf_diff_moist_time,rhobfi,zflim,plttime_var)
# thlvpfi_diff_dry = vint(thlvpf_diff_dry_time,rhobfi,zflim,plttime_var)


# Fit the haeting fluctuation's evolution
[exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime_var]*3600, 
                                       thlvpfi_moist,
                                       p0=[1,0])

[exp_dry,fac_dry], cov = curve_fit(lambda x, a, b: b * x**a, 
                                       time[plttime_var]*3600, 
                                       thlvpfi_dry,
                                       p0=[-3,0])

# And differentiate to estimate its tendency
thlvpfi_tend_moist = fac_moist*exp_moist*(time[plttime_var]*3600)**(exp_moist-1)
thlvpfi_tend_dry = fac_dry*exp_dry*(time[plttime_var]*3600)**(exp_dry-1)

# Estimate residual
thlvpfi_resid_moist = thlvpfi_tend_moist + thlvpfi_prod_moist + thlvpfi_vdiv_moist + thlvpfi_hdiv_moist + thlvpfi_subs_moist
thlvpfi_resid_dry = thlvpfi_tend_dry + thlvpfi_prod_dry + thlvpfi_vdiv_dry + thlvpfi_hdiv_dry + thlvpfi_subs_dry

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
axs[0].plot(time[plttime_var],thlvpfi_tend_moist,c='midnightblue')
# axs[0].plot(time[plttime_var],qtpfi_prod_moist,c='darkseagreen')
axs[0].plot(time[plttime_var],-thlvpfi_prod_moist,c='maroon')
axs[0].plot(time[plttime_var],-thlvpfi_vdiv_moist,c='peachpuff')
axs[0].plot(time[plttime_var],-thlvpfi_hdiv_moist,c='olive')
axs[0].plot(time[plttime_var],-thlvpfi_subs_moist,c='skyblue')
axs[0].plot(time[plttime_var],thlvpfi_resid_moist,c='slategray')
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist region')

axs[1].plot(time[plttime_var],qtpfi_tend_dry,c='midnightblue',label=r"$\frac{\partial\langle\tilde{\theta_{lv}'}\rangle}{\partial t}$")
# axs[1].plot(time[plttime_var],qtpfi_prod_dry,c='darkseagreen',label=r"$F_{\langle\tilde{q_t'}\rangle}$")
axs[1].plot(time[plttime_var],-thlvpfi_prod_dry,c='maroon',label=r"$-\tilde{w'}\frac{\partial \overline{\theta_{lv}}}{\partial z}$")
axs[1].plot(time[plttime_var],-thlvpfi_vdiv_dry,c='peachpuff',label=r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'\theta_{lv}'}-\overline{w'\theta_{lv}'}\right)\right)$")
axs[1].plot(time[plttime_var],-thlvpfi_hdiv_dry,c='olive',label=r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'\theta_{lv}'}\right)$")
axs[1].plot(time[plttime_var],-thlvpfi_subs_dry,c='skyblue',label=r"$-\overline{w_{LS}}\frac{\partial \tilde{\theta_{lv}'}}{\partial z}$")
axs[1].plot(time[plttime_var],thlvpfi_resid_dry,c='slategray',label=r"Residual")
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry region')

axs[0].set_ylabel('Large-scale heating rate [K/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

#%% Fluxes (no sgs yet)

# Time to average over
tpltmin = 11.
tpltmax = 12.

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

thlmn_av = np.mean(thl_av_time[itpltmin:itpltmax,:],axis=0)

wthlpfmn_moist = np.mean(wthlpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthlpfmn_dry = np.mean(wthlpf_dry_time[itpltmin:itpltmax,:],axis=0)

wqtpfmn_moist = np.mean(wqtpf_moist_time[itpltmin:itpltmax,:],axis=0)
wqtpfmn_dry = np.mean(wqtpf_dry_time[itpltmin:itpltmax,:],axis=0)

wqlpfmn_moist = np.mean(wqlpf_moist_time[itpltmin:itpltmax,:],axis=0)
wqlpfmn_dry = np.mean(wqlpf_dry_time[itpltmin:itpltmax,:],axis=0)

wthlvpmn_av = np.mean(wthlvp_av_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_moist = np.mean(wthlvpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_dry = np.mean(wthlvpf_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_l_moist = np.mean(wthlvpf_l_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_l_dry = np.mean(wthlvpf_l_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_c_moist = np.mean(wthlvpf_c_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_c_dry = np.mean(wthlvpf_c_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_r_moist = np.mean(wthlvpf_r_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_r_dry = np.mean(wthlvpf_r_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvppmn_moist = np.mean(wthlvpp_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvppmn_dry = np.mean(wthlvpp_dry_time[itpltmin:itpltmax,:],axis=0)

# Plot with moist/dry for w_thlv (Total fluxes)
plt.plot(wthlvpfmn_moist,zflim,label=r"$\widetilde{w'''\theta_{lv}'''}$, moist")
plt.plot(wthlvpfmn_dry,zflim,label=r"$\widetilde{w'''\theta_{lv}'''}$, dry")
plt.plot(wthlvpmn_av,zflim,label=r"$\overline{w'\theta_{lv}'}$")
plt.ylabel('z [m]')
plt.xlabel(r"$w'\theta_{lv}'$ [Km/s]")
plt.legend()
plt.show()

# Plot with moist for wthl, wqt contribution to wthlv (full fluxes)
plt.plot(wthlvpfmn_moist,zflim,label=r"$\widetilde{w'''\theta_{lv}'''}$")
plt.plot(wthlpfmn_moist,zflim,label=r"$\widetilde{w'''\theta_{l}'''}$")
plt.plot(0.608*thlmn_av*wqtpfmn_moist,zflim,label=r"$0.608\overline{\theta_l}\widetilde{w'''q_t'''}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'''\theta_{lv}'''}$, moist region [Km/s]")
plt.legend()
plt.show()

plt.plot(wthlvpfmn_dry,zflim,label=r"$\widetilde{w'''\theta_{lv}'''}$")
plt.plot(wthlpfmn_dry,zflim,label=r"$\widetilde{w'''\theta_{l}'''}$")
plt.plot(0.608*thlmn_av*wqtpfmn_dry,zflim,label=r"$0.608\overline{\theta_l}\widetilde{w'''q_t'''}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'''\theta_{lv}'''}$, dry region [Km/s]")
plt.legend()
plt.show()

# Plot with moist anomalous wthlv # as a function of anomalous wthl and wqt
plt.plot(wthlvpfmn_moist-wthlvpmn_av,zflim,label=r"$\widetilde{w'''\theta_{lv}'''} - \overline{w'\theta_{lv}'}$")
# plt.plot(wthl_r_moist_tot-wthl_av[izmin+2:izmax+2],zflim,label=r"$\widetilde{w'''\theta_{l}'''} - \overline{w'\theta_{l}'}$")
# plt.plot(0.608*thl_av[izmin+2:izmax+2]*(wqt_r_moist_tot-wqt_av[izmin+2:izmax+2]),zflim,label=r"$\widetilde{w'''q_t'''} - \overline{w'q_t'}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'''\theta_{lv}'''} - \overline{w'\theta_{lv}'}$, moist region [Km/s]")
plt.legend()
plt.show()

# And in dry
plt.plot(wthlvpfmn_dry-wthlvpmn_av,zflim,label=r"$\widetilde{w'''\theta_{lv}'''} - \overline{w'\theta_{lv}'}$")
# plt.plot(wthl_r_dry_tot-wthl_av[izmin+2:izmax+2],zflim,label=r"$\widetilde{w'''\theta_{l}'''} - \overline{w'\theta_{l}'}$")
# plt.plot(0.608*thl_av[izmin+2:izmax+2]*(wqt_r_dry_tot-wqt_av[izmin+2:izmax+2]),zflim,label=r"$\widetilde{w'''q_t'''} - \overline{w'q_t'}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'''\theta_{lv}'''} - \overline{w'\theta_{lv}'}$, dry region [Km/s]")
plt.legend()
plt.show()

# Plot wthlv in moist/dry
plt.plot(wthlvppmn_moist,zflim,label=r"Moist")
plt.plot(wthlvppmn_dry,zflim,label=r"Dry")
plt.plot(wthlvpmn_av,zflim,label=r"Average")
plt.ylabel('z [m]')
plt.xlabel(r"$(w'''\theta_{lv}''')''$ [Km/s]")
plt.legend()
plt.show()

# Plot wthlv anomaly in large/small and moist/dry
plt.plot(wthlvpfmn_moist-wthlvpmn_av,zflim,label=r"Moist, large")
plt.plot(wthlvpfmn_dry-wthlvpmn_av,zflim,label=r"Dry, large")
plt.plot(wthlvppmn_moist,zflim,linestyle='--',c='C0',label=r"Moist, small")
plt.plot(wthlvppmn_dry,zflim,linestyle='--',c='C1',label=r"Dry, small")
plt.ylabel('z [m]')
plt.xlabel(r"$w'''\theta_{lv}''' - \overline{w'''\theta_{lv}'''}$ [Km/s]")
plt.legend()
plt.show()

# Plot scale decomposition of wthlvf
plt.plot(wthlvpfmn_moist,zflim,label=r"$\widetilde{w'\theta_{lv}'}$")
plt.plot(wthlvpfmn_l_moist,zflim,label=r"$\widetilde{\widetilde{w'}\widetilde{\theta_{lv}'}}$")
plt.plot(wthlvpfmn_c_moist,zflim,label=r"$\widetilde{\widetilde{w'}\theta_{lv}'''} + \widetilde{w'''\widetilde{\theta_{lv}'}}$")
plt.plot(wthlvpfmn_r_moist,zflim,label=r"$\widetilde{w'''\theta_{lv}''''}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'\theta_{lv}'}$ [Km/s]")
plt.xlim((-0.045,0.015))
plt.legend()
plt.show()

plt.plot(wthlvpfmn_dry,zflim,label=r"$\widetilde{w'\theta_{lv}'}$")
plt.plot(wthlvpfmn_l_dry,zflim,label=r"$\widetilde{\widetilde{w'}\widetilde{\theta_{lv}'}}$")
plt.plot(wthlvpfmn_c_dry,zflim,label=r"$\widetilde{\widetilde{w'}\theta_{lv}'''} + \widetilde{w'''\widetilde{\theta_{lv}'}}$")
plt.plot(wthlvpfmn_r_dry,zflim,label=r"$\widetilde{w'''\theta_{lv}''''}$")
plt.ylabel('z [m]')
plt.xlabel(r"$\widetilde{w'\theta_{lv}'}$ [Km/s]")
plt.xlim((-0.045,0.015))
plt.legend()
plt.show()

#%% Flux in time

tpltmin = 6.
tpltmax = 18.
dit = 1.0 # Rounds to closest multiple of dt in time

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

# Average flux
fig,axs = plt.subplots(ncols=1,figsize=(5,5))
for i in range(len(plttime_var)):
    col = plt.cm.cubehelix(i/len(plttime_var))
    
    wthlvpmn_av = np.mean(wthlvp_av_time[plttime_var[i:i+2],:],axis=0)

    axs.plot(wthlvpmn_av, zflim, color=col,linestyle='-', label='t=%.2f'%time[plttime_var[i]])
    axs.axvline(0,color='gray',linestyle='dotted')
    axs.set_xlabel(r"$\overline{w'\theta_{lv}'}$ [Km/s]")
    axs.set_xlim((-2.5e-2,1e-2))
    axs.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs.set_ylabel('z [m]')
axs.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

# Flux anomaly
wthlvpf_moist_anom = wthlvpf_moist_time - wthlvp_av_time
wthlvpf_dry_anom = wthlvpf_dry_time - wthlvp_av_time

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
for i in range(len(plttime_var)):
    col = plt.cm.cubehelix(i/len(plttime_var))
    
    wthlvpfmn_moist_anom = np.mean(wthlvpf_moist_anom[plttime_var[i:i+2],:],axis=0)
    wthlvpfmn_dry_anom = np.mean(wthlvpf_dry_anom[plttime_var[i:i+2],:],axis=0)
     
    axs[0].plot(wthlvpfmn_moist_anom, zflim, color=col,linestyle='-')
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].set_xlabel(r"$\widetilde{w'''\theta_{lv}'''} - \overline{w'\theta_{lv}'}$, moist region [Km/s]")
    axs[0].set_xlim((-5e-2,0.5e-3))
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[1].plot(wthlvpfmn_dry_anom, zflim, color=col,linestyle='-', label='t=%.2f'%time[plttime_var[i]])
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].set_xlabel(r"$\widetilde{w'\theta_{lv}'} - \overline{w'\theta_{lv}'}$, dry region [Km/s]")
    axs[1].set_xlim((-0.5e-3,2.5e-2))
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

#%% Relation qtpf - wthlvpf_anom

tpltmin = 6.
tpltmax = 18.

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

# FIXME just copied from namoptions now
wthl0 = 8e-3
wqt0 = 5.2e-5
thl0 = 299.1
grav = 9.81

wstar = (grav/thl0*(wthl0+0.608*thl0*wqt0)*600)**(1/3)

# Flux anomaly
wthlvpf_moist_anom = wthlvpf_moist_time - wthlvp_av_time
wthlvpf_dry_anom = wthlvpf_dry_time - wthlvp_av_time

# qtpf
qtpf_moist_mod = 0.608*thl_av_time*qtpf_moist_time
qtpf_dry_mod = 0.608*thl_av_time*qtpf_dry_time

# Time filter
wthlvpf_moist_anom = wthlvpf_moist_anom[itpltmin:itpltmax,:]
wthlvpf_dry_anom = wthlvpf_dry_anom[itpltmin:itpltmax,:]
qtpf_moist_mod = qtpf_moist_mod[itpltmin:itpltmax,:]
qtpf_dry_mod = qtpf_dry_mod[itpltmin:itpltmax,:]

plt.scatter(qtpf_moist_mod.flatten(),wthlvpf_moist_anom.flatten(),c='C1',s=0.1)
plt.scatter(qtpf_dry_mod.flatten(),wthlvpf_dry_anom.flatten(),c='C1',s=0.1)
plt.ylabel(r"Average $\widetilde{w'\theta_{lv}'} - \overline{w'\theta_{lv}'}$ in moist and dry regions")
plt.xlabel(r"$0.608\overline{\theta_l}\widetilde{q_t'}$")
plt.legend()

#%% Mean profiles

tpltmin = 0.
tpltmax = 20.
dit = 4 # Rounds to closest multiple of dt in time
fac=1e5 # Factor for plotting tendency magnitudes
tav=0.25 # time (hrs) around plotting time over which to average flux divergence

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig = plt.figure(figsize=(5,5))
ax = plt.gca()
for i in range(len(plttime_var)):
    
    imin = np.argmin(wthlvp_av_time[plttime_var[i],:])
    ilab = r"$\widetilde{w'\theta_{lv}'}$ min" if i==0 else None
    # wthlv_min = np.min(wthlvp_av_time[plttime_var[i],:])
    # zmin_time = zf[imin]
    
    iqlbase = np.where(np.abs(qlpf_moist_time[plttime_var[i],:]) > 0)[0][0]
    qlbaselab = r"Cloud base" if i==0 else None
    
    iqltop = np.where(np.abs(qlpf_moist_time[plttime_var[i],:]) > 0)[0][-1]
    if iqltop > izmax-3:
        iqltop = izmax-3
    qltoplab = r"Cloud top" if i==0 else None
    
    # Subsidence at cloud base
    qtavp_subs = (-wfls[izmin+1:izmax-1]*Gamma_qt[plttime_var[i],:])
    thlvavp_subs = (-wfls[izmin+1:izmax-1]*Gamma_thlv[plttime_var[i],:])
    
    # Large-scale drying at cloud base
    qtavp_larq = (qtavp_ls[izmin:izmax])
    thlvavp_larq = (0.608*thl_av_time[plttime_var[i],:]*qtavp_ls[izmin:izmax])

    # Large-scale cooling at cloud base
    thlvavp_lart = (thlavp_ls[izmin:izmax])
    
    # Flux divergence, averaged over tav around the current time
    it1d = np.argmin(np.abs(time[plttime_var[i]]-time1d/3600))
    itmin = np.argmin(np.abs(time1d - (time1d[it1d]-tav*3600)))
    itmax = np.argmin(np.abs(time1d - (time1d[it1d]+tav*3600)))
    ddz_wthlv_av = np.mean(ddz_wthlv_av_time[itmin:itmax],axis=0)
    ddz_wqt_av = np.mean(ddz_wqt_av_time[itmin:itmax],axis=0)
    
    col = plt.cm.cubehelix(i/len(plttime_var))
    ax.plot(thlv_av_time[plttime_var[i],:], qt_av_time[plttime_var[i],:],
            label='t=%.2f'%time[plttime_var[i]],color=col)
    ax.scatter(thlv_av_time[plttime_var[i],imin],qt_av_time[plttime_var[i],imin],
               color=col,zorder=100,label=ilab)
    ax.scatter(thlv_av_time[plttime_var[i],iqlbase],qt_av_time[plttime_var[i],iqlbase],
               marker='s',color=col,zorder=100,label=qlbaselab)
    ax.scatter(thlv_av_time[plttime_var[i],iqltop],qt_av_time[plttime_var[i],iqltop],
               marker='^',color=col,zorder=100,label=qltoplab)
    ax.set_xlabel(r"$\overline{\theta_{lv}} [K]$")
    ax.set_ylabel(r"$\overline{q_t}$ [kg/kg]")
    
    heights_budg = [iqlbase, imin, iqltop]
    heights_budg = np.arange(iqlbase, imin, 4)
    if i == len(plttime_var)-1:
        for k in range(len(heights_budg)):

            subslab = r"subsidence" if k==0 else None  
            moislab = r"large-scale moistening" if k==0 else None
            heatlab = r"large-scale heating" if k==0 else None
            fluxlab = r"Flux convergence" if k==0 else None

            ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_subs[heights_budg[k]]],
                    [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*qtavp_subs[heights_budg[k]]],
                    color='olive',label=subslab)
            
            ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_larq[heights_budg[k]]],
                    [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*qtavp_larq[heights_budg[k]]],
                    color='palevioletred',label=moislab)
            
            ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_lart[heights_budg[k]]],
                    [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]],
                    color='peru',label=heatlab)
            
            ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*ddz_wthlv_av[heights_budg[k]]],
                    [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*ddz_wqt_av[heights_budg[k]]],
                    color='maroon',label=fluxlab)

ax.legend(loc='upper right',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)
    
    
#%% Profiles of mean flux divergence

tpltmin = 6.
tpltmax = 16.
dit = 2.0 # Rounds to closest multiple of dt in time
tav = 1.0

itpltmin = np.where(time1d/3600>=tpltmin)[0][0]
itpltmax = np.where(time1d/3600<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig,axs = plt.subplots(ncols=4,sharey=True,figsize=(10,5),squeeze=False)
for i in range(len(plttime_var)):
    col = plt.cm.cubehelix(i/len(plttime_var))
    
    itmin = np.argmin(np.abs(time1d - (time1d[plttime_var[i]]-tav*3600)))
    itmax = np.argmin(np.abs(time1d - (time1d[plttime_var[i]]+tav*3600)))
    
    thlv_av = np.mean(thlv_av_time[itmin:itmax],axis=0)
    qt_av = np.mean(qt_av_time[itmin:itmax],axis=0)
    ddz_wthlv_av = np.mean(ddz_wthlv_av_time[itmin:itmax],axis=0)
    ddz_wqt_av = np.mean(ddz_wqt_av_time[itmin:itmax],axis=0)
     
    axs[0,0].plot(thlv_av, zflim, color=col,linestyle='-')
    axs[0,0].set_xlabel(r"$\overline{\theta_{lv}}$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[0,1].plot(qt_av, zflim, color=col,linestyle='-')
    axs[0,1].set_xlabel(r"$\overline{q_t}$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[0,2].plot(ddz_wthlv_av, zflim[1:-1], color=col,linestyle='-')
    axs[0,2].axvline(0,color='gray',linestyle='dotted')
    axs[0,2].set_xlabel(r"$\frac{\partial}{\partial z}\left(\overline{w'\theta_{lv}'}\right)$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[0,3].plot(ddz_wqt_av, zflim[1:-1], color=col,linestyle='-',label='t=%.2f'%(time1d[plttime_var[i]]/3600))
    axs[0,3].axvline(0,color='gray',linestyle='dotted')
    axs[0,3].set_xlabel(r"$\frac{\partial}{\partial z}\left(\overline{w'q_t'}\right)$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0,0].set_ylabel('z [m]')
axs[0,3].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)



