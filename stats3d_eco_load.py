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
from dataloader import DataLoaderDALES, DataLoaderMicroHH
from functions import vint

# Run specifics
lp = '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/BOMEXStability/bomex200_e12/ppagg_ql'
lp = '/Users/martinjanssens/Documents/Wageningen/EUREC4A/moisture_circulation/eurec4a_mean/ppagg_new'
sp = lp+'/../figs'
mod = 'dales'

if mod == 'dales':
    dl = DataLoaderDALES(lp+'/..')
elif mod == 'microhh':
    dl = DataLoaderMicroHH(lp+'/..')
    
time1d = dl.time1d
rhobf = dl.rhobf

# Larger-scale processes
zf_inp = dl.zf_inp
wfls = dl.wfls
dqdt_ls = dl.dqdt_ls
dthldt_ls = dl.dthldt_ls

time = np.load(lp+'/time.npy')
plttime = np.load(lp+'/plttime.npy')
zflim = np.load(lp+'/zf.npy')

dzh = np.diff(zflim)[0] # FIXME only valid in lower part of domain

izmin = np.where(zflim[0] >= zf_inp)[0][0]
izmax = np.where(zflim[-1] <= zf_inp)[0][0]+1

rhobfi = rhobf[0,izmin:izmax] # Won't really change much through time, so ok to take 0 value

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
qtpf_micr_moist_time = np.load(lp+'/qtpf_micr_moist_time.npy')
qtpf_micr_dry_time = np.load(lp+'/qtpf_micr_dry_time.npy')

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
thlvpf_diff_moist_time = np.load(lp+'/thlvpf_diff_moist_time.npy')
thlvpf_diff_dry_time = np.load(lp+'/thlvpf_diff_dry_time.npy')
thlvpf_radi_moist_time = np.load(lp+'/thlvpf_radi_moist_time.npy')
thlvpf_radi_dry_time = np.load(lp+'/thlvpf_radi_dry_time.npy')
thlvpf_micr_moist_time = np.load(lp+'/thlvpf_micr_moist_time.npy')
thlvpf_micr_dry_time = np.load(lp+'/thlvpf_micr_dry_time.npy')

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

qlpf_vdiv_moist_time = np.load(lp+'/qlpf_vdiv_moist_time.npy')
qlpf_vdiv_dry_time = np.load(lp+'/qlpf_vdiv_dry_time.npy')

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

wthlp_av_time = np.load(lp+'/wthlp_av_time.npy')
wthlpf_moist_time = np.load(lp+'/wthlpf_moist_time.npy')
wthlpf_dry_time = np.load(lp+'/wthlpf_dry_time.npy')

wqtp_av_time = np.load(lp+'/wqtp_av_time.npy')
wqtpf_moist_time = np.load(lp+'/wqtpf_moist_time.npy')
wqtpf_dry_time = np.load(lp+'/wqtpf_dry_time.npy')

wqlp_av_time = np.load(lp+'/wqlp_av_time.npy')
wqlpf_moist_time = np.load(lp+'/wqlpf_moist_time.npy')
wqlpf_dry_time = np.load(lp+'/wqlpf_dry_time.npy')
wqlpf_l_moist_time = np.load(lp+'/wqlpf_l_moist_time.npy')
wqlpf_l_dry_time = np.load(lp+'/wqlpf_l_dry_time.npy')
wqlpf_c_moist_time = np.load(lp+'/wqlpf_c_moist_time.npy')
wqlpf_c_dry_time = np.load(lp+'/wqlpf_c_dry_time.npy')
wqlpf_r_moist_time = np.load(lp+'/wqlpf_r_moist_time.npy')
wqlpf_r_dry_time = np.load(lp+'/wqlpf_r_dry_time.npy')

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

# Flux anomaly
wthlvpf_moist_anom = wthlvpf_moist_time - wthlvp_av_time
wthlvpf_dry_anom = wthlvpf_dry_time - wthlvp_av_time

# Buoyancy approximation
thvpf_moist_time = thlvpf_moist_time + 7*thl_av_time*qlpf_moist_time
thvpf_dry_time = thlvpf_dry_time + 7*thl_av_time*qlpf_dry_time

# Mean ql (we don't have this from stats3d)
ql_av_1d = dl.load_qlav(izmin, izmax)

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

wthlvpf_tend_moist_time = np.zeros(wthlvpf_prod_moist_time.shape)
wthlvpf_tend_dry_time = np.zeros(wthlvpf_prod_moist_time.shape)
wthlvpf_tend_moist_time[1:,:] = tderive(wthlvpf_moist_anom, time)
wthlvpf_tend_dry_time[1:,:] = tderive(wthlvpf_dry_anom, time)


## Reconstruct slab-mean budget terms

thl_av_1d = dl.load_thlav(izmin, izmax)
qt_av_1d = dl.load_qtav(izmin, izmax)
thlv_av_1d = thl_av_1d*(1 + 0.608*qt_av_1d)

# Tendencies
ddt_thlv_av_time = tderive(thlv_av_1d, time1d/3600)
ddt_qt_av_time = tderive(qt_av_1d, time1d/3600)

# Flux divergence (approximately, i.e. ignoring rho)
wthl_av = dl.load_wthlav(izmin, izmax)
wqt_av = dl.load_wqtav(izmin, izmax)
wthlv_av = wthl_av + 0.608*thl_av_1d*wqt_av

ddz_wthlv_av_time = ((wthlv_av[:,1:] - wthlv_av[:,:-1])/dzh)
ddz_wqt_av_time = ((wqt_av[:,1:] - wqt_av[:,:-1])/dzh)

ddz_wthlv_av_time = (ddz_wthlv_av_time[:,1:] + ddz_wthlv_av_time[:,:-1])*0.5
ddz_wqt_av_time = (ddz_wqt_av_time[:,1:] + ddz_wqt_av_time[:,:-1])*0.5

# Subsidence
Gamma_thlv_1d = (thlv_av_1d[:,1:] - thlv_av_1d[:,:-1])/dzh
Gamma_thlv_1d = (Gamma_thlv_1d[:,1:] + Gamma_thlv_1d[:,:-1])/2.

Gamma_qt_1d = (qt_av_1d[:,1:] - qt_av_1d[:,:-1])/dzh
Gamma_qt_1d = (Gamma_qt_1d[:,1:] + Gamma_qt_1d[:,:-1])/2.

wfls_dthlvdz_av_time = wfls[izmin+1:izmax-1]*Gamma_thlv_1d
wfls_dqtdz_av_time = wfls[izmin+1:izmax-1]*Gamma_qt_1d

# Large scale warming
dqdt_ls = dqdt_ls[izmin:izmax]
dthldt_ls = dthldt_ls[izmin:izmax]
dthlvdt_ls = dthldt_ls + 0.608*thl_av_1d*dqdt_ls

#%% Plotprofiles of  mesoscale-filtered variables in time
tpltmin = 2.
tpltmax = 14.
dit = 2.0 # Rounds to closest multiple of dt in time
dtav = 2.0 # Around each plotted time step
alpha = 0.5
lw=2

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
idtav  = int(round(dtav/2/(time[1]-time[0])))
idtav1 = int(round(dtav/2/(time1d[1]-time1d[0])*3600))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

def add_twinx(ax, norm, offset, label, return_axs=False):
    ax2 = ax.twiny()
    ax2.xaxis.set_ticks_position('top')
    ax2.xaxis.set_label_position('top')
    ax2.spines['top'].set_position(('axes',offset))
    ax2.set_xlim((ax.get_xlim()[0]/norm,ax.get_xlim()[1]/norm))
    ax2.set_xlabel(label)
    if return_axs:
        return ax2

ax2offs=1.0
fig,axs = plt.subplots(ncols=6,sharey=True,figsize=(14,5))
for i in range(len(plttime_var)):
    ti = time[plttime_var[i]]
    
    colm = plt.cm.Blues(i/len(plttime_var))
    cold = plt.cm.Reds(i/len(plttime_var))
    colc = plt.cm.Greys(i/len(plttime_var))
    
    it1d = np.argmin(abs(ti-time1d/3600))
    
    ql_avi = np.mean(ql_av_1d[it1d-idtav1:it1d+idtav1,:],axis=0)
    z_cb = zflim[ql_avi>0][1]
    
    z_ib = zflim[np.argmin(np.mean(wthlv_av[it1d-idtav1:it1d+idtav1],axis=0))]
    
    z_ct = zflim[ql_avi>0][-10]
    
    axs[0].plot(np.mean(qtpf_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim, 
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[0].plot(np.mean(qtpf_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[0].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[0].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[0].annotate('a)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[0].set_xlabel(r"$q_{t_m}'$ [kg/kg]")
        axs[0].set_xlim((-6e-4,6e-4))
        axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        qt2max = np.max(np.sqrt(dl.load_qt2av(izmin,izmax)[it1d]))
        add_twinx(axs[0], qt2max, ax2offs, r"$q_{t_m}'/\max \overline{\sqrt{q_t'^2}}$")

    axs[1].plot(np.mean(qlpf_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[1].plot(np.mean(qlpf_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[1].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[1].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[1].annotate('b)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[1].set_xlabel(r"$q_{l_m}'$ [kg/kg]")
        axs[1].set_xlim((-9e-6,9e-6))
        axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        add_twinx(axs[1], np.max(ql_av_1d), ax2offs, r"$q_{l_m}'/\max \overline{q_l}$")

    axs[2].plot(np.mean(wff_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[2].plot(np.mean(wff_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[2].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[2].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[2].annotate('c)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[2].set_xlabel(r"$w_m'$ [m/s]")
        axs[2].set_xlim((-1.7e-2,1.7e-2))
        axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        w2max = np.max(np.sqrt(dl.load_w2tav(izmin,izmax)[it1d]))
        add_twinx(axs[2], w2max, ax2offs, r"$w_m'/\max \overline{\sqrt{w^2}}$")

    axs[3].plot(np.mean(thlpf_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[3].plot(np.mean(thlpf_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[3].axvline(0,color='gray',linestyle='dotted')
    axs[3].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[3].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[3].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[3].annotate('d)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[3].set_xlabel(r"$\theta_{l_m}'$ [K]")
        axs[3].set_xlim((-1.2e-1,1.2e-1))
        axs[3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        thl2max = np.max(np.sqrt(dl.load_thl2av(izmin,izmax)[it1d]))
        axs32 = add_twinx(axs[3], thl2max, ax2offs, r"$\theta_{l_m}'/\max \overline{\sqrt{\theta_l'^2}}$", return_axs=True)
        # axs32.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[4].plot(np.mean(thvpf_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[4].plot(np.mean(thvpf_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[4].axvline(0,color='gray',linestyle='dotted')
    axs[4].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[4].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[4].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[4].annotate('e)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[4].set_xlabel(r"$\theta_{v_m}'$ [K]")
        axs[4].set_xlim((-2.6e-2,2.6e-2))
        axs[4].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        thv2max = np.max(np.sqrt(dl.load_thv2av(izmin,izmax)[it1d]))
        axs42 = add_twinx(axs[4], thv2max, ax2offs, r"$\theta_{v_m}'/\max \overline{\sqrt{\theta_v'^2}}$", return_axs=True)
        # axs42.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[5].plot(np.mean(thlvpf_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                label='%.2f'%ti,color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[5].plot(np.mean(thlvpf_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                label=' ',color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[5].axvline(0,color='gray',linestyle='dotted')
    axs[5].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[5].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[5].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    if i == len(plttime_var)-1:
        axs[5].annotate('f)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        axs[5].set_xlabel(r"$\theta_{lv_m}'$ [K]")
        axs[5].set_xlim((-4e-2,4e-2))
        axs[5].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
        axs52 = add_twinx(axs[5], thv2max, ax2offs, r"$\theta_{lv_m}'/\max \overline{\sqrt{\theta_{lv}'^2}}$", return_axs=True)
        # axs52.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
axs[5].annotate(r"Cloud base",(4.5e-2,z_cb),annotation_clip=False)
axs[5].annotate(r"Inversion base",(4.5e-2,z_ib),annotation_clip=False)
axs[5].annotate(r"Cloud top",(4.5e-2,z_ct),annotation_clip=False)
handles, labels = axs[5].get_legend_handles_labels()
handm = handles[::2];  labsm = labels[::2]
handd = handles[1::2]; labsd = labels[1::2]
handles = np.concatenate((handm,handd))
labels  = np.concatenate((labsm,labsd))
axs[5].legend(handles, labels, loc='best',bbox_to_anchor=(1.8,1),
              ncol=2,title='moist  time [hr]   dry')
plt.savefig(sp+'/vars_meso_evo.pdf', bbox_inches='tight')

#%% Plot profiles of small-scale-filtered variables in time
tpltmin = 6.
tpltmax = 24.
dit = 1.0 # Rounds to closest multiple of dt in time
dtav = 1.0 # Around each plotted time step
alpha = 0.5
lw=2

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
idtav  = int(round(dtav/2/(time[1]-time[0])))
idtav1 = int(round(dtav/2/(time1d[1]-time1d[0])*3600))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig,axs = plt.subplots(ncols=4,sharey=True,figsize=(12,5))
for i in range(len(plttime_var)):
    ti = time[plttime_var[i]]
    
    colm = plt.cm.Blues(i/len(plttime_var))
    cold = plt.cm.Reds(i/len(plttime_var))
    colc = plt.cm.Greys(i/len(plttime_var))
    
    it1d = np.argmin(abs(ti-time1d/3600))
    
    ql_avi = np.mean(ql_av_1d[it1d-idtav1:it1d+idtav1],axis=0)
    z_cb = zflim[ql_avi>0][1]
    
    z_ib = zflim[np.argmin(np.mean(wthlv_av[it1d-idtav1:it1d+idtav1],axis=0))]
    
    z_ct = zflim[ql_avi>0][-10]
    
    axs[0].plot(np.mean(qlpp_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[0].plot(np.mean(qlpp_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[0].axvline(0,color='gray',linestyle='dotted')
    axs[0].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[0].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[0].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    axs[0].annotate('b)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
    axs[0].set_xlabel(r"$q_{l_s}'$")
    # axs[0].set_xlim((-9e-6,9e-6))
    axs[0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[1].plot(np.mean(wfp_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[1].plot(np.mean(wfp_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[1].axvline(0,color='gray',linestyle='dotted')
    axs[1].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[1].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[1].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    axs[1].annotate('c)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
    axs[1].set_xlabel(r"$w_s'$")
    # axs[1].set_xlim((-1.7e-2,1.7e-2))
    axs[1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[2].plot(np.mean(thlpp_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[2].plot(np.mean(thlpp_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[2].axvline(0,color='gray',linestyle='dotted')
    axs[2].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[2].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[2].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    axs[2].annotate('d)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
    axs[2].set_xlabel(r"$\theta_{l_s}'$")
    # axs[2].set_xlim((-1.2e-1,1.2e-1))
    axs[2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[3].plot(np.mean(thlvpp_moist_time[plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                label='%.2f'%ti,color=colm,linestyle='-',alpha=alpha,lw=lw)
    axs[3].plot(np.mean(thlvpp_dry_time  [plttime_var[i]-idtav:plttime_var[i]+idtav,:],axis=0), zflim,
                label=' ',color=cold,linestyle='-',alpha=alpha,lw=lw)
    axs[3].axvline(0,color='gray',linestyle='dotted')
    axs[3].axhline(z_cb,color=colc,linestyle='-',alpha=alpha)
    axs[3].axhline(z_ib,color=colc,linestyle='-',alpha=alpha)
    axs[3].axhline(z_ct,color=colc,linestyle='-',alpha=alpha)
    axs[3].annotate('f)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
    axs[3].set_xlabel(r"$\theta_{lv_s}'$")
    # axs[3].set_xlim((-4e-2,4e-2))
    axs[3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0].set_ylabel('z [m]')
# axs[3].annotate(r"Cloud base",(4.5e-2,z_cb),annotation_clip=False)
# axs[3].annotate(r"Inversion base",(4.5e-2,z_ib),annotation_clip=False)
# axs[3].annotate(r"Cloud top",(4.5e-2,z_ct),annotation_clip=False)
handles, labels = axs[3].get_legend_handles_labels()
handm = handles[::2];  labsm = labels[::2]
handd = handles[1::2]; labsd = labels[1::2]
handles = np.concatenate((handm,handd))
labels  = np.concatenate((labsm,labsd))
axs[3].legend(handles, labels, loc='best',bbox_to_anchor=(1.8,1),
              ncol=2,title='moist  time [hr]   dry')
plt.savefig(sp+'/vars_small_evo.pdf', bbox_inches='tight')

#%% Average qtpf budget contributions over time dimension
tpltmin = 6.
tpltmax = 20.

# Budget terms
# terms = [r"$\frac{\partial\langle\tilde{q_t'}\rangle}{\partial t}$",
#          r"$-\tilde{w'}\frac{\partial \overline{q_t}}{\partial z}$",
#          r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'''q_t'''}-\overline{w'q_t'}\right)\right)$",
#          r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'q_t'}\right)$",
#          r"$-\overline{w_{LS}}\frac{\partial \tilde{q_t'}}{\partial z}$",
#          r"$\widetilde{\frac{\partial}{\partial x_j}\left(K_h\frac{\partial q_t'}{\partial x_j}\right)}+\widetilde{\frac{\partial}{\partial x_j}\left(K_h'\frac{\partial \overline{q_t}}{\partial x_j}\right)}$"
#          ]

terms = ['Tendency                               ',
         'Gradient production',
         'Vertical flux convergence',
         'Horizontal flux convergence',
         'Subsidence',
         'SFS diffusion',
         'Residual',
         'Precipitation'
         ]

colors = ['black',
          'cadetblue',
          'lightsteelblue',
          'olivedrab',
          'sienna',
          'goldenrod',
          'lightgray',
          'maroon',]

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

qtpfmn_tend_moist = np.mean(qtpf_tend_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_moist_wex = np.mean(qtpf_prod_moist_wex_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_vdiv_moist = np.mean(qtpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_hdiv_moist = np.mean(qtpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_subs_moist = np.mean(qtpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_diff_moist = np.mean(qtpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_micr_moist = np.mean(qtpf_micr_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_budg_moist = (-qtpfmn_prod_moist_wex[1:-1] - qtpfmn_vdiv_moist[1:-1]
                     -qtpfmn_hdiv_moist[1:-1] - qtpfmn_subs_moist[1:-1]
                     +qtpfmn_diff_moist + qtpfmn_micr_moist[2:-2])
qtpfmn_resi_moist = qtpfmn_tend_moist[1:-1] - qtpfmn_budg_moist

# The residual is mostly due to integration error of vertical transport
# -> Include the residual in this term
qtpfmn_tend_moist = qtpfmn_tend_moist[1:-1] - qtpfmn_resi_moist

qtpfmn_tend_dry = np.mean(qtpf_tend_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_wex = np.mean(qtpf_prod_dry_wex_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_vdiv_dry = np.mean(qtpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_hdiv_dry = np.mean(qtpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_subs_dry = np.mean(qtpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_diff_dry = np.mean(qtpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_micr_dry = np.mean(qtpf_micr_dry_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_budg_dry = (-qtpfmn_prod_dry_wex[1:-1] - qtpfmn_vdiv_dry[1:-1]
                     -qtpfmn_hdiv_dry[1:-1] - qtpfmn_subs_dry[1:-1]
                     +qtpfmn_diff_dry + qtpfmn_micr_dry[2:-2])
qtpfmn_resi_dry = qtpfmn_tend_dry[1:-1] - qtpfmn_budg_dry
qtpfmn_tend_dry = qtpfmn_tend_dry[1:-1] - qtpfmn_resi_dry


alpha = 0.75
lw = 2

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
# fig.suptitle(colors)
axs[0].plot(qtpfmn_tend_moist, zflim[2:-2],c=colors[0],alpha=alpha,lw=lw)
axs[0].plot(-qtpfmn_prod_moist_wex, zflim[1:-1],c=colors[1],alpha=alpha,lw=lw)
axs[0].plot(-qtpfmn_vdiv_moist, zflim[1:-1],c=colors[2],alpha=alpha,lw=lw)
axs[0].plot(-qtpfmn_hdiv_moist, zflim[1:-1],c=colors[3],alpha=alpha,lw=lw)
axs[0].plot(-qtpfmn_subs_moist, zflim[1:-1],c=colors[4],alpha=alpha,lw=lw)
axs[0].plot(qtpfmn_diff_moist, zflim[2:-2],c=colors[5],alpha=alpha,lw=lw)
axs[0].plot(qtpfmn_micr_moist, zflim,c=colors[7],alpha=alpha,lw=lw)
# axs[0].plot(qtpfmn_resi_moist, zflim[2:-2],c='gray')
axs[0].set_xlabel(r"Contribution to $q_{t_m}'$ tendency [kg/kg/s]")
axs[0].set_xlim((-7.5e-8,7.5e-8))
axs[0].set_title('Moist')

axs[1].plot(qtpfmn_tend_dry, zflim[2:-2],c=colors[0],label=terms[0],alpha=alpha,lw=lw)
axs[1].plot(-qtpfmn_prod_dry_wex, zflim[1:-1],c=colors[1],label=terms[1],alpha=alpha,lw=lw)
axs[1].plot(-qtpfmn_vdiv_dry, zflim[1:-1],c=colors[2],label=terms[2],alpha=alpha,lw=lw)
axs[1].plot(-qtpfmn_hdiv_dry, zflim[1:-1],c=colors[3],label=terms[3],alpha=alpha,lw=lw)
axs[1].plot(-qtpfmn_subs_dry, zflim[1:-1],c=colors[4],label=terms[4],alpha=alpha,lw=lw)
axs[1].plot(qtpfmn_diff_dry, zflim[2:-2],c=colors[5],label=terms[5],alpha=alpha,lw=lw)
axs[1].plot(qtpfmn_micr_dry, zflim,c=colors[7],label=terms[7],alpha=alpha,lw=lw)
# axs[1].plot(qtpfmn_resi_dry, zflim[2:-2],c='gray',label='Residual')
axs[1].set_xlabel(r"Contribution to $q_{t_m}'$ tendency [kg/kg/s]")
axs[1].set_xlim((-7.5e-8,7.5e-8))
axs[1].set_title('Dry')

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='upper left',bbox_to_anchor=(1,1))

plt.savefig(sp+'/qtpf_budget.pdf',bbox_inches='tight')

#%% Average thlvpf budget contributions over time dimension
tpltmin = 6.
tpltmax = 20.

terms = ['Tendency                               ',
         'Gradient production',
         'Vertical flux convergence',
         'Horizontal flux convergence',
         'Subsidence',
         'SFS diffusion',
         'Residual',
         'Precipitation',
         'Radiation'
         ]

colors = ['black',
          'cadetblue',
          'lightsteelblue',
          'olivedrab',
          'sienna',
          'goldenrod',
          'lightgray',
          'maroon',
          'midnightblue']

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

thlvpfmn_tend_moist = np.mean(thlvpf_tend_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_prod_moist = np.mean(thlvpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_vdiv_moist = np.mean(thlvpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_hdiv_moist = np.mean(thlvpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_subs_moist = np.mean(thlvpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_diff_moist = np.mean(thlvpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_micr_moist = np.mean(thlvpf_micr_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_radi_moist = np.mean(thlvpf_radi_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_budg_moist = (-thlvpfmn_prod_moist[1:-1] - thlvpfmn_vdiv_moist[1:-1]
                       -thlvpfmn_hdiv_moist[1:-1] - thlvpfmn_subs_moist[1:-1]
                       +thlvpfmn_diff_moist + thlvpfmn_micr_moist[1:-1]
                       +thlvpfmn_radi_moist[2:-2])
thlvpfmn_resi_moist = thlvpfmn_tend_moist[1:-1] - thlvpfmn_budg_moist
thlvpfmn_tend_moist = thlvpfmn_tend_moist[1:-1] - thlvpfmn_resi_moist

thlvpfmn_tend_dry = np.mean(thlvpf_tend_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_prod_dry = np.mean(thlvpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_vdiv_dry = np.mean(thlvpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_hdiv_dry = np.mean(thlvpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_subs_dry = np.mean(thlvpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_diff_dry = np.mean(thlvpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_micr_dry = np.mean(thlvpf_micr_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_radi_dry = np.mean(thlvpf_radi_dry_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_budg_dry = (-thlvpfmn_prod_dry[1:-1] - thlvpfmn_vdiv_dry[1:-1]
                     -thlvpfmn_hdiv_dry[1:-1] - thlvpfmn_subs_dry[1:-1]
                     +thlvpfmn_diff_dry + thlvpfmn_micr_dry[1:-1]
                     +thlvpfmn_radi_dry[2:-2])
thlvpfmn_resi_dry = thlvpfmn_tend_dry[1:-1] - thlvpfmn_budg_dry
thlvpfmn_tend_dry = thlvpfmn_tend_dry[1:-1] - thlvpfmn_resi_dry

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(thlvpfmn_tend_moist, zflim[2:-2],c=colors[0],alpha=alpha,lw=lw)
axs[0].plot(-thlvpfmn_prod_moist, zflim[1:-1],c=colors[1],alpha=alpha,lw=lw)
axs[0].plot(-thlvpfmn_vdiv_moist, zflim[1:-1],c=colors[2],alpha=alpha,lw=lw)
axs[0].plot(-thlvpfmn_hdiv_moist, zflim[1:-1],c=colors[3],alpha=alpha,lw=lw)
axs[0].plot(-thlvpfmn_subs_moist, zflim[1:-1],c=colors[4],alpha=alpha,lw=lw)
axs[0].plot( thlvpfmn_diff_moist, zflim[2:-2],c=colors[5],alpha=alpha,lw=lw)
axs[0].plot( thlvpfmn_micr_moist, zflim[1:-1],c=colors[7],alpha=alpha,lw=lw)
axs[0].plot( thlvpfmn_radi_moist, zflim,c=colors[8],alpha=alpha,lw=lw)
# axs[0].plot( thlvpfmn_resi_moist, zflim[2:-2],c='gray')
axs[0].set_xlabel(r"Contribution to $\theta_{lv_m}'$ tendency [K/s]")
axs[0].set_xlim((-5.5e-5,5.5e-5))
axs[0].set_title('Moist')

axs[1].plot(thlvpfmn_tend_dry, zflim[2:-2],c=colors[0],label=terms[0],alpha=alpha,lw=lw)
axs[1].plot(-thlvpfmn_prod_dry, zflim[1:-1],c=colors[1],label=terms[1],alpha=alpha,lw=lw)
axs[1].plot(-thlvpfmn_vdiv_dry, zflim[1:-1],c=colors[2],label=terms[2],alpha=alpha,lw=lw)
axs[1].plot(-thlvpfmn_hdiv_dry, zflim[1:-1],c=colors[3],label=terms[3],alpha=alpha,lw=lw)
axs[1].plot(-thlvpfmn_subs_dry, zflim[1:-1],c=colors[4],label=terms[4],alpha=alpha,lw=lw)
axs[1].plot (thlvpfmn_diff_dry, zflim[2:-2],c=colors[5],label=terms[5],alpha=alpha,lw=lw)
axs[1].plot (thlvpfmn_micr_dry, zflim[1:-1],c=colors[7],label=terms[7],alpha=alpha,lw=lw)
axs[1].plot (thlvpfmn_radi_dry, zflim,c=colors[8],label=terms[8],alpha=alpha,lw=lw)
# axs[1].plot( thlvpfmn_resi_dry, zflim[2:-2],c='gray',label='Residual')
axs[1].set_xlabel(r"Contribution to $\theta_{lv_m}'$ tendency [K/s]")
axs[1].set_xlim((-5.5e-5,5.5e-5))
axs[1].set_title('Dry')

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

plt.savefig(sp+'/thlvpf_budget.pdf',bbox_inches='tight')

#%% Average thlvpp budget contributions over time dimension
tpltmin = 10.
tpltmax = 16.

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

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

#%% Average wthlvpf budget contributions over time dimension
tpltmin = 6.
tpltmax = 16.


terms = ['Tendency                               ',
         'Gradient production',
         'Vertical flux convergence',
         'Horizontal flux convergence',
         'Buoyancy',
         'Pressure gradient',
         'Subsidence',
         'SFS diffusion'
         ]

colors = ['black',
          'cadetblue',
          'lightsteelblue',
          'olivedrab',
          'darkolivegreen',
          'rosybrown',
          'sienna',
          'goldenrod']

alpha=0.75
lw=2

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

wthlvpfmn_tend_moist = np.mean(wthlvpf_tend_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_prod_moist = np.mean(wthlvpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_vdiv_moist = np.mean(wthlvpf_vdiv_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_hdiv_moist = np.mean(wthlvpf_hdiv_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_buoy_moist = np.mean(wthlvpf_buoy_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_pres_moist = np.mean(wthlvpf_pres_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_subs_moist = np.mean(wthlvpf_subs_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_diff_moist = np.mean(wthlvpf_diff_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_budg_moist = (-wthlvpfmn_prod_moist[1:-2] - wthlvpfmn_vdiv_moist[1:-1]
                        -wthlvpfmn_hdiv_moist[1:-2] - wthlvpfmn_subs_moist[1:-2]
                        +wthlvpfmn_buoy_moist[1:-2] - wthlvpfmn_pres_moist[1:-2])
                        # +wthlvpfmn_diff_moist)
wthlvpfmn_resi_moist = wthlvpfmn_tend_moist[1:-2] - wthlvpfmn_budg_moist
wthlvpfmn_vdiv_moist= wthlvpfmn_vdiv_moist[1:-1] - wthlvpfmn_resi_moist

wthlvpfmn_tend_dry = np.mean(wthlvpf_tend_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_prod_dry = np.mean(wthlvpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_vdiv_dry = np.mean(wthlvpf_vdiv_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_hdiv_dry = np.mean(wthlvpf_hdiv_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_buoy_dry = np.mean(wthlvpf_buoy_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_pres_dry = np.mean(wthlvpf_pres_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_subs_dry = np.mean(wthlvpf_subs_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_diff_dry = np.mean(wthlvpf_diff_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_budg_dry = (-wthlvpfmn_prod_dry[1:-2] - wthlvpfmn_vdiv_dry[1:-1]
                      -wthlvpfmn_hdiv_dry[1:-2] - wthlvpfmn_subs_dry[1:-2]
                      +wthlvpfmn_buoy_dry[1:-2] - wthlvpfmn_pres_dry[1:-2])
                      # +wthlvpfmn_diff_dry)
wthlvpfmn_resi_dry = wthlvpfmn_tend_dry[1:-2] - wthlvpfmn_budg_dry
wthlvpfmn_vdiv_dry = wthlvpfmn_vdiv_dry[1:-1] - wthlvpfmn_resi_dry

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
axs[0].plot( wthlvpfmn_tend_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[0])
axs[0].plot(-wthlvpfmn_prod_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[1])
axs[0].plot(-wthlvpfmn_vdiv_moist, zflim[2:-3], alpha=alpha, lw=lw, c=colors[2])
axs[0].plot(-wthlvpfmn_hdiv_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[3])
axs[0].plot( wthlvpfmn_buoy_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[4])
axs[0].plot(-wthlvpfmn_pres_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[5])
axs[0].plot(-wthlvpfmn_subs_moist, zflim[1:-1], alpha=alpha, lw=lw, c=colors[6])
# axs[0].plot( wthlvpfmn_diff_moist, zflim[2:-3],c='skyblue')
# axs[0].plot( wthlvpfmn_resi_moist, zflim[2:-3],c='lightgray')
axs[0].set_xlabel(r"Contribution to $\left( w_s'\theta_{lv_s}'\right)_m - \overline{w_s'\theta_{lv_s}'}$ tendency, [K m/s$^2$]")
axs[0].set_xlim((-3.2e-4,3.2e-4))
axs[0].set_title('Moist')


axs[1].plot( wthlvpfmn_tend_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[0],label=terms[0])
axs[1].plot(-wthlvpfmn_prod_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[1],label=terms[1])
axs[1].plot(-wthlvpfmn_vdiv_dry, zflim[2:-3], alpha=alpha, lw=lw, c=colors[2],label=terms[2])
axs[1].plot(-wthlvpfmn_hdiv_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[3],label=terms[3])
axs[1].plot( wthlvpfmn_buoy_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[4],label=terms[4])
axs[1].plot(-wthlvpfmn_pres_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[5],label=terms[5])
axs[1].plot(-wthlvpfmn_subs_dry, zflim[1:-1], alpha=alpha, lw=lw, c=colors[6],label=terms[6])
# axs[1].plot( wthlvpfmn_diff_dry, zflim[2:-3],c='skyblue',label='SFS diffusion')
# axs[1].plot( wthlvpfmn_resi_dry, zflim[2:-3],c='lightgray',label='Residual')
axs[1].set_xlabel(r"Contribution to $\left( w_s'\theta_{lv_s}'\right)_m - \overline{w_s'\theta_{lv_s}'}$ tendency, [K m/s$^2$]")
axs[1].set_xlim((-3.2e-4,3.2e-4))
axs[1].set_title('Dry')


axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig(sp+'/wthlvpf_budget.pdf',bbox_inches='tight')
plt.show()

#%% WTG-based model of moisture variance production

tpltmin = 6.
tpltmax = 16.

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

colors = ['black',
          'maroon',
          'peru',
          'olive',
          'seagreen',
          ]
lw=2
alpha=0.5

# Vertical velocities

# Exact
wffmn_moist = np.mean(wff_moist_time[itpltmin:itpltmax],axis=0)
wffmn_dry = np.mean(wff_dry_time[itpltmin:itpltmax],axis=0)

# w model with actual thlvpf_vdiv (including budget residual)
thlvpf_budg_moist = (-thlvpf_prod_moist_time[:,1:-1] - thlvpf_vdiv_moist_time[:,1:-1]
                      -thlvpf_hdiv_moist_time[:,1:-1] - thlvpf_subs_moist_time[:,1:-1]
                      +thlvpf_diff_moist_time)
thlvpf_resi_moist = thlvpf_tend_moist_time[:,1:-1] -  thlvpf_budg_moist
thlvpf_vdiv_moist = thlvpf_vdiv_moist_time[:,1:-1] - thlvpf_resi_moist

thlvpf_budg_dry = (-thlvpf_prod_dry_time[:,1:-1] - thlvpf_vdiv_dry_time[:,1:-1]
                      -thlvpf_hdiv_dry_time[:,1:-1] - thlvpf_subs_dry_time[:,1:-1]
                      +thlvpf_diff_dry_time)
thlvpf_resi_dry = thlvpf_tend_dry_time[:,1:-1] -  thlvpf_budg_dry
thlvpf_vdiv_dry = thlvpf_vdiv_dry_time[:,1:-1] - thlvpf_resi_dry


wff_moist_wtg = -thlvpf_vdiv_moist/Gamma_thlv[:,1:-1]
wff_dry_wtg = -thlvpf_vdiv_dry/Gamma_thlv[:,1:-1]

wffmn_moist_wtg = np.mean(wff_moist_wtg[itpltmin:itpltmax,:],axis=0)
wffmn_dry_wtg = np.mean(wff_dry_wtg[itpltmin:itpltmax,:],axis=0)

# w model with simple reliance on qtpf
wff_moist_mod = -qtpf_prod_moist_time/Gamma_qt
wff_dry_mod = -qtpf_prod_dry_time/Gamma_qt

wffmn_moist_mod = np.mean(wff_moist_mod[itpltmin:itpltmax,:],axis=0)
wffmn_dry_mod = np.mean(wff_dry_mod[itpltmin:itpltmax,:],axis=0)

# w model with thlpf_vdiv
thlpf_vdiv_moist_time = thlvpf_vdiv_moist_time - 0.608*thl_av_time[:,1:-1]*qtpf_vdiv_moist_time
thlpf_vdiv_dry_time = thlvpf_vdiv_dry_time - 0.608*thl_av_time[:,1:-1]*qtpf_vdiv_dry_time

Gamma_thl = (thl_av_time[:,1:] - thl_av_time[:,:-1])/dzh
Gamma_thl_f = (Gamma_thl[:,1:] + Gamma_thl[:,:-1])*0.5

wff_moist_thl = -thlpf_vdiv_moist_time/Gamma_thl_f
wff_dry_thl = -thlpf_vdiv_dry_time/Gamma_thl_f

wffmn_moist_thl = np.mean(wff_moist_thl[itpltmin:itpltmax,:],axis=0)
wffmn_dry_thl = np.mean(wff_dry_thl[itpltmin:itpltmax,:],axis=0)

# w model with qlpf_vdiv
wff_moist_ql = 7*thl_av_time[:,1:-1]*qlpf_vdiv_moist_time/Gamma_thl_f
wff_dry_ql = 7*thl_av_time[:,1:-1]*qlpf_vdiv_dry_time/Gamma_thl_f

wffmn_moist_ql =  np.mean(wff_moist_ql[itpltmin:itpltmax,:],axis=0)
wffmn_dry_ql =  np.mean(wff_dry_ql[itpltmin:itpltmax,:],axis=0)

# Moisture variance production
qtpfmn_prod_moist_wex = np.mean(qtpf_prod_moist_wex_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_wex = np.mean(qtpf_prod_dry_wex_time[itpltmin:itpltmax,:],axis=0)

qtpfmn_prod_moist = np.mean(-qtpf_prod_moist_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry = np.mean(-qtpf_prod_dry_time[itpltmin:itpltmax,:],axis=0)

qtpf_prod_moist_wtg_time = wff_moist_wtg*Gamma_qt[:,1:-1]
qtpf_prod_dry_wtg_time = wff_dry_wtg*Gamma_qt[:,1:-1]

qtpfmn_prod_moist_wtg = np.mean(qtpf_prod_moist_wtg_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_wtg = np.mean(qtpf_prod_dry_wtg_time[itpltmin:itpltmax,:],axis=0)

qtpf_prod_moist_thl_time = wff_moist_thl*Gamma_qt
qtpf_prod_dry_thl_time = wff_dry_thl*Gamma_qt

qtpfmn_prod_moist_thl = np.mean(qtpf_prod_moist_thl_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_thl = np.mean(qtpf_prod_dry_thl_time[itpltmin:itpltmax,:],axis=0)

qtpf_prod_moist_ql_time = wff_moist_ql*Gamma_qt
qtpf_prod_dry_ql_time = wff_dry_ql*Gamma_qt

qtpfmn_prod_moist_ql =  np.mean(qtpf_prod_moist_ql_time[itpltmin:itpltmax,:],axis=0)
qtpfmn_prod_dry_ql =  np.mean(qtpf_prod_dry_ql_time[itpltmin:itpltmax,:],axis=0)

# w plot
fig,axs = plt.subplots(nrows=2,ncols=2,sharey=True,figsize=(10,10),squeeze=False)
axs[0,0].plot(wffmn_moist, zflim, c=colors[0],linewidth=lw,alpha=alpha)
axs[0,0].plot(wffmn_moist_wtg, zflim[2:-2], c=colors[1],linewidth=lw,alpha=alpha)
axs[0,0].plot(wffmn_moist_thl, zflim[1:-1], c=colors[2])
axs[0,0].plot(wffmn_moist_ql, zflim[1:-1], c=colors[4])
# axs[0,0].plot(wffmn_moist_mod, zflim[1:-1], c='black',linestyle='-.')
axs[0,0].axvline(0,color='gray',linestyle='dotted')
axs[0,0].set_xlim((-0.012,0.012))
axs[0,0].set_xlabel(r"$w_m'$ [m/s]")
axs[0,0].set_title(r"Moist")
axs[0,0].annotate('a)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[0,0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0,1].plot(wffmn_dry, zflim,  c=colors[0],linewidth=lw,alpha=alpha, label=r"Ground truth")
axs[0,1].plot(wffmn_dry_wtg, zflim[2:-2],  c=colors[1],linewidth=lw,alpha=alpha, label=r"$\theta_{lv}$-based WTG model")
axs[0,1].plot(wffmn_dry_thl, zflim[1:-1],  c=colors[2],linewidth=lw,alpha=alpha, label=r"$\theta_{l}$-based WTG model")
axs[0,1].plot(wffmn_dry_ql, zflim[1:-1],  c=colors[4],linewidth=lw,alpha=alpha, label=r"$q_l$-based WTG model")
# axs[0,1].plot(wffmn_dry_mod, zflim[1:-1], c='black',linestyle='-.')
axs[0,1].axvline(0,color='gray',linestyle='dotted')
axs[0,1].set_xlim((-0.012,0.012))
axs[0,1].set_xlabel(r"$w_m'$ [m/s]")
axs[0,1].set_title(r"Dry")
axs[0,1].annotate('b)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[0,1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[0,0].set_ylabel(r'Height [m]')
axs[0,1].legend(loc='upper left',bbox_to_anchor=(1,1))

# Moisture variance production plot
axs[1,0].plot(qtpfmn_prod_moist_wex, zflim[1:-1], c=colors[0],linewidth=lw,alpha=alpha)
axs[1,0].plot(qtpfmn_prod_moist_wtg, zflim[2:-2], c=colors[1],linewidth=lw,alpha=alpha)
axs[1,0].plot(qtpfmn_prod_moist_thl, zflim[1:-1], c=colors[2],linewidth=lw,alpha=alpha)
axs[1,0].plot(qtpfmn_prod_moist_ql,  zflim[1:-1], c=colors[4],linewidth=lw,alpha=alpha)
# axs[1,0].plot(qtpfmn_prod_moist, zflim[1:-1], c='black',linestyle='-.')
axs[1,0].axvline(0,color='gray',linestyle='dotted')
axs[1,0].set_xlim((-9e-8,9e-8))
axs[1,0].set_xlabel(r"$w_m'\Gamma_{q_t}$ [kg/kg/s]")
axs[1,0].annotate('c)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[1,0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[1,1].plot(qtpfmn_prod_dry_wex, zflim[1:-1], c=colors[0],linewidth=lw,alpha=alpha)
axs[1,1].plot(qtpfmn_prod_dry_wtg, zflim[2:-2], c=colors[1],linewidth=lw,alpha=alpha)
axs[1,1].plot(qtpfmn_prod_dry_thl, zflim[1:-1], c=colors[2],linewidth=lw,alpha=alpha)
axs[1,1].plot(qtpfmn_prod_dry_ql,  zflim[1:-1], c=colors[4],linewidth=lw,alpha=alpha)
axs[1,1].axvline(0,color='gray',linestyle='dotted')
axs[1,1].set_xlim((-9e-8,9e-8))
axs[1,1].set_xlabel(r"$w_m'\Gamma_{q_t}$ [kg/kg/s]")
axs[1,1].annotate('d)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[1,1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

axs[1,0].set_ylabel(r'Height [m]')

plt.savefig(sp+'/wpf_qtpfprod_wtg.pdf',bbox_inches='tight')


#%% Vertically integrated statistics
tpltmin = 2.
tpltmax = 36.
dit = 0.5 # Rounds to closest multiple of dt in time
dtav = 1.0 # Average around each plotted time step
alpha = 0.5
lw=2

terms = ['Tendency',
         'Gradient production',
         'Vertical flux convergence',
         'Horizontal flux convergence',
         'Subsidence',
         'SFS diffusion',
         'Residual',
         'Precipitation',
         'Radiation'
         ]

colors = ['black',
          'cadetblue',
          'lightsteelblue',
          'olivedrab',
          'sienna',
          'goldenrod',
          'lightgray',
          'maroon',
          'midnightblue']

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
idtav  = int(round(dtav/2/(time[1]-time[0])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)
  

qtpfi_moist = vint(qtpf_moist_time,rhobfi,zflim,plttime_var)
qtpfi_dry = vint(qtpf_dry_time,rhobfi,zflim,plttime_var)

# Tendency
qtpfi_tend_moist = vint(qtpf_tend_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
qtpfi_tend_dry = vint(qtpf_tend_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

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
qtpfi_diff_moist = vint(qtpf_diff_moist_time,rhobfi[2:-2],zflim[2:-2],plttime_var)
qtpfi_diff_dry = vint(qtpf_diff_dry_time,rhobfi[2:-2],zflim[2:-2],plttime_var)

# Moistening through precipitation
qtpfi_micr_moist = vint(qtpf_micr_moist_time,rhobfi,zflim,plttime_var)
qtpfi_micr_dry = vint(qtpf_micr_dry_time,rhobfi,zflim,plttime_var)

# Estimate residual
qtpfi_resid_moist = qtpfi_tend_moist + qtpfi_prod_wex_moist + qtpfi_vdiv_moist + qtpfi_hdiv_moist + qtpfi_subs_moist - qtpfi_diff_moist
qtpfi_resid_dry = qtpfi_tend_dry + qtpfi_prod_wex_dry + qtpfi_vdiv_dry + qtpfi_hdiv_dry + qtpfi_subs_dry - qtpfi_diff_dry

# And include it in the horizontal advection, which is the worst-estimated
# qtpfi_hdiv_moist-=qtpfi_resid_moist
# qtpfi_hdiv_dry-=qtpfi_resid_dry

# Temporal plot
fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,10/3))
axs[0].plot(time[plttime_var],qtpfi_tend_moist,c=colors[0],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],-qtpfi_prod_wex_moist,c=colors[1],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],-qtpfi_vdiv_moist,c=colors[2],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],-qtpfi_hdiv_moist,c=colors[3],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],-qtpfi_subs_moist,c=colors[4],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],qtpfi_diff_moist,c=colors[5],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],qtpfi_resid_moist,c=colors[6],alpha=alpha,lw=lw)
axs[0].plot(time[plttime_var],qtpfi_micr_moist,c=colors[7],alpha=alpha,lw=lw)
axs[0].set_xlabel('Time [hr]')
axs[0].set_title('Moist')

axs[1].plot(time[plttime_var],qtpfi_tend_dry,c=colors[0],label=terms[0],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],-qtpfi_prod_wex_dry,c=colors[1],label=terms[1],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],-qtpfi_vdiv_dry,c=colors[2],label=terms[2],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],-qtpfi_hdiv_dry,c=colors[3],label=terms[3],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],-qtpfi_subs_dry,c=colors[4],label=terms[4],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],qtpfi_diff_dry,c=colors[5],label=terms[5],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],qtpfi_resid_dry,c=colors[6],label=terms[6],alpha=alpha,lw=lw)
axs[1].plot(time[plttime_var],qtpfi_micr_dry,c=colors[7],label=terms[7],alpha=alpha,lw=lw)
axs[1].set_xlabel('Time [hr]')
axs[1].set_title('Dry')

axs[0].set_ylabel('Mesoscale moistening rate [kg/m$^2$/s]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))
plt.savefig(sp+'/qtpf_budget_int.pdf',bbox_inches='tight')

# Model evaluation
# fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
# axs[0].plot(time[plttime_var],qtpfi_tend_moist,c='midnightblue')
# axs[0].plot(time[plttime_var],-qtpfi_prod_wex_moist,c='darkseagreen')
# axs[0].plot(time[plttime_var],qtpfi_prod_moist,c='darkseagreen',linestyle='dotted')
# axs[0].set_xlabel('Time [hr]')
# axs[0].set_title('Moist region')

# axs[1].plot(time[plttime_var],qtpfi_tend_dry,c='midnightblue',label=terms[0])
# axs[1].plot(time[plttime_var],-qtpfi_prod_wex_dry,c='darkseagreen',label=terms[1])
# axs[1].plot(time[plttime_var],qtpfi_prod_dry,c='darkseagreen',linestyle='dotted',label='Modelled '+terms[1])
# axs[1].set_xlabel('Time [hr]')
# axs[1].set_title('Dry region')

# axs[0].set_ylabel('Large-scale moistening rate [kg/kg/s]')
# axs[1].legend(loc='best',bbox_to_anchor=(1,1))

# ### And now for thlv
# thlvpfi_moist = vint(thlvpf_moist_time,rhobfi,zflim,plttime_var)
# thlvpfi_dry = vint(thlvpf_dry_time,rhobfi,zflim,plttime_var)

# # thlv gradient production
# thlvpfi_prod_moist = vint(thlvpf_prod_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
# thlvpfi_prod_dry = vint(thlvpf_prod_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# # heating through anomalous vertical fluxes
# # FIXME offset zf in integration by 1 from field
# thlvpfi_vdiv_moist = vint(thlvpf_vdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
# thlvpfi_vdiv_dry = vint(thlvpf_vdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# # Heating through horizontal advection
# thlvpfi_hdiv_moist = vint(thlvpf_hdiv_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
# thlvpfi_hdiv_dry = vint(thlvpf_hdiv_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# # Heating through subsidence
# thlvpfi_subs_moist = vint(thlvpf_subs_moist_time,rhobfi[1:-1],zflim[1:-1],plttime_var)
# thlvpfi_subs_dry = vint(thlvpf_subs_dry_time,rhobfi[1:-1],zflim[1:-1],plttime_var)

# # Heating through SFS diffusion
# # thlvpfi_diff_moist = vint(thlvpf_diff_moist_time,rhobfi,zflim,plttime_var)
# # thlvpfi_diff_dry = vint(thlvpf_diff_dry_time,rhobfi,zflim,plttime_var)


# # Fit the haeting fluctuation's evolution
# [exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * x**a, 
#                                        time[plttime_var]*3600, 
#                                        thlvpfi_moist,
#                                        p0=[1,0])

# [exp_dry,fac_dry], cov = curve_fit(lambda x, a, b: b * x**a, 
#                                        time[plttime_var]*3600, 
#                                        thlvpfi_dry,
#                                        p0=[-3,0])

# # And differentiate to estimate its tendency
# thlvpfi_tend_moist = fac_moist*exp_moist*(time[plttime_var]*3600)**(exp_moist-1)
# thlvpfi_tend_dry = fac_dry*exp_dry*(time[plttime_var]*3600)**(exp_dry-1)

# # Estimate residual
# thlvpfi_resid_moist = thlvpfi_tend_moist + thlvpfi_prod_moist + thlvpfi_vdiv_moist + thlvpfi_hdiv_moist + thlvpfi_subs_moist
# thlvpfi_resid_dry = thlvpfi_tend_dry + thlvpfi_prod_dry + thlvpfi_vdiv_dry + thlvpfi_hdiv_dry + thlvpfi_subs_dry

# fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(15,5))
# axs[0].plot(time[plttime_var],thlvpfi_tend_moist,c='midnightblue')
# # axs[0].plot(time[plttime_var],qtpfi_prod_moist,c='darkseagreen')
# axs[0].plot(time[plttime_var],-thlvpfi_prod_moist,c='maroon')
# axs[0].plot(time[plttime_var],-thlvpfi_vdiv_moist,c='peachpuff')
# axs[0].plot(time[plttime_var],-thlvpfi_hdiv_moist,c='olive')
# axs[0].plot(time[plttime_var],-thlvpfi_subs_moist,c='skyblue')
# axs[0].plot(time[plttime_var],thlvpfi_resid_moist,c='slategray')
# axs[0].set_xlabel('Time [hr]')
# axs[0].set_title('Moist region')

# axs[1].plot(time[plttime_var],qtpfi_tend_dry,c='midnightblue',label=r"$\frac{\partial\langle\tilde{\theta_{lv}'}\rangle}{\partial t}$")
# # axs[1].plot(time[plttime_var],qtpfi_prod_dry,c='darkseagreen',label=r"$F_{\langle\tilde{q_t'}\rangle}$")
# axs[1].plot(time[plttime_var],-thlvpfi_prod_dry,c='maroon',label=r"$-\tilde{w'}\frac{\partial \overline{\theta_{lv}}}{\partial z}$")
# axs[1].plot(time[plttime_var],-thlvpfi_vdiv_dry,c='peachpuff',label=r"$-\frac{1}{\rho_0}\frac{\partial}{\partial z}\left(\rho_0\left(\widetilde{w'\theta_{lv}'}-\overline{w'\theta_{lv}'}\right)\right)$")
# axs[1].plot(time[plttime_var],-thlvpfi_hdiv_dry,c='olive',label=r"$-\frac{\partial}{\partial x_{hj}}\left(\widetilde{u_{hj}'\theta_{lv}'}\right)$")
# axs[1].plot(time[plttime_var],-thlvpfi_subs_dry,c='skyblue',label=r"$-\overline{w_{LS}}\frac{\partial \tilde{\theta_{lv}'}}{\partial z}$")
# axs[1].plot(time[plttime_var],thlvpfi_resid_dry,c='slategray',label=r"Residual")
# axs[1].set_xlabel('Time [hr]')
# axs[1].set_title('Dry region')

# axs[0].set_ylabel('Large-scale heating rate [K/s]')
# axs[1].legend(loc='best',bbox_to_anchor=(1,1))

#%% Fluxes and fluctuations of thv

# Time to average over
tpltmin = 2.
tpltmax = 6.

terms0 = [r"$\theta_{v_m}'$",
          r"$\theta_{lv_m}'$",
          r"$\theta_{l_m}'$",
          r"$0.608\overline{\theta_l}q_{t_m}'$",
          r"$7\overline{\theta_l}q_{l_m}'$",
         ]

terms1 = [r"$\left(w'\theta_v'\right)_m$",
          r"$\left(w'\theta_{lv}'\right)_m$",
          r"$\left(w'\theta_l'\right)_m$",
          r"$0.608\overline{\theta_l}\left(w'q_t'\right)_m$",
          r"$7\overline{\theta_l}\left(w'q_l'\right)_m$",
         ]

colors = ['black',
          'maroon',
          'peru',
          'olive',
          'seagreen',
          ]

alpha = 0.6
lw = 2

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

itpltmin1d = np.where(time1d>=tpltmin*3600)[0][0]
itpltmax1d = np.where(time1d<tpltmax*3600)[0][-1]+1

# Contributions to thv
thvpfmn_moist = np.mean(thvpf_moist_time[itpltmin:itpltmax,:],axis=0)
thvpfmn_dry = np.mean(thvpf_dry_time[itpltmin:itpltmax,:],axis=0)

thlvpfmn_moist = np.mean(thlvpf_moist_time[itpltmin:itpltmax,:],axis=0)
thlvpfmn_dry = np.mean(thlvpf_dry_time[itpltmin:itpltmax,:],axis=0)

thlpfmn_moist = np.mean(thlpf_moist_time[itpltmin:itpltmax,:],axis=0)
thlpfmn_dry = np.mean(thlpf_dry_time[itpltmin:itpltmax,:],axis=0)

a2qtpfmn_moist = np.mean(0.608*thl_av_time[itpltmin:itpltmax,:]*qtpf_moist_time[itpltmin:itpltmax,:],axis=0)
a2qtpfmn_dry = np.mean(0.608*thl_av_time[itpltmin:itpltmax,:]*qtpf_dry_time[itpltmin:itpltmax,:],axis=0)

a3qlpfmn_moist = np.mean(7*thl_av_time[itpltmin:itpltmax,:]*qlpf_moist_time[itpltmin:itpltmax,:],axis=0)
a3qlpfmn_dry = np.mean(7*thl_av_time[itpltmin:itpltmax,:]*qlpf_dry_time[itpltmin:itpltmax,:],axis=0)

# Contributions to wthv
wthvpfmn_moist = np.mean(wthlvpf_moist_time[itpltmin:itpltmax,:]+
                         7*thl_av_time[itpltmin:itpltmax,:]*wqlpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthvpfmn_dry = np.mean(wthlvpf_dry_time[itpltmin:itpltmax,:]+
                         7*thl_av_time[itpltmin:itpltmax,:]*wqlpf_dry_time[itpltmin:itpltmax,:],axis=0)
wthvpmn_av = np.mean(dl.load_wthvrav(izmin, izmax)[itpltmin1d:itpltmax1d,:],axis=0)

wthlvpfmn_moist = np.mean(wthlvpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_dry = np.mean(wthlvpf_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpmn_av = np.mean(wthlvp_av_time[itpltmin:itpltmax,:],axis=0)

wthlpfmn_moist = np.mean(wthlpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthlpfmn_dry = np.mean(wthlpf_dry_time[itpltmin:itpltmax,:],axis=0)
wthlpmn_av = np.mean(dl.load_wthlrav(izmin, izmax)[itpltmin1d:itpltmax1d,:],axis=0)

a2wqtpfmn_moist = np.mean(0.608*thl_av_time[itpltmin:itpltmax,:]*wqtpf_moist_time[itpltmin:itpltmax,:],axis=0)
a2wqtpfmn_dry = np.mean(0.608*thl_av_time[itpltmin:itpltmax,:]*wqtpf_dry_time[itpltmin:itpltmax,:],axis=0)
a2wqtpmn_av = np.mean(0.608*thl_av_1d[itpltmin1d:itpltmax1d,:]*dl.load_wqtrav(izmin,izmax)[itpltmin1d:itpltmax1d,:],axis=0)

a3wqlpfmn_moist = np.mean(7*thl_av_time[itpltmin:itpltmax,:]*wqlpf_moist_time[itpltmin:itpltmax,:],axis=0)
a3wqlpfmn_dry = np.mean(7*thl_av_time[itpltmin:itpltmax,:]*wqlpf_dry_time[itpltmin:itpltmax,:],axis=0)
a3wqlpmn_av = np.mean(7*thl_av_1d[itpltmin1d:itpltmax1d,:]*dl.load_wqlrav(izmin,izmax)[itpltmin1d:itpltmax1d,:],axis=0)

# thvpf plot
fig,axs = plt.subplots(nrows=1,ncols=2,sharey=True,figsize=(10,5), squeeze=False)
axs[0,0].plot(thvpfmn_moist, zflim,c=colors[0],alpha=alpha,lw=lw)
axs[0,0].plot(thlvpfmn_moist, zflim,c=colors[1],alpha=alpha,lw=lw)
axs[0,0].plot(thlpfmn_moist, zflim,c=colors[2],alpha=alpha,lw=lw)
axs[0,0].plot(a2qtpfmn_moist, zflim,c=colors[3],alpha=alpha,lw=lw)
axs[0,0].plot(a3qlpfmn_moist, zflim,c=colors[4],alpha=alpha,lw=lw)
axs[0,0].set_xlabel(r"Contribution to $\theta_{v_m}'$ [K]")
axs[0,0].set_ylabel(r'Height [m]')
axs[0,0].axvline(0,color='gray',linestyle='dotted')
axs[0,0].set_xlim((-0.072,0.072))
axs[0,0].set_title('Moist')
# axs[0,0].annotate('a)', (0.05,0.9), xycoords='axes fraction', fontsize=14)

axs[0,1].plot(thvpfmn_dry, zflim,c=colors[0],alpha=alpha,lw=lw,label=terms0[0])
axs[0,1].plot(thlvpfmn_dry, zflim,c=colors[1],alpha=alpha,lw=lw,label=terms0[1])
axs[0,1].plot(thlpfmn_dry, zflim,c=colors[2],alpha=alpha,lw=lw,label=terms0[2])
axs[0,1].plot(a2qtpfmn_dry, zflim,c=colors[3],alpha=alpha,lw=lw,label=terms0[3])
axs[0,1].plot(a3qlpfmn_dry, zflim,c=colors[4],alpha=alpha,lw=lw,label=terms0[4])
axs[0,1].set_xlabel(r"Contribution to $\theta_{v_m}'$ [K]")
axs[0,1].axvline(0,color='gray',linestyle='dotted')
axs[0,1].set_xlim((-0.072,0.072))
axs[0,1].set_title('Dry')
# axs[0,1].annotate('b)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[0,1].legend(bbox_to_anchor=(1,1),loc='upper left')

plt.savefig(sp+'/thv_decomposition.pdf',bbox_inches='tight')

# wthvpf plot
fig,axs = plt.subplots(nrows=1,ncols=2,sharey=True,figsize=(10,5), squeeze=False)
axs[0,0].plot(wthvpfmn_moist, zflim,c=colors[0],alpha=alpha,lw=lw)
axs[0,0].plot(wthlvpfmn_moist, zflim,c=colors[1],alpha=alpha,lw=lw)
axs[0,0].plot(wthlpfmn_moist, zflim,c=colors[2],alpha=alpha,lw=lw)
axs[0,0].plot(a2wqtpfmn_moist, zflim,c=colors[3],alpha=alpha,lw=lw)
axs[0,0].plot(a3wqlpfmn_moist, zflim,c=colors[4],alpha=alpha,lw=lw)
axs[0,0].plot(wthvpmn_av, zflim,c=colors[0],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,0].plot(wthlvpmn_av, zflim,c=colors[1],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,0].plot(wthlpmn_av, zflim,c=colors[2],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,0].plot(a2wqtpmn_av, zflim,c=colors[3],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,0].plot(a3wqlpmn_av, zflim,c=colors[4],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,0].set_xlabel(r"Contribution to $\left(w'\theta_{v}'\right)_m$ [Km/s]")
axs[0,0].set_ylabel(r'Height [m]')
axs[0,0].axvline(0,color='gray',linestyle='dotted')
axs[0,0].set_xlim((-0.047,0.047))
axs[0,0].set_title('Moist')
# axs[1,0].annotate('c)', (0.05,0.9), xycoords='axes fraction', fontsize=14)

axs[0,1].plot(wthvpfmn_dry, zflim,c=colors[0],alpha=alpha,lw=lw,label=terms1[0])
axs[0,1].plot(wthlvpfmn_dry, zflim,c=colors[1],alpha=alpha,lw=lw,label=terms1[1])
axs[0,1].plot(wthlpfmn_dry, zflim,c=colors[2],alpha=alpha,lw=lw,label=terms1[2])
axs[0,1].plot(a2wqtpfmn_dry, zflim,c=colors[3],alpha=alpha,lw=lw,label=terms1[3])
axs[0,1].plot(a3wqlpfmn_dry, zflim,c=colors[4],alpha=alpha,lw=lw,label=terms1[4])
axs[0,1].plot(wthvpmn_av, zflim,c=colors[0],alpha=alpha,lw=lw-0.5,linestyle='--',label='Slab average')
axs[0,1].plot(wthlvpmn_av, zflim,c=colors[1],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,1].plot(wthlpmn_av, zflim,c=colors[2],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,1].plot(a2wqtpmn_av, zflim,c=colors[3],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,1].plot(a3wqlpmn_av, zflim,c=colors[4],alpha=alpha,lw=lw-0.5,linestyle='--')
axs[0,1].set_xlabel(r"Contribution to $\left(w'\theta_{v}'\right)_m$ [Km/s]")
axs[0,1].axvline(0,color='gray',linestyle='dotted')
axs[0,1].set_xlim((-0.047,0.047))
axs[0,1].set_title('Dry')
# axs[0,1].annotate('d)', (0.05,0.9), xycoords='axes fraction', fontsize=14)
axs[0,1].legend(bbox_to_anchor=(1,1),loc='upper left')

plt.savefig(sp+'/wthv_decomposition.pdf',bbox_inches='tight')


#%% Flux scale decomposition

# Time to average over
tpltmin = 6.
tpltmax = 10.

labs = [r"$(w'\theta_{lv}')_m$",
        r"$(w_m'\theta_{lv_m}')_m$",
        r"$(w_m'\theta_{lv_s}')_m + (w_s'\theta_{lv_m}')_m$",
        r"$(w_s'\theta_{lv_s}')_m$",
        r"$(w'q_l')_m$",
        r"$(w_m'q_{l_m}')_m$",
        r"$(w_m'q_{l_s}')_m + (w_s'q_{l_m}')_m$",
        r"$(w_s'q_{l_s}')_m$"
        ]

col = ['maroon',
       'indianred',
       'darksalmon',
       'papayawhip',
       'seagreen',
       'mediumseagreen',
       'yellowgreen',
       'palegoldenrod'
       ]

alpha = 1
lw = 2

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1

wthlvpfmn_moist = np.mean(wthlvpf_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_dry = np.mean(wthlvpf_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_l_moist = np.mean(wthlvpf_l_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_l_dry = np.mean(wthlvpf_l_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_c_moist = np.mean(wthlvpf_c_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_c_dry = np.mean(wthlvpf_c_dry_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_r_moist = np.mean(wthlvpf_r_moist_time[itpltmin:itpltmax,:],axis=0)
wthlvpfmn_r_dry = np.mean(wthlvpf_r_dry_time[itpltmin:itpltmax,:],axis=0)

wqlpfmn_moist   = -7*np.mean((thl_av_time*wqlpf_moist_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_dry     = -7*np.mean((thl_av_time*wqlpf_dry_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_l_moist = -7*np.mean((thl_av_time*wqlpf_l_moist_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_l_dry   = -7*np.mean((thl_av_time*wqlpf_l_dry_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_c_moist = -7*np.mean((thl_av_time*wqlpf_c_moist_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_c_dry   = -7*np.mean((thl_av_time*wqlpf_c_dry_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_r_moist = -7*np.mean((thl_av_time*wqlpf_r_moist_time)[itpltmin:itpltmax,:],axis=0)
wqlpfmn_r_dry   = -7*np.mean((thl_av_time*wqlpf_r_dry_time)[itpltmin:itpltmax,:],axis=0)

fig,axs = plt.subplots(nrows=1,ncols=2,sharey=True,figsize=(10,5))
axs[0].plot(wthlvpfmn_moist,zflim,alpha=alpha,lw=lw,c=col[0])
axs[0].plot(wthlvpfmn_l_moist,zflim,alpha=alpha,lw=lw,c=col[1])
axs[0].plot(wthlvpfmn_c_moist,zflim,alpha=alpha,lw=lw,c=col[2])
axs[0].plot(wthlvpfmn_r_moist,zflim,alpha=alpha,lw=lw,c=col[3])
# axs[0].plot(wqlpfmn_moist,zflim,alpha=alpha,lw=lw,c=col[4])
# axs[0].plot(wqlpfmn_l_moist,zflim,alpha=alpha,lw=lw,c=col[5])
# axs[0].plot(wqlpfmn_c_moist,zflim,alpha=alpha,lw=lw,c=col[6])
# axs[0].plot(wqlpfmn_r_moist,zflim,alpha=alpha,lw=lw,c=col[7])
axs[0].set_title('Moist')
axs[0].set_ylabel('z [m]')
axs[0].set_xlabel(r"$(w'\theta_{lv}')_m$ [Km/s]")
axs[0].set_xlim((-0.045,0.015))

axs[1].plot(wthlvpfmn_dry,zflim,alpha=alpha,lw=lw,c=col[0],label=labs[0])
axs[1].plot(wthlvpfmn_l_dry,zflim,alpha=alpha,lw=lw,c=col[1],label=labs[1])
axs[1].plot(wthlvpfmn_c_dry,zflim,alpha=alpha,lw=lw,c=col[2],label=labs[2])
axs[1].plot(wthlvpfmn_r_dry,zflim,alpha=alpha,lw=lw,c=col[3],label=labs[3])
# axs[1].plot(wqlpfmn_dry,zflim,alpha=alpha,lw=lw,c=col[4],label=labs[4])
# axs[1].plot(wqlpfmn_l_dry,zflim,alpha=alpha,lw=lw,c=col[5],label=labs[5])
# axs[1].plot(wqlpfmn_c_dry,zflim,alpha=alpha,lw=lw,c=col[6],label=labs[6])
# axs[1].plot(wqlpfmn_r_dry,zflim,alpha=alpha,lw=lw,c=col[7],label=labs[7])
axs[1].set_title('Dry')
axs[1].set_xlabel(r"$(w'\theta_{lv}')_m$ [Km/s]")
axs[1].set_xlim((-0.045,0.015))
axs[1].legend(bbox_to_anchor=(1,1),loc='upper left',ncol=1)

plt.savefig(sp+'/wthvm_scale_decomposition.pdf',bbox_inches='tight')

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
tpltmax = 24.
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
tpltmax = 24.
C = 0.4 # Constant of proportionality
fs = 14 # fontsize

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
plttime_plt = np.arange(itpltmax-itpltmin)

itpltmin1d = np.where(time1d/3600>=tpltmin)[0][0]
itpltmax1d = np.where(time1d/3600<tpltmax)[0][-1]+1

# OPTION TO USE A VELOCITY SCALING
# just copied from DALES namoptions now
wthl0 = 8e-3
wqt0 = 5.2e-5
thl0 = 299.1
grav = 9.81

# Using wstar
wstar = (grav/thl0*(wthl0+0.608*thl0*wqt0)*500)**(1/3)

# Calculate from data
zwmax = 1500
izwmax = np.where(zflim>zwmax)[0][0]
w2 = dl.load_w2tav(izmin, izmax)
wstar = np.sqrt(np.mean(w2))

# gammas
Gamrat = Gamma_qt/Gamma_thlv
Gamratz = (Gamrat[:,1:] - Gamrat[:,:-1])/(zflim[1] - zflim[0])

# Time filter
wthlvpf_moist_anom_pl = -wthlvpf_moist_anom[itpltmin:itpltmax,:]
wthlvpf_dry_anom_pl = -wthlvpf_dry_anom[itpltmin:itpltmax,:]
qtpf_moist_pl = qtpf_moist_time[itpltmin:itpltmax,:]
qtpf_dry_pl = qtpf_dry_time[itpltmin:itpltmax,:]
Gamratz_pl = np.mean(Gamratz[itpltmin:itpltmax,:],axis=0)

# Model
qtpf_mod = np.linspace(qtpf_dry_pl.min(),qtpf_moist_pl.max(),10)
wthlvpf_anom_mod = C*wstar*thl0*qtpf_mod

# Vertically integrated model
qtpfi_moist = vint(qtpf_moist_time[itpltmin:itpltmax,:],rhobfi, zflim,plttime=plttime_plt)
wthlvpf_anomi_moist = -vint(wthlvpf_moist_anom[itpltmin:itpltmax,2:-1]*Gamratz[itpltmin:itpltmax,:],rhobfi[2:-1], zflim[2:-1], plttime=plttime_plt)
qtpfi_dry = vint(qtpf_dry_time[itpltmin:itpltmax,:],rhobfi, zflim,plttime=plttime_plt)
wthlvpf_anomi_dry = -vint(wthlvpf_dry_anom[itpltmin:itpltmax,2:-1]*Gamratz[itpltmin:itpltmax,:],rhobfi[2:-1], zflim[2:-1], plttime=plttime_plt)

# Cloud layer-mean Gamratz
Gamratz_mn = np.mean(Gamratz_pl[20:40])
qtpfi_mod = np.linspace(qtpfi_dry.min(),qtpfi_moist.max(),10)
tau = 1. / (C*thl0*wstar*Gamratz_mn)
wthlvpf_anomi_mod = qtpfi_mod  / tau


fig,axs=plt.subplots(ncols=2,figsize=(10,5),squeeze=False)
axs[0,0].scatter(qtpf_moist_pl.flatten(),wthlvpf_moist_anom_pl.flatten(),c='k',s=0.1)
axs[0,0].scatter(qtpf_dry_pl.flatten(),wthlvpf_dry_anom_pl.flatten(),c='k',s=0.1)
axs[0,0].plot(qtpf_mod,wthlvpf_anom_mod,'k')
axs[0,0].set_ylabel(r"$-F_{{\theta_{lv}}_m'}$ [K m/s]",fontsize=fs)
axs[0,0].set_xlabel(r"$q_{t_m}'$ [kg/kg]",fontsize=fs)
axs[0,0].annotate('a)', (0.05,0.92), xycoords='axes fraction', fontsize=fs)

axs[0,1].scatter(qtpfi_moist,wthlvpf_anomi_moist,s=0.5,c='k')#time[plttime_plt],s=0.5)
axs[0,1].scatter(qtpfi_dry,wthlvpf_anomi_dry,s=0.5,c='k')#time[plttime_plt],s=0.5)
axs[0,1].plot(qtpfi_mod,wthlvpf_anomi_mod,'k')
axs[0,1].set_ylabel(r"$-\left\langle F_{{\theta_{lv}}_m'} \frac{\partial}{\partial z}\left(\frac{\Gamma_{q_t}}{\Gamma_{\theta_{lv}}}\right) \right\rangle$ [kg/$m^2$/s]", fontsize=fs)
axs[0,1].set_xlabel(r"$\left\langle q_{t_m}'\right\rangle$ [kg/m$^2$]", fontsize=fs)
axs[0,1].annotate(r"$\tau_{q_{t_m}'} =$ %.2f hr"%(tau/3600), (0.55,0.06), xycoords='axes fraction', fontsize=14)
axs[0,1].annotate('b)', (0.05,0.92), xycoords='axes fraction', fontsize=14)
plt.tight_layout()
plt.savefig(sp+'/qtpf_qtpfprod_model.png',dpi=300,bbox_inches='tight')
plt.show()


#%% Mean profiles

tpltmin = 1.
tpltmax = 21.
dit = 5 # Rounds to closest multiple of dt in time
fac=1e5 # Factor for plotting tendency magnitudes
tav=2 # time (hrs) around plotting time over which to average flux divergence
lw=2
alpha=0.5

terms = ['Tendency',
         'Vertical flux convergence',
         'Subsidence',
         'Large-scale moistening',
         'Large-scale heating'
         ]

colors = ['black',
          'lightsteelblue',
          'sienna',
          'cornflowerblue',
          'palevioletred'
          ]

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig = plt.figure(figsize=(5,5))
ax = plt.gca()
for i in range(len(plttime_var)):
    
    imin = np.argmin(wthlvp_av_time[plttime_var[i],:])
    ilab = r"Inversion base" if i==0 else None

    it1d = np.argmin(np.abs(time[plttime_var[i]]-time1d/3600))
    itmin = np.argmin(np.abs(time1d - (time1d[it1d]-tav*3600)))
    itmax = np.argmin(np.abs(time1d - (time1d[it1d]+tav*3600)))
    
    iqlbase = np.where(np.abs(ql_av_1d[it1d,:]) > 0)[0][0]
    qlbaselab = r"Cloud base" if i==0 else None
    
    iqltop = np.where(np.abs(ql_av_1d[it1d,:]) > 0)[0][-1]
    if iqltop > izmax-3:
        iqltop = izmax-3
    qltoplab = r"Cloud top" if i==0 else None
    
    # Tendency, averaged over tav around the current time
    ddt_thlv_av = np.mean(ddt_thlv_av_time[itmin:itmax],axis=0)
    ddt_qt_av = np.mean(ddt_qt_av_time[itmin:itmax],axis=0)
    
    # Flux divergence, averaged over tav around the current time
    ddz_wthlv_av = np.mean(-ddz_wthlv_av_time[itmin:itmax],axis=0)
    ddz_wqt_av = np.mean(-ddz_wqt_av_time[itmin:itmax],axis=0)
    
    # Subsidence
    qtavp_subs = (-wfls[izmin+1:izmax-1]*Gamma_qt[plttime_var[i],:])
    thlvavp_subs = (-wfls[izmin+1:izmax-1]*Gamma_thlv[plttime_var[i],:])
    
    # Large-scale drying
    qtavp_larq = dqdt_ls
    thlvavp_larq = 0.608*thl_av_time[plttime_var[i],:]*dqdt_ls

    # Large-scale cooling
    thlvavp_lart = dthldt_ls
    
    col = plt.cm.Greys((i+1)/len(plttime_var))
    ax.plot(thlv_av_time[plttime_var[i],:], qt_av_time[plttime_var[i],:],
            label='t=%.0f hr'%time[plttime_var[i]],color=col,lw=lw,alpha=alpha)
    ax.scatter(thlv_av_time[plttime_var[i],imin],qt_av_time[plttime_var[i],imin],
               color=col,zorder=100,label=ilab)
    ax.scatter(thlv_av_time[plttime_var[i],iqlbase],qt_av_time[plttime_var[i],iqlbase],
               marker='s',color=col,zorder=100,label=qlbaselab)
    ax.scatter(thlv_av_time[plttime_var[i],iqltop],qt_av_time[plttime_var[i],iqltop],
               marker='^',color=col,zorder=100,label=qltoplab)
    ax.set_xlabel(r"$\overline{\theta_{lv}} [K]$")
    ax.set_ylabel(r"$\overline{q_t}$ [kg/kg]")
    
    heights_budg = [iqlbase, imin, iqltop-25]
    # heights_budg = np.arange(iqlbase, iqltop,10)
    if i == len(plttime_var)-1:
        for k in range(len(heights_budg)):

            tendlab = terms[0] if k==0 else None
            fluxlab = terms[1] if k==0 else None
            subslab = terms[2] if k==0 else None  
            moislab = terms[3] if k==0 else None
            heatlab = terms[4] if k==0 else None
                        
            ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*ddt_thlv_av[heights_budg[k]]],
                    [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*ddt_qt_av[heights_budg[k]]],
                    color=colors[0],label=tendlab,lw=lw,alpha=alpha,linestyle='--')
            
            # ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*ddz_wthlv_av[heights_budg[k]]],
            #         [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*ddz_wqt_av[heights_budg[k]]],
            #         color=colors[1],label=fluxlab,lw=lw,alpha=alpha)

            # ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_subs[heights_budg[k]]],
            #         [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*qtavp_subs[heights_budg[k]]],
            #         color=colors[2],label=subslab,lw=lw,alpha=alpha)
            
            # ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_larq[heights_budg[k]]],
            #         [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]+fac*qtavp_larq[heights_budg[k]]],
            #         color=colors[3],label=moislab,lw=lw,alpha=alpha)
            
            # ax.plot([thlv_av_time[plttime_var[i],heights_budg[k]],thlv_av_time[plttime_var[i],heights_budg[k]]+fac*thlvavp_lart[heights_budg[k]]],
            #         [qt_av_time[plttime_var[i],heights_budg[k]],qt_av_time[plttime_var[i],heights_budg[k]]],
            #         color=colors[4],label=heatlab,lw=lw,alpha=alpha)
            
ax.legend(loc='upper left',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

plt.savefig(sp+'/conserved.pdf',bbox_inches='tight')

    
#%% Profiles of mean flux divergence

tpltmin = 6.
tpltmax = 18.
dit = 2 # Rounds to closest multiple of dt in time
tav = 1.0
alpha=0.5

itpltmin = np.where(time[plttime]>=tpltmin)[0][0]
itpltmax = np.where(time[plttime]<tpltmax)[0][-1]+1
idtplt = int(round(dit/(time[plttime[1]]-time[plttime[0]])))
plttime_var = np.arange(itpltmin,itpltmax,idtplt)

fig,axs = plt.subplots(ncols=4,sharey=True,figsize=(10,5),squeeze=False)
for i in range(len(plttime_var)):
    col = plt.cm.Greys((i+1)/len(plttime_var))
    
    itmin = np.argmin(np.abs(time1d - (time1d[plttime_var[i]]-tav*3600)))
    itmax = np.argmin(np.abs(time1d - (time1d[plttime_var[i]]+tav*3600)))
    
    # State
    thlv_av = np.mean(thlv_av_time[itmin:itmax],axis=0)
    qt_av = np.mean(qt_av_time[itmin:itmax],axis=0)
    
    # Flux divergence
    ddz_wthlv_av = np.mean(ddz_wthlv_av_time[itmin:itmax],axis=0)
    ddz_wqt_av = np.mean(ddz_wqt_av_time[itmin:itmax],axis=0)
         
    axs[0,0].plot(thlv_av, zflim, color=col,linestyle='-')
    axs[0,0].set_xlabel(r"$\overline{\theta_{lv}}$")
    # axs[0,0].set_xlim((0,6e-4))
    # axs[0,0].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[0,1].plot(qt_av, zflim, color=col,linestyle='-')
    axs[0,1].set_xlabel(r"$\overline{q_t}$")
    # axs[0,0].set_xlim((0,6e-4))
    # axs[0,1].ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    axs[0,2].plot(ddz_wthlv_av, zflim[1:-1], color=col,linestyle='-')
    axs[0,2].axvline(0,color='gray',linestyle='dotted')
    axs[0,2].set_xlabel(r"$\frac{\partial}{\partial z}\left(\overline{w'\theta_{lv}'}\right)$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,2].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    axs[0,3].plot(ddz_wqt_av, zflim[1:-1], color=col,linestyle='-',label='t=%.2f'%(time[plttime_var[i]]))
    axs[0,3].axvline(0,color='gray',linestyle='dotted')
    axs[0,3].set_xlabel(r"$\frac{\partial}{\partial z}\left(\overline{w'q_t'}\right)$")
    # axs[0,0].set_xlim((0,6e-4))
    axs[0,3].ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    
    if i == len(plttime_var) - 1:
        thlvz = np.diff(thlv_av)
        qtz = np.diff(qt_av)
        gamrat_z = np.diff(qtz/thlvz)
        gamrat_z[zflim[1:-1]<500] = 0
        gamrat_z[zflim[1:-1]>2500] = 0
        
        wthlvzz = np.diff(ddz_wthlv_av)
        wthlvzzz = np.diff(wthlvzz)
        wthlvzzz[zflim[2:-2]>2500] = 0
        wthlvzzz[zflim[2:-2]<500] = 0
        
        wqtzz = np.diff(ddz_wqt_av)
        wqtzzz = np.diff(wqtzz)
        
        axs[0,0].fill_betweenx(zflim[1:-1], thlv_av.min(), thlv_av.max(),
                               where=(gamrat_z>0),
                               facecolor='none',
                               edgecolor='silver',
                               hatch=r"\ ",
                               linewidth=0.0)
        axs[0,0].annotate('a)', (0.1,0.92), xycoords='axes fraction', fontsize=14)

        axs[0,1].fill_betweenx(zflim[1:-1], qt_av.min(), qt_av.max(),
                               where=(gamrat_z>0),
                               facecolor='none',
                               edgecolor='silver',
                               hatch=r"/",
                               linewidth=0.0)
        axs[0,1].annotate('b)', (0.1,0.92), xycoords='axes fraction', fontsize=14)

        axs[0,2].fill_betweenx(zflim[2:-2], ddz_wthlv_av.min(), ddz_wthlv_av.max(),
                               where=(wthlvzzz<-1e-9),
                               facecolor='none',
                               edgecolor='silver',
                               hatch=r"\ ",
                               linewidth=0.0)
        axs[0,2].annotate('c)', (0.1,0.92), xycoords='axes fraction', fontsize=14)

        axs[0,3].fill_betweenx(zflim[2:-2], ddz_wqt_av.min(), ddz_wqt_av.max(),
                               where=(wqtzzz<-2e-10),
                               facecolor='none',
                               edgecolor='silver',
                               hatch=r"/",
                               linewidth=0.0)
        axs[0,3].annotate('d)', (0.1,0.92), xycoords='axes fraction', fontsize=14)
        

axs[0,0].set_ylabel('z [m]')
axs[0,3].legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime_var)//13+1)

plt.savefig(sp+'/convexity.pdf',bbox_inches='tight')

#%% Time-averaged, slab-averaged budget

tpltmin = 6.
tpltmax = 10.

terms = ['Tendency',
         'Vertical flux convergence',
         'Subsidence',
         'Large-scale forcing',
         'residual'
         ]

colors = ['black',
          'lightsteelblue',
          'sienna',
          'cornflowerblue',
          'lightgray']

itpltmin = np.where(time1d/3600>=tpltmin)[0][0]
itpltmax = np.where(time1d/3600<tpltmax)[0][-1]+1

qtav_tend = np.mean(ddt_qt_av_time[itpltmin:itpltmax,:],axis=0)
qtav_vdiv = np.mean(ddz_wqt_av_time[itpltmin:itpltmax,:],axis=0)
qtav_subs = np.mean(wfls_dqtdz_av_time[itpltmin:itpltmax,:],axis=0)
qtav_larg = dqdt_ls[1:-1]
qtav_budg = -qtav_vdiv - qtav_subs + qtav_larg
qtav_resi = qtav_tend - qtav_budg

# The residual is mostly due to integration error of vertical transport
# -> Include the residual in this term
qtav_vdiv = qtav_vdiv - qtav_resi

thlvav_tend = np.mean(ddt_thlv_av_time[itpltmin:itpltmax,:],axis=0)
thlvav_vdiv = np.mean(ddz_wthlv_av_time[itpltmin:itpltmax,:],axis=0)
thlvav_subs = np.mean(wfls_dthlvdz_av_time[itpltmin:itpltmax,:],axis=0)
thlvav_larg = np.mean(dthlvdt_ls[itpltmin:itpltmax,1:-1],axis=0)
thlvav_budg = -thlvav_vdiv - thlvav_subs + thlvav_larg
thlvav_resi = thlvav_tend - thlvav_budg

thlvav_vdiv = thlvav_vdiv - thlvav_resi

alpha = 0.75
lw = 2

fig,axs = plt.subplots(ncols=2,sharey=True,figsize=(10,5))
# fig.suptitle(colors)
axs[0].plot( qtav_tend, zflim[1:-1],c=colors[0],alpha=alpha,lw=lw)
axs[0].plot(-qtav_vdiv, zflim[1:-1],c=colors[1],alpha=alpha,lw=lw)
axs[0].plot(-qtav_subs, zflim[1:-1],c=colors[2],alpha=alpha,lw=lw)
axs[0].plot( qtav_larg, zflim[1:-1],c=colors[3],alpha=alpha,lw=lw)
# axs[0].plot( qtav_resi, zflim[1:-1],c=colors[-1],alpha=alpha,lw=lw)
axs[0].set_xlabel(r"Contribution to $\overline{q_t}$ tendency [kg/kg/s]")

axs[1].plot( thlvav_tend, zflim[1:-1],c=colors[0],label=terms[0],alpha=alpha,lw=lw)
axs[1].plot(-thlvav_vdiv, zflim[1:-1],c=colors[1],label=terms[1],alpha=alpha,lw=lw)
axs[1].plot(-thlvav_subs, zflim[1:-1],c=colors[2],label=terms[2],alpha=alpha,lw=lw)
axs[1].plot( thlvav_larg, zflim[1:-1],c=colors[3],label=terms[3],alpha=alpha,lw=lw)
# axs[1].plot( thlvav_resi, zflim[1:-1],c=colors[-1],label='Residual')
axs[1].set_xlabel(r"Contribution to $\overline{\theta_{lv}}$ tendency [K/s]")
# axs[1].set_xlim((-7.5e-8,7.5e-8))
# axs[1].set_title('Dry')

axs[0].set_ylabel(r'Height [m]')
axs[1].legend(loc='best',bbox_to_anchor=(1,1))

plt.savefig(sp+'/slab_av_budget.pdf',bbox_inches='tight')
