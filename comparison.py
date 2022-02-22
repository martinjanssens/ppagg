#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:35:15 2021

@author: janssens
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy.optimize import curve_fit
from skimage.measure import block_reduce
from functions import tderive, zderivef, vint
from dataloader import DataLoaderDALES, DataLoaderMicroHH

# lps = ['/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/BOMEXStability/bomex100_e12/ppagg_new',
#        '/Users/martinjanssens/Documents/Wageningen/Patterns-in-satellite-images/BOMEXStability/bomex200_from100/ppagg_merged']
lps = ['/scratch-shared/janssens/bomex100_e12/ppagg',
        '/scratch-shared/janssens/bomex200_from100_12hr/ppagg_merged',
        '/scratch-shared/janssens/bomex100a5_from100_12hr/ppagg_merged',
        # '/scratch-shared/janssens/bomex200_fiso_from100_12hr/ppagg_merged',
        '/scratch-shared/janssens/bomex200_f200_from100_12hr/ppagg_merged',]
labs = [r'$\Delta x = 100m$',
        r'$\Delta x = 200m$',
        r'$\Delta x = 200m$, a5',
        r'$\Delta x = 200m$, f200']
mods = ['dales','dales','dales','dales']
# lps = ['/scratch-shared/janssens/bomex200_e12/ppagg',
       # '/scratch-shared/janssens/tmp.bomex/bomex_200m/ppagg',
       # '/scratch-shared/janssens/bomex100_e12/ppagg',
       # '/scratch-shared/janssens/tmp.bomex/bomex_100m/ppagg']
# labs = [r'DALES, $\Delta x = 200m$',
        # r'MicroHH, $\Delta x = 200m$',
        # r'DALES, $\Delta x = 100m$',
        # r'MicroHH, $\Delta x = 100m$']
# mods = ['dales','microhh','dales','microhh']
sp = lps[-1]+'/../figs'

# lps = ['/scratch-shared/janssens/bomex200aswitch/a2/ppagg',
#        '/scratch-shared/janssens/bomex200aswitch/a5_froma2_12hr/ppagg']
# labs = [r'2nd order advection',
#         r'5th order advection']

# lps = ['/scratch-shared/janssens/bomex200aswitch/a5/ppagg',
#        '/scratch-shared/janssens/bomex200aswitch/a2_froma5_12hr/ppagg']
# labs = [r'5th order advection',
#         r'2nd order advection']

# Loading loop
ld = []
for i in range(len(lps)):
    pp3d_out = {}
    ld.append(pp3d_out)
    lp = lps[i]
    mod = mods[i]

    if mod == 'dales':
        dl = DataLoaderDALES(lp+'/..')
    elif mod == 'microhh':
        dl = DataLoaderMicroHH(lp+'/..')
    
    # time  = np.ma.getdata(ds.variables['time'][:]) / 3600
    ld[i]['time'] = np.load(lp+'/time.npy')
    ld[i]['zf']   = dl.zf_inp
    
    ld[i]['time1d'] = dl.time1d
    ld[i]['rhobf'] = dl.rhobf
    
    dzh = np.diff(ld[i]['zf'])[0] # FIXME only valid in lower part of domain
    
    # Larger-scale subsidence
    ld[i]['wfls'] = dl.wfls
    
    ld[i]['plttime'] = np.load(lp+'/plttime.npy')
    ld[i]['zflim'] = np.load(lp+'/zf.npy')
    
    ld[i]['izmin'] = np.where(ld[i]['zflim'][0] == ld[i]['zf'])[0][0]
    ld[i]['izmax'] = np.where(ld[i]['zflim'][-1] == ld[i]['zf'])[0][0]+1
    
    ld[i]['qtpf_moist_time'] = np.load(lp+'/qtpf_moist_time.npy')
    ld[i]['qtpf_dry_time'] = np.load(lp+'/qtpf_dry_time.npy')
    ld[i]['qtpf_prod_moist_time'] = np.load(lp+'/qtpf_prod_moist_time.npy')
    ld[i]['qtpf_prod_dry_time'] = np.load(lp+'/qtpf_prod_dry_time.npy')
    ld[i]['qtpf_prod_wex_moist_time'] = np.load(lp+'/qtpf_prod_moist_wex_time.npy')
    ld[i]['qtpf_prod_wex_dry_time'] = np.load(lp+'/qtpf_prod_dry_wex_time.npy')
    ld[i]['qtpf_vdiv_moist_time'] = np.load(lp+'/qtpf_vdiv_moist_time.npy')
    ld[i]['qtpf_vdiv_dry_time'] = np.load(lp+'/qtpf_vdiv_dry_time.npy')
    ld[i]['qtpf_hdiv_moist_time'] = np.load(lp+'/qtpf_hdiv_moist_time.npy')
    ld[i]['qtpf_hdiv_dry_time'] = np.load(lp+'/qtpf_hdiv_dry_time.npy')
    ld[i]['qtpf_subs_moist_time'] = np.load(lp+'/qtpf_subs_moist_time.npy')
    ld[i]['qtpf_subs_dry_time'] = np.load(lp+'/qtpf_subs_dry_time.npy')
    
    ld[i]['thlvpf_moist_time'] = np.load(lp+'/thlvpf_moist_time.npy')
    ld[i]['thlvpf_dry_time'] = np.load(lp+'/thlvpf_dry_time.npy')
    ld[i]['thlvpf_prod_moist_time'] = np.load(lp+'/thlvpf_prod_moist_time.npy')
    ld[i]['thlvpf_prod_dry_time'] = np.load(lp+'/thlvpf_prod_dry_time.npy')
    ld[i]['thlvpf_vdiv_moist_time'] = np.load(lp+'/thlvpf_vdiv_moist_time.npy')
    ld[i]['thlvpf_vdiv_dry_time'] = np.load(lp+'/thlvpf_vdiv_dry_time.npy')
    ld[i]['thlvpf_hdiv_moist_time'] = np.load(lp+'/thlvpf_hdiv_moist_time.npy')
    ld[i]['thlvpf_hdiv_dry_time'] = np.load(lp+'/thlvpf_hdiv_dry_time.npy')
    ld[i]['thlvpf_subs_moist_time'] = np.load(lp+'/thlvpf_subs_moist_time.npy')
    ld[i]['thlvpf_subs_dry_time'] = np.load(lp+'/thlvpf_subs_dry_time.npy')
    
    ld[i]['thlvpp_moist_time'] = np.load(lp+'/thlvpp_moist_time.npy')
    ld[i]['thlvpp_dry_time'] = np.load(lp+'/thlvpp_dry_time.npy')
    ld[i]['thlvpp_prod_moist_time'] = np.load(lp+'/thlvpp_prod_moist_time.npy')
    ld[i]['thlvpp_prod_dry_time'] = np.load(lp+'/thlvpp_prod_dry_time.npy')
    ld[i]['thlvpp_vdiv_moist_time'] = np.load(lp+'/thlvpp_vdiv_moist_time.npy')
    ld[i]['thlvpp_vdiv_dry_time'] = np.load(lp+'/thlvpp_vdiv_dry_time.npy')
    ld[i]['thlvpp_hdiv_moist_time'] = np.load(lp+'/thlvpp_hdiv_moist_time.npy')
    ld[i]['thlvpp_hdiv_dry_time'] = np.load(lp+'/thlvpp_hdiv_dry_time.npy')
    ld[i]['thlvpp_subs_moist_time'] = np.load(lp+'/thlvpp_subs_moist_time.npy')
    ld[i]['thlvpp_subs_dry_time'] = np.load(lp+'/thlvpp_subs_dry_time.npy')
    
    ld[i]['wthlvpf_prod_moist_time'] = np.load(lp+'/wthlvpf_prod_moist_time.npy')
    ld[i]['wthlvpf_prod_dry_time'] =  np.load(lp+'/wthlvpf_prod_dry_time.npy')
    ld[i]['wthlvpf_vdiv_moist_time'] =  np.load(lp+'/wthlvpf_vdiv_moist_time.npy')
    ld[i]['wthlvpf_vdiv_dry_time'] = np.load(lp+'/wthlvpf_vdiv_dry_time.npy')
    ld[i]['wthlvpf_hdiv_moist_time'] = np.load(lp+'/wthlvpf_hdiv_moist_time.npy')
    ld[i]['wthlvpf_hdiv_dry_time'] = np.load(lp+'/wthlvpf_hdiv_dry_time.npy')
    ld[i]['wthlvpf_buoy_moist_time'] = np.load(lp+'/wthlvpf_buoy_moist_time.npy')
    ld[i]['wthlvpf_buoy_dry_time'] = np.load(lp+'/wthlvpf_buoy_dry_time.npy')
    ld[i]['wthlvpf_pres_moist_time'] = np.load(lp+'/wthlvpf_pres_moist_time.npy')
    ld[i]['wthlvpf_pres_dry_time'] = np.load(lp+'/wthlvpf_pres_dry_time.npy')
    ld[i]['wthlvpf_subs_moist_time'] = np.load(lp+'/wthlvpf_subs_moist_time.npy')
    ld[i]['wthlvpf_subs_dry_time'] = np.load(lp+'/wthlvpf_subs_dry_time.npy')
    ld[i]['wthlvpf_diff_moist_time'] = np.load(lp+'/wthlvpf_diff_moist_time.npy')
    ld[i]['wthlvpf_diff_dry_time'] = np.load(lp+'/wthlvpf_diff_dry_time.npy')
    
    ld[i]['thl_av_time'] = np.load(lp+'/thl_av_time.npy')
    ld[i]['thlv_av_time'] = np.load(lp+'/thlv_av_time.npy')
    ld[i]['qt_av_time'] = np.load(lp+'/qt_av_time.npy')
    
    ld[i]['thlpf_moist_time'] = np.load(lp+'/thlpf_moist_time.npy')
    ld[i]['thlpf_dry_time'] = np.load(lp+'/thlpf_dry_time.npy')
    ld[i]['wff_moist_time'] = np.load(lp+'/wff_moist_time.npy')
    ld[i]['wff_dry_time'] = np.load(lp+'/wff_dry_time.npy')
    ld[i]['qlpf_moist_time'] = np.load(lp+'/qlpf_moist_time.npy') 
    ld[i]['qlpf_dry_time'] = np.load(lp+'/qlpf_dry_time.npy')
    
    ld[i]['thlpp_moist_time'] = np.load(lp+'/thlpp_moist_time.npy')
    ld[i]['thlpp_dry_time'] = np.load(lp+'/thlpp_dry_time.npy')
    ld[i]['wfp_moist_time'] = np.load(lp+'/wfp_moist_time.npy')
    ld[i]['wfp_dry_time'] = np.load(lp+'/wfp_dry_time.npy')
    ld[i]['qlpp_moist_time'] = np.load(lp+'/qlpp_moist_time.npy') 
    ld[i]['qlpp_dry_time'] = np.load(lp+'/qlpp_dry_time.npy')
    
    ld[i]['wthlpf_moist_time'] = np.load(lp+'/wthlpf_moist_time.npy')
    ld[i]['wthlpf_dry_time'] = np.load(lp+'/wthlpf_dry_time.npy')
    
    ld[i]['wqtpf_moist_time'] = np.load(lp+'/wqtpf_moist_time.npy')
    ld[i]['wqtpf_dry_time'] = np.load(lp+'/wqtpf_dry_time.npy')
    
    ld[i]['wqlpf_moist_time'] = np.load(lp+'/wqlpf_moist_time.npy')
    ld[i]['wqlpf_dry_time'] = np.load(lp+'/wqlpf_dry_time.npy')
    
    ld[i]['wthlvp_av_time'] = np.load(lp+'/wthlvp_av_time.npy')
    ld[i]['wthlvpf_moist_time'] = np.load(lp+'/wthlvpf_moist_time.npy')
    ld[i]['wthlvpf_dry_time'] = np.load(lp+'/wthlvpf_dry_time.npy')
    ld[i]['wthlvpf_l_moist_time'] = np.load(lp+'/wthlvpf_l_moist_time.npy')
    ld[i]['wthlvpf_l_dry_time'] = np.load(lp+'/wthlvpf_l_dry_time.npy')
    ld[i]['wthlvpf_c_moist_time'] = np.load(lp+'/wthlvpf_c_moist_time.npy')
    ld[i]['wthlvpf_c_dry_time'] = np.load(lp+'/wthlvpf_c_dry_time.npy')
    ld[i]['wthlvpf_r_moist_time'] = np.load(lp+'/wthlvpf_r_moist_time.npy')
    ld[i]['wthlvpf_r_dry_time'] = np.load(lp+'/wthlvpf_r_dry_time.npy')
    ld[i]['wthlvpp_moist_time'] = np.load(lp+'/wthlvpp_moist_time.npy')
    ld[i]['wthlvpp_dry_time'] = np.load(lp+'/wthlvpp_dry_time.npy')
    
    ld[i]['wthlvpf_anom_moist_time'] = ld[i]['wthlvpf_moist_time'] - ld[i]['wthlvp_av_time']
    ld[i]['wthlvpf_anom_dry_time'] = ld[i]['wthlvpf_dry_time'] - ld[i]['wthlvp_av_time']

    ld[i]['Gamma_qt_av_time'] = zderivef(ld[i]['qt_av_time'],dzh)
    ld[i]['Gamma_thlv_av_time'] = zderivef(ld[i]['thlv_av_time'],dzh)
    ld[i]['Gamrat_av_time'] = ld[i]['Gamma_qt_av_time']/ld[i]['Gamma_thlv_av_time']
    ld[i]['Gamrat_av_time'][np.abs(ld[i]['Gamrat_av_time'])>0.03] = np.nan
    ## Reconstruct slab-mean budget terms
    ## FIXME Not working yet, would need support for time1d vs time and handling different time dimension sizes in restart and original
    # thl_av_1d = ds1['thl'][:,ld[i]['izmin']:ld[i]['izmax']]
    # qt_av_1d = ds1['qt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # thlv_av_1d = thl_av_1d*(1 + 0.608*qt_av_1d)
    
    # # Tendencies
    # ld[i]['ddt_thlv_av_time'] = tderive(thlv_av_1d, ld[i]['time1d']/3600)
    # ld[i]['ddt_qt_av_time'] = tderive(qt_av_1d, ld[i]['time1d']/3600)
    
    # # Flux divergence (approximately, i.e. ignoring rho)
    # ld[i]['wthl_av_time'] = ds1['wthlt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # ld[i]['wqt_av_time'] = ds1['wqtt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # ld[i]['wthlv_av_time'] = ld[i]['wthl_av_time'] + 0.608*thl_av_1d*ld[i]['wqt_av_time']
    
    # ld[i]['ddz_wthlv_av_time'] = zderivef(ld[i]['wthlv_av_time'],dzh)
    # ld[i]['ddz_wqt_av_time'] = zderivef(ld[i]['wqt_av_time'],dzh)

#%% Function for plotting chosen variable comparison

def plot_comparison(ld,pltvars,varlab,tpltmin,tpltmax,dit,tav,lines,
                    sharex=True,alpha=0.75,lw=2):
    ndt = int((tpltmax-tpltmin)/dit)
    nvar = len(pltvars)

    fig,axs = plt.subplots(nrows=ndt,ncols=nvar,figsize=(2.5*nvar,3*ndt+0.25),
                           sharex=sharex,sharey=True,squeeze=False)
    col_av = 'k'
    col_moist = plt.cm.RdYlBu(0.99)
    col_dry = plt.cm.RdYlBu(0)
    lns= []; lbs = []
    for l in range(len(lps)):
    
        itpltmin = np.where(ld[l]['time'][ld[l]['plttime']]>=tpltmin)[0][0]
        itpltmax = np.where(ld[l]['time'][ld[l]['plttime']]<tpltmax)[0][-1]+1
        idtplt = int(round(dit/(ld[l]['time'][ld[l]['plttime'][1]]-ld[l]['time'][ld[l]['plttime'][0]])))
        plttime_var = np.arange(itpltmin,itpltmax,idtplt)
        
        pltvars_moist = []
        pltvars_dry = []
        pltvars_av = []
        for p in range(nvar):
            if '_av' in pltvars[p]:
                pltvars_av.append(ld[l][pltvars[p]+'_time'])
            else:
                pltvars_moist.append(ld[l][pltvars[p]+'_moist_time'])
                pltvars_dry.append(ld[l][pltvars[p]+'_dry_time'])
        
        for i in range(len(plttime_var)):
            
            ti = ld[l]['time'][plttime_var[i]]
            itmn_min = np.where(ld[l]['time'][ld[l]['plttime']] >= ti-tav)[0][0]
            itmn_max = np.where(ld[l]['time'][ld[l]['plttime']] <= ti+tav)[0][-1]
            
            for p in range(nvar):
                if '_av' in pltvars[p]:
                    
                    pltvar_av_mn = np.mean(ld[l][pltvars[p]+'_time'][itmn_min:itmn_max,:],axis=0)
                    
                    diffz = len(ld[l]['zflim']) - len(pltvar_av_mn)
                    if diffz != 0:
                        zplt = ld[l]['zflim'][diffz//2:-diffz//2]
                    else:
                        zplt = ld[l]['zflim']
                    
                    lab = labs[l]+', slab-averaged'
                    ln = axs[i,p].plot(pltvar_av_mn, zplt, 
                                  color=col_av, linestyle=lines[l],
                                  alpha=alpha, lw=lw)
                    if i == 0:
                        lns.append(ln[0])
                        lbs.append(lab)
                    # axs[i,p].set_xlim((-0.05,0.05))
                else:
                    pltvar_moist_mn = np.mean(ld[l][pltvars[p]+'_moist_time'][itmn_min:itmn_max,:],axis=0)
                    pltvar_dry_mn = np.mean(ld[l][pltvars[p]+'_dry_time'][itmn_min:itmn_max,:],axis=0)
                    
                    diffz = len(ld[l]['zflim']) - len(pltvar_moist_mn)
                    if diffz != 0:
                        zplt = ld[l]['zflim'][diffz//2:-diffz//2]
                    else:
                        zplt = ld[l]['zflim']
                    
                    # Ugly but works if you start with moist/dry var:
                    labm = labs[l]+', moist'
                    labd = labs[l]+', dry'
                    lnm = axs[i,p].plot(pltvar_moist_mn, zplt, 
                                  color=col_moist, linestyle=lines[l],
                                  alpha=alpha, lw=lw)
                    lnd = axs[i,p].plot(pltvar_dry_mn, zplt, 
                                  color=col_dry, linestyle=lines[l], 
                                  alpha=alpha, lw=lw)
                    if i == 0 and p == 0:
                        lns.append(lnm[0])
                        lns.append(lnd[0])
                        lbs.append(labm)
                        lbs.append(labd)
    
        for i in range(len(plttime_var)):
            axs[i,0].set_ylabel('Height [m]')
            for p in range(nvar):
                axs[i,p].set_title('%.0f'%(ld[l]['time'][plttime_var[i]]-tav)+'-'+
                                   '%.0f hr'%(ld[l]['time'][plttime_var[i]]+tav))
    
    for p in range(nvar):
        axs[-1,p].set_xlabel(varlab[p])
    # axs[-1,-1].legend(loc='best',bbox_to_anchor=(1,-0.25),ncol=2)
    fig.legend(lns, lbs, bbox_to_anchor=(0.9,0.075),ncol=len(lps))

#%% Plot variables

# Minus 'moist_time' or 'dry_time'
pltvars = ['qtpf','wthlvpf_anom', 'Gamrat_av','qtpf_prod_wex']
varlab = [r"${q_{t_m}'}$ [kg/kg]", 
          r"$F_{{\theta_{lv}'}_m}$ [K m/s]", 
          r"$\Gamma_{q_t}/\Gamma_{\theta_{lv}}$ [kg/kg/K]",
          r"$w_m'\Gamma_{q_t}$ [kg/kg/s]"]

# pltvars = ['qtpf_prod_wex','qtpf_vdiv', 'qtpf_hdiv']
# varlab = [r"Gradient production", 
#           r"Vertical transport",
#           r"Horizontal transport"]

# pltvars = ['thlvpf_prod','thlvpf_vdiv', 'thlvpf_hdiv']
# varlab = [r"Gradient production", 
#           r"Vertical transport",
#           r"Horizontal transport"]

# pltvars = ['thlvpp']
# varlab = [r"$\theta_{lv}'''$"]

lines = ['-','--',':','-.']

tpltmin = 12
tpltmax = 20
dit = 2.0 # Rounds to closest multiple of dt in time
tav = 1.0 # Averaging time centred around current time


plot_comparison(ld,pltvars,varlab,tpltmin,tpltmax,dit,tav,lines,sharex='col')
plt.savefig(sp+'/comparison_vars.pdf',bbox_inches='tight')

#%% qtpf important budget terms
pltvars = ['qtpf_prod_wex','qtpf_vdiv', 'qtpf_hdiv']
varlab = [r"Gradient production", 
          r"Vertical transport",
          r"Horizontal transport"]

# pltvars = ['thlvpf_prod','thlvpf_vdiv', 'thlvpf_hdiv']
# varlab = [r"Gradient production", 
#           r"Vertical transport",
#           r"Horizontal transport"]

lines = ['-','--',':','-.']

tpltmin = 6
tpltmax = 24
dit = 6.0 # Rounds to closest multiple of dt in time
tav = 4.0 # Averaging time centred around current time

plot_comparison(ld,pltvars,varlab,tpltmin,tpltmax,dit,tav,lines)
plt.savefig(sp+'/comparison_qtpf.pdf',bbox_inches='tight')

#%% Fluxes
pltvars = ['wthlvpf_r','wqlpf','wqtpf']
varlab = [r"$\left(w_s'\theta_{lv_s}'\right)$",
          r"$\left(w_s'q_l'\right)$",
          r"$\left(w_s'q_t'\right)$",]

lines = ['-','--',':','-.']

tpltmin = 13
tpltmax = 19
dit = 2.0 # Rounds to closest multiple of dt in time
tav = 1.0 # Averaging time centred around current time

plot_comparison(ld,pltvars,varlab,tpltmin,tpltmax,dit,tav,lines,sharex='col')

#%% wthlv budget
pltvars = ['wthlvpf_prod','wthlvpf_vdiv','wthlvpf_diff']
varlab =  [r"Gradient production", 
          r"Vertical transport",
          r"Horizontal transport"]
lines = ['-','--',':','-.']

tpltmin = 13
tpltmax = 19
dit = 2.0 # Rounds to closest multiple of dt in time
tav = 1.0 # Averaging time centred around current time

plot_comparison(ld,pltvars,varlab,tpltmin,tpltmax,dit,tav,lines,sharex='col')


#%% Plot 1d comparison of vertically integrated mesoscale moisture fluctuation

ls = ['-','--',':','-.']

tmin = 6.
tmax = [18.,
        12.,
        36.,
        36.,
        36.]

alpha=0.5
lw=2
col_moist = plt.cm.RdYlBu(0.99)
col_dry = plt.cm.RdYlBu(0)

f1 = plt.figure(figsize=(5,10/3)); axs1 = plt.gca()
for i in range(len(lps)):
    time = ld[i]['time']
    itpltmin = np.where(time>=tmin)[0][0]
    itpltmax = np.where(time<tmax[i])[0][-1]+1
    plttime_var = np.arange(itpltmin,itpltmax,1)
    
    zflim = ld[i]['zflim']
    rhobfi = ld[i]['rhobf'][0,ld[i]['izmin']:ld[i]['izmax']]
    qtpfmi = ld[i]['qtpf_moist_time']
    qtpfdi = ld[i]['qtpf_dry_time']

    twppf_moist = vint(qtpfmi,rhobfi,zflim,plttime_var)
    twppf_dry = vint(qtpfdi,rhobfi,zflim,plttime_var)

    axs1.plot(time[plttime_var],twppf_moist,c=col_moist,linestyle=ls[i],lw=lw,alpha=alpha,label=labs[i])
    axs1.plot(time[plttime_var],twppf_dry,c=col_dry,linestyle=ls[i],lw=lw,alpha=alpha)
axs1.set_xlabel('Time [hr]')
axs1.set_ylabel(r"$TWP_m'$ [kg/m$^2$]")
axs1.legend(loc='upper left',bbox_to_anchor=(1,1))
plt.savefig(sp+'/twp_evo_num.pdf',bbox_inches='tight')

#%%
plt.plot(np.mean(wff_moist_time[35:39,:],axis=0),zflim)
plt.plot(np.mean(wff_moist_time100[23:25,:],axis=0),zflim)
plt.plot(np.mean(wff_dry_time[35:39,:],axis=0),zflim)
plt.plot(np.mean(wff_dry_time100[23:25,:],axis=0),zflim)

plt.plot(np.mean(wthlvpf_moist_time[35:39,:],axis=0),zflim,c='C0')
plt.plot(np.mean(wthlvpf_moist_time100[23:25,:],axis=0),zflim,c='C0',linestyle='--')
plt.plot(np.mean(wthlpf_moist_time[35:39,:],axis=0),zflim,c='C1')
plt.plot(np.mean(wthlpf_moist_time100[23:25,:],axis=0),zflim,c='C1',linestyle='--')
plt.plot(np.mean(0.608*thl_av_time[35:39,:]*wqtpf_moist_time[35:39,:],axis=0),zflim,c='C1')
plt.plot(np.mean(0.608*thl_av_time100[23:25,:]*wqtpf_moist_time100[23:25,:],axis=0),zflim,c='C1',linestyle='--')

# 12 hour restart
its=47
ite=51
bsl = '100m'
opt = '200m'
plt.plot(np.mean(wff_moist_time[its:ite,:],axis=0),zflim,c='C0',label=bsl+', moist')
plt.plot(np.mean(wff_moist_time100[its:ite,:],axis=0),zflim,c='C0',linestyle='--',label=opt+', moist')
plt.plot(np.mean(wff_dry_time[its:ite,:],axis=0),zflim,c='C1',label=bsl+', dry')
plt.plot(np.mean(wff_dry_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--',label=opt+', dry')
plt.legend()
plt.show()

plt.plot(np.mean(wthlvpf_moist_time[its:ite,:]-wthlvp_av_time[its:ite,:],axis=0),zflim,c='C0',label=bsl+', moist')
plt.plot(np.mean(wthlvpf_moist_time100[its:ite,:]-wthlvp_av_time100[its:ite,:],axis=0),zflim,c='C0',linestyle='--',label=opt+', moist')
plt.plot(np.mean(wthlvpf_dry_time[its:ite,:]-wthlvp_av_time[its:ite,:],axis=0),zflim,c='C1',label=bsl+', dry')
plt.plot(np.mean(wthlvpf_dry_time100[its:ite,:]-wthlvp_av_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--',label=opt+', dry')
plt.legend()
plt.show()

plt.plot(np.mean(wthlvpf_moist_time[its:ite,:],axis=0),zflim,c='C0')
plt.plot(np.mean(wthlvpf_moist_time100[its:ite,:],axis=0),zflim,c='C0',linestyle='--')
plt.plot(np.mean(wthlpf_moist_time[its:ite,:],axis=0),zflim,c='C1')
plt.plot(np.mean(wthlpf_moist_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--')
plt.plot(np.mean(0.608*thl_av_time[its:ite,:]*wqtpf_moist_time[its:ite,:],axis=0),zflim,c='C1')
plt.plot(np.mean(0.608*thl_av_time100[its:ite,:]*wqtpf_moist_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--')


#%% Spectra
k1d100 = k1d
spec_qt_mn100 = spec_qt_mn
spec_thl_mn100 = spec_thl_mn
spec_thlv_mn100 = spec_thlv_mn
spec_w_mn100 = spec_w_mn
spec_wqt_mn100 = spec_wqt_mn
spec_wthl_mn100 = spec_wthl_mn
spec_wthlv_mn100 = spec_wthlv_mn

labbsl = r"100m"
labopt = r"200m from 100m"

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_qt_mn[-1,izpl,:],label=labbsl)
ax.loglog(k1d100,spec_qt_mn100[-1,izpl,:],label=labopt)
ax.set_ylabel(r"$k\widehat{q_t}'^2$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_thlv_mn[-1,izpl,:],label=labbsl)
ax.loglog(k1d100,spec_thlv_mn100[-1,izpl,:],label=labopt)
ax.set_ylabel(r"$k\widehat{\theta_{lv}}'^2$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_w_mn[-1,izpl,:],label=labbsl)
ax.loglog(k1d100,spec_w_mn100[-1,izpl,:],label=labopt)
ax.set_ylabel(r"$k\widehat{w}'^2$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()

fig = plt.figure(); ax = plt.gca()
ax.loglog(k1d,spec_wthlv_mn[-1,izpl,:],label=labbsl)
ax.loglog(k1d100,spec_wthlv_mn100[-1,izpl,:],label=labopt)
ax.set_ylabel(r"$k\widehat{w'\theta_{lv}'}^2$")
ax.set_xlabel(r"Wavenumber [1/m]")
ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

ax2 = ax.twiny()
fig.subplots_adjust(bottom=0.22)
ax2.xaxis.set_ticks_position('bottom')
ax2.xaxis.set_label_position('bottom')
ax2.spines['bottom'].set_position(('axes',-0.22))
ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
ax2.set_xscale('log')
ax2.set_xlabel('Wavelength [m]')
plt.show()
