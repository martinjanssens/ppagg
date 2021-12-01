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


lps = ['/scratch-shared/janssens/bomex100_e12/ppagg_merged',
        '/scratch-shared/janssens/bomex200_from100_12hr/ppagg_merged']
labs = [r'$\Delta x = 100m$',
        r'$\Delta x = 200m$']

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
    
    ds = nc.Dataset(lp+'/../fielddump.001.nc')
    ds1= nc.Dataset(lp+'/../profiles.001.nc')
    ds0= nc.Dataset(lp+'/../tmser.001.nc')
    ilp = np.loadtxt(lp+'/../lscale.inp.001')
    
    # time  = np.ma.getdata(ds.variables['time'][:]) / 3600
    ld[i]['time'] = np.load(lp+'/time.npy')
    ld[i]['zf']    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
    
    ld[i]['time1d'] = np.ma.getdata(ds1.variables['time'][:])
    ld[i]['rhobf'] = np.ma.getdata(ds1.variables['rhobf'][:])
    
    dzh = np.diff(ld[i]['zf'])[0] # FIXME only valid in lower part of domain
    
    # Larger-scale subsidence
    ld[i]['wfls'] = ilp[:,3]
    
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

#%% Plot variables

# Minus 'moist_time' or 'dry_time'
pltvars = ['qtpf','wthlvpf_anom', 'wff']
varlab = [r"$\widetilde{q_t'}$", 
          r"$\widetilde{w'\theta_{lv}'}-\overline{w'\theta_{lv}'}$", 
          r"$\widetilde{w'}$",]

pltvars = ['qtpf_prod_wex','qtpf_vdiv', 'qtpf_hdiv']
varlab = [r"Gradient production", 
          r"Vertical transport",
          r"Horizontal transport"]

# pltvars = ['thlvpf_prod','thlvpf_vdiv', 'thlvpf_hdiv']
# varlab = [r"Gradient production", 
#           r"Vertical transport",
#           r"Horizontal transport"]

pltvars = ['thlvpp']
varlab = [r"$\theta_{lv}'''$"]

lines = ['-','--']

tpltmin = 12.5
tpltmax = 14.5
dit = 1.0 # Rounds to closest multiple of dt in time
tav = 0.5 # Averaging time centred around current time
ndt = int((tpltmax-tpltmin)/dit)
nvar = len(pltvars)

fig,axs = plt.subplots(nrows=ndt,ncols=nvar,figsize=(4*nvar,3*ndt+0.25),
                       sharex=True,sharey=True,squeeze=False)

col_moist = plt.cm.RdYlBu(0.99)
col_dry = plt.cm.RdYlBu(0)
for l in range(len(lps)):

    itpltmin = np.where(ld[l]['time'][ld[l]['plttime']]>=tpltmin)[0][0]
    itpltmax = np.where(ld[l]['time'][ld[l]['plttime']]<tpltmax)[0][-1]+1
    idtplt = int(round(dit/(ld[l]['time'][ld[l]['plttime'][1]]-ld[l]['time'][ld[l]['plttime'][0]])))
    plttime_var = np.arange(itpltmin,itpltmax,idtplt)
    
    pltvars_moist = []
    pltvars_dry = []
    for p in range(nvar):
        pltvars_moist.append(ld[l][pltvars[p]+'_moist_time'])
        pltvars_dry.append(ld[l][pltvars[p]+'_dry_time'])
    
    for i in range(len(plttime_var)):
        
        ti = ld[l]['time'][plttime_var[i]]
        itav_min = np.where(ld[l]['time'][ld[l]['plttime']] >= ti-tav)[0][0]
        itav_max = np.where(ld[l]['time'][ld[l]['plttime']] <= ti+tav)[0][-1]
        
        for p in range(nvar):
            pltvar_moist_av = np.mean(pltvars_moist[p][itav_min:itav_max,:],axis=0)
            pltvar_dry_av = np.mean(pltvars_dry[p][itav_min:itav_max,:],axis=0)
            
            if len(ld[l]['zflim']) != len(pltvar_moist_av):
                zplt = ld[l]['zflim'][1:-1]
            else:
                zplt = ld[l]['zflim']
            
            axs[i,p].plot(pltvar_moist_av, zplt, color=col_moist, linestyle=lines[l], label=labs[l]+', moist')
            axs[i,p].plot(pltvar_dry_av, zplt, color=col_dry, linestyle=lines[l], label=labs[l]+', dry')

    for i in range(len(plttime_var)):
        axs[i,0].set_ylabel('Height [m]')
        for p in range(nvar):
            axs[i,p].set_title('%.1f hr'%ld[l]['time'][plttime_var[i]])

for p in range(nvar):
    axs[-1,p].set_xlabel(varlab[p])
axs[-1,-1].legend(loc='best',bbox_to_anchor=(1,-0.25),ncol=2)

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