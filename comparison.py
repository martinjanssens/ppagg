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
    
    dzh = np.diff(zf)[0] # FIXME only valid in lower part of domain
    
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
    ld[i]['qtpf_prod_moist_wex_time'] = np.load(lp+'/qtpf_prod_moist_wex_time.npy')
    ld[i]['qtpf_prod_dry_wex_time'] = np.load(lp+'/qtpf_prod_dry_wex_time.npy')
    ld[i]['qtpf_vdiv_moist_time'] = np.load(lp+'/qtpf_vdiv_moist_time.npy')
    ld[i]['qtpf_vdiv_dry_time'] = np.load(lp+'/qtpf_vdiv_dry_time.npy')
    ld[i]['qtpf_hdiv_moist_time'] = np.load(lp+'/qtpf_hdiv_moist_time.npy')
    ld[i]['qtpf_hdiv_dry_time'] = np.load(lp+'/qtpf_hdiv_dry_time.npy')
    ld[i]['qtpf_subs_moist_time'] = np.load(lp+'/qtpf_subs_moist_time.npy')
    ld[i]['qtpf_subs_dry_time'] = np.load(lp+'/qtpf_subs_dry_time.npy')
    
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

## Just lines to copy into a console with preloaded 100m res run for now
runcell(0, '/nfs/home4/janssens/scripts/pp3d/stats3d_eco_load.py')
runcell('Relation qtpf - wthlvpf_anom', '/nfs/home4/janssens/scripts/pp3d/stats3d_eco_load.py')
qtpf_dry_time100 = qtpf_dry_time
wthlvpf_moist_anom100 = wthlvpf_moist_anom
wthlvpf_dry_anom100 = wthlvpf_dry_anom
thlvpf_dry_time100 = thlvpf_dry_time
wff_dry_time100 = wff_dry_time
qlpf_dry_time100 = qlpf_dry_time
wthlvp_av_time100 = wthlvp_av_time
qtpf_moist_time100 = qtpf_moist_time
thlvpf_moist_time100 = thlvpf_moist_time
wff_moist_time100 = wff_moist_time
thlpf_moist_time100 = thlpf_moist_time
qlpf_moist_time100 = qlpf_moist_time
wqlpf_moist_time100 = wqlpf_moist_time
wqlpf_dry_time100 = wqlpf_moist_time
wthlpf_moist_time100 = wthlpf_moist_time
wthlpf_dry_time100 = wthlpf_dry_time
wthlvpf_moist_time100 = wthlvpf_moist_time
wthlvpf_dry_time100 = wthlvpf_dry_time
wqtpf_moist_time100 = wqtpf_moist_time
wqtpf_dry_time100 = wqtpf_moist_time
thl_av_time100 = thl_av_time

# old indices
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