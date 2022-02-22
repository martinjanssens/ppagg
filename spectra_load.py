#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 13:44:52 2021

@author: janssens
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from skimage.measure import block_reduce
import gc
from functions import *
import argparse
from matplotlib.patches import FancyArrowPatch
from scipy.optimize import curve_fit


lps = ['/scratch-shared/janssens/bomex100_e12/spectra',
        '/scratch-shared/janssens/bomex200_from100_12hr/spectra',
        # '/scratch-shared/janssens/bomex100a5_from100_12hr/spectra',
        '/scratch-shared/janssens/bomex200_fiso_from100_12hr/spectra',
        '/scratch-shared/janssens/bomex200_f200_from100_12hr/spectra']
labs = [r'$\Delta x = 100m$',
        r'$\Delta x = 200m$',
        # r'$\Delta x = 100m$, a5',
        r'$\Delta x = 200m$, fiso',
        r'$\Delta x = 200m$, f200']

# lps = ['/scratch-shared/janssens/bomex200_e12/spectra',
#        '/scratch-shared/janssens/tmp.bomex/bomex_200m/spectra',
#        '/scratch-shared/janssens/bomex100_e12/spectra',
#        '/scratch-shared/janssens/tmp.bomex/bomex_100m/spectra']
# labs = [r'DALES, $\Delta x = 200m$',
#         r'MicroHH, $\Delta x = 200m$',
#         r'DALES, $\Delta x = 100m$',
#         r'MicroHH, $\Delta x = 100m$']

def _add_twinx(fig, ax, offset=0.22):
    ax2 = ax.twiny()
    fig.subplots_adjust(bottom=offset)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_label_position('bottom')
    ax2.spines['bottom'].set_position(('axes',-offset))
    ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
    ax2.set_xscale('log')
    ax2.set_xlabel('Wavelength [m]')

def plot_spectrum(k1d, spec, lab, plttime, time):
    fig = plt.figure(); ax = plt.gca()
    for i in range(len(plttime)):
        col = plt.cm.cubehelix(i/len(plttime))
        ax.loglog(k1d,spec[i,izpl,:],c=col,label='t=%.2f'%time[plttime[i]])
    # ax.set_ylim((1e-4,1e2))
    ax.set_ylabel(lab)
    ax.set_xlabel(r"Wavenumber [1/m]")
    ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)
    
    _add_twinx(fig, ax)
    
# Loading loop
ld = []
sps = []
for i in range(len(lps)):
    spec_out = {}
    ld.append(spec_out)
    lp = lps[i]
    sps.append(lp+'/figs')
    
    ld[i]['time'] = np.load(lp+'/time_spec.npy')
    ld[i]['plttime'] = np.load(lp+'/plttime_spec.npy')
    ld[i]['zf'] = np.load(lp+'/zf_spec.npy')
    ld[i]['k1d'] = np.load(lp+'/k1d.npy')

    ld[i]['spec_qt'] = np.load(lp+'/spec_qt.npy')
    ld[i]['spec_thl'] = np.load(lp+'/spec_thl.npy')
    ld[i]['spec_thlv'] = np.load(lp+'/spec_thlv.npy')
    ld[i]['spec_w'] = np.load(lp+'/spec_w.npy')
    ld[i]['spec_ql'] = np.load(lp+'/spec_ql.npy')

    ld[i]['spec_wqt'] = np.load(lp+'/spec_wqt.npy')
    ld[i]['spec_wthl'] = np.load(lp+'/spec_wthl.npy')
    ld[i]['spec_wthlv'] = np.load(lp+'/spec_wthlv.npy')
    ld[i]['spec_wql'] = np.load(lp+'/spec_wql.npy')

#%% Plot time evolution of spectra per loaded run 

itav = 4 # number of time steps to average over -> len(plttime) MUST BE MULTIPLE OF THIS
izpl = 9 # (Closest) height to plot spectra at

for i in range(len(lps)):
    
    #Average over time
    spec_qt_mn = block_reduce(ld[i]['spec_qt'],(itav,1,1),func=np.mean)
    spec_thl_mn = block_reduce(ld[i]['spec_thl'],(itav,1,1),func=np.mean)
    spec_thlv_mn = block_reduce(ld[i]['spec_thlv'],(itav,1,1),func=np.mean)
    spec_w_mn = block_reduce(ld[i]['spec_w'],(itav,1,1),func=np.mean)
    spec_ql_mn = block_reduce(ld[i]['spec_ql'],(itav,1,1),func=np.mean)
    
    spec_wqt_mn = block_reduce(ld[i]['spec_wqt'],(itav,1,1),func=np.mean)
    spec_wthl_mn = block_reduce(ld[i]['spec_wthl'],(itav,1,1),func=np.mean)
    spec_wthlv_mn = block_reduce(ld[i]['spec_wthlv'],(itav,1,1),func=np.mean)
    spec_wql_mn = block_reduce(ld[i]['spec_wql'],(itav,1,1),func=np.mean)
    
    plttime_mn = ld[i]['plttime'][::itav]
    
    # Plotting
    # Variances
    plot_spectrum(ld[i]['k1d'], spec_qt_mn, r"$k\widehat{q}_t'^2$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_qt.pdf', bbox_inches='tight')

    plot_spectrum(ld[i]['k1d'], spec_thl_mn, r"$k\widehat{\theta}_l'^2$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_thl.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_thlv_mn, r"$k\widehat{\theta}_{lv}'^2$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_thlv.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_w_mn, r"$k\widehat{w}'^2$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_w.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_ql_mn, r"$k\widehat{q_l}'^2$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_ql.pdf', bbox_inches='tight')
    
    # Fluxes
    plot_spectrum(ld[i]['k1d'], spec_wqt_mn, r"$k\widehat{wq}_t'$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_wqt.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_wthl_mn, r"$k\widehat{w\theta}_l'$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_wthl.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_wthlv_mn, r"$k\widehat{w\theta}_{lv}'$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_wthlv.pdf', bbox_inches='tight')
    
    plot_spectrum(ld[i]['k1d'], spec_wql_mn, r"$k\widehat{wq}_{l}'$", plttime_mn, ld[i]['time'])
    plt.savefig(sps[i]+'/spec_wql.pdf', bbox_inches='tight')

#%% Plot multiple runs  at same time

tplt = 13.
tav = 1. # Hours after tplt
zplt = 1500.
klp = 4
kar = 255

pltvars = ['spec_qt',
           'spec_thlv',
           'spec_w']
varlab = [r"$k\widehat{q}_t'^2$",
          r"$k\widehat{\theta}_{lv}'^2$",
          r"$k\widehat{w}'^2$"]
lines = ['-','--',':','-.']
dashes=[(1,0),(3,6),(1,1),(2,2)]

iz = np.argmin(np.abs(ld[i]['zf'] - zplt))

fig,axs = plt.subplots(nrows=len(pltvars),ncols=1,figsize=(5,len(pltvars)*4+0.4),sharex=True)

lns = []; lbs = []
are = np.zeros(len(pltvars))
for i in range(len(lps)):
    
    itpltmin = np.argmin(np.abs(ld[i]['time'][ld[i]['plttime']] - tplt))
    itpltmax = np.argmin(np.abs(ld[i]['time'][ld[i]['plttime']] - (tplt + tav)))
    itplt = ld[i]['plttime'][itpltmin:itpltmax]
    
    k1d_plt = ld[i]['k1d']
    
    
    for j in range(len(pltvars)):
        
        spec_plt = np.mean(ld[i][pltvars[j]][itpltmin:itpltmax,iz], axis=0)
        
        ln = axs[j].loglog(k1d_plt[::2],spec_plt[::2],c='k',linestyle=lines[i], dashes=dashes[i])
        
        if i == 0:
            axs[j].set_ylabel(varlab[j])
            axs[j].axvline(k1d_plt[klp],c='Gray')
            
            b0 = np.max(spec_plt)*1e-3
            axs[j].loglog(k1d_plt[100:], b0 * k1d_plt[100:] ** (-5/3), c="Gray")
            if j == len(pltvars)-1:
                axs[j].set_xlabel(r"Wavenumber [1/m]")
                _add_twinx(fig, axs[j], offset=0.4)
                axs[j].annotate(r'$k_m$', (0.25,-0.15),xycoords='axes fraction')        
        if j == 0:
            lns.append(ln[0])
            lbs.append(labs[i])
        if i > 0:
            # axs[j].axvline(k1d_plt[kar],spec_plt[kar],are[j],c='Gray')
            # print(are[j])
            # print((k1d_plt[kar],spec_plt[kar]), (k1d_plt[kar],are[j]))
            ar = FancyArrowPatch(posA=(k1d_plt[kar],spec_plt[kar]), 
                                 posB=(k1d_plt[kar],are[j]),
                                 arrowstyle='<|-|>', color='0.5',
                                 mutation_scale=10,linestyle=lines[i])
            axs[j].add_artist(ar)
        are[j] = spec_plt[kar]
fig.legend(lns, lbs, bbox_to_anchor=(0.9,0.3),ncol=2)    
plt.savefig(sps[-1]+'/spectra_comparison.pdf', bbox_inches='tight')

#%% Plot in same spectrum FIXME is not yet implemented

# fig = plt.figure(); ax = plt.gca()
# ax.loglog(k1d,spec_wthlv_mn[-1,izpl,:],label=r"$w\theta_{lv}$")
# ax.loglog(k1d,spec_wthlv_t_mn[-1,izpl,:],label=r"$w\theta_{lv}'$")
# ax.loglog(k1d,spec_wthlv_r_mn[-1,izpl,:],label=r"$w'''\theta_{lv}'''$")
# # ax.set_ylim((1e-4,1e2))
# ax.set_ylabel(r"$k\widehat{w\theta}_{lv}'$")
# ax.set_xlabel(r"Wavenumber [1/m]")
# ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

# ax2 = ax.twiny()
# fig.subplots_adjust(bottom=0.22)
# ax2.xaxis.set_ticks_position('bottom')
# ax2.xaxis.set_label_position('bottom')
# ax2.spines['bottom'].set_position(('axes',-0.22))
# ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
# ax2.set_xscale('log')
# ax2.set_xlabel('Wavelength [m]')
# plt.show()

# #%% wthlv scale decomposition (run cells above first)
# klp = 4

# spec_wthlv_l = np.zeros((len(plttime),len(zflim),N2))
# spec_wthlv_c = np.zeros((len(plttime),len(zflim),N2))
# spec_wthlv_r = np.zeros((len(plttime),len(zflim),N2))

# # Mask for low-[ass filtering
# circ_mask = np.zeros((xf.size,xf.size))
# rad = getRad(circ_mask)
# circ_mask[rad<=klp] = 1

# for i in range(len(plttime)):
    
#     # 3D fields
#     qt  = np.ma.getdata(ds.variables['qt'][plttime[i],izmin:izmax,:,:])
#     wh = np.ma.getdata(ds.variables['w'][plttime[i],izmin:izmax+1,:,:])
#     thl =  np.ma.getdata(ds.variables['thl'][plttime[i],izmin:izmax,:,:])
    
#     thl_av = np.mean(thl,axis=(1,2))
#     thlv = thl + 0.608*thl_av[:,np.newaxis,np.newaxis]*qt
#     thlv = thlv - np.mean(thlv,axis=(1,2))
    
#     wf = (wh[1:,:,:] + wh[:-1,:,:])*0.5
#     del wh

#     # Low-pass filter (and identify high-pass filtered remainder)    
#     wff = lowPass(wf, circ_mask)
#     wfp = wf - wff
#     del wf
                
#     thlvf = lowPass(thlv, circ_mask)
#     thlvp = thlv - thlvf
#     del thlv
    
#     gc.collect()
    
#     for iz in range(len(zflim)):
#         # k1d,spec_wthlv_l[i,iz,:] = compute_spectrum(wff[iz,:,:]*thlvf[iz,:,:], dx)
#         # k1d,spec_wthlv_c[i,iz,:] = compute_spectrum(wff[iz,:,:]*thlvp[iz,:,:]+
#         #                                             wfp[iz,:,:]*thlvf[iz,:,:], dx)
#         # k1d,spec_wthlv_r[i,iz,:] = compute_spectrum(wfp[iz,:,:]*thlvp[iz,:,:], dx)

#         k1d,spec_wthlv_l[i,iz,:] = compute_spectrum(wff[iz,:,:], dx,thlvf[iz,:,:])
#         k1d,spec_wthlv_c[i,iz,:] = compute_spectrum(wff[iz,:,:], dx, thlvp[iz,:,:])
#         _,spec_wthlv_c2 = compute_spectrum(wfp[iz,:,:], dx, thlvf[iz,:,:])
#         spec_wthlv_c[i,iz,:] += spec_wthlv_c2
#         k1d,spec_wthlv_r[i,iz,:] = compute_spectrum(wfp[iz,:,:], dx, thlvp[iz,:,:])


# spec_wthlv_l_mn = block_reduce(spec_wthlv_l,(itav,1,1),func=np.mean)
# spec_wthlv_c_mn = block_reduce(spec_wthlv_c,(itav,1,1),func=np.mean)
# spec_wthlv_r_mn = block_reduce(spec_wthlv_r,(itav,1,1),func=np.mean)
# sumtest = spec_wthlv_l_mn + spec_wthlv_c_mn + spec_wthlv_r_mn

# fig = plt.figure(); ax = plt.gca()
# ax.loglog(k1d,spec_wthlv_mn[-1,izpl,:],label=r"$w'\theta_{lv}'$")
# ax.loglog(k1d,spec_wthlv_l_mn[-1,izpl,:],label=r"$\widetilde{w'}\widetilde{\theta_{lv}'}$")
# ax.loglog(k1d,spec_wthlv_c_mn[-1,izpl,:],label=r"$\widetilde{w'}\theta_{lv}'''+w'''\widetilde{\theta_{lv}'}$")
# ax.loglog(k1d,spec_wthlv_r_mn[-1,izpl,:],label=r"$w'''\theta_{lv}'''$")
# # ax.loglog(k1d,sumtest[-1,izpl,:],label=r"sum")
# # ax.set_ylim((1e-4,1e2))
# ax.set_ylabel(r"$k\widehat{w\theta}_{lv}'$")
# ax.set_xlabel(r"Wavenumber [1/m]")
# ax.legend(loc='best',bbox_to_anchor=(1,1),ncol=len(plttime)//13+1)

# ax2 = ax.twiny()
# fig.subplots_adjust(bottom=0.22)
# ax2.xaxis.set_ticks_position('bottom')
# ax2.xaxis.set_label_position('bottom')
# ax2.spines['bottom'].set_position(('axes',-0.22))
# ax2.set_xlim((2*np.pi/ax.get_xlim()[0],2*np.pi/ax.get_xlim()[1]))
# ax2.set_xscale('log')
# ax2.set_xlabel('Wavelength [m]')
# plt.show()
