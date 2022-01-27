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

lp = '/scratch-shared/janssens/bomex100_e12'

time = np.load(lp+'/time_spec.npy')
plttime = np.load(lp+'/plttime_spec.npy')
zf = np.load(lp+'/zf_spec.npy')
k1d = np.load(lp+'/k1d.npy')

spec_qt = np.load(lp+'/spec_qt.npy')
spec_thl = np.load(lp+'/spec_thl.npy')
spec_thlv = np.load(lp+'/spec_thlv.npy')
spec_w = np.load(lp+'/spec_w.npy')
spec_ql = np.load(lp+'/spec_ql.npy')

spec_wqt = np.load(lp+'/spec_wqt.npy')
spec_wthl = np.load(lp+'/spec_wthl.npy')
spec_wthlv = np.load(lp+'/spec_wthlv.npy')
spec_wql = np.load(lp+'/spec_wql.npy')

def plot_spectrum(k1d, spec, lab, plttime):
    fig = plt.figure(); ax = plt.gca()
    for i in range(len(plttime)):
        col = plt.cm.cubehelix(i/len(plttime))
        ax.loglog(k1d,spec[i,izpl,:],c=col,label='t=%.2f'%time[plttime[i]])
    # ax.set_ylim((1e-4,1e2))
    ax.set_ylabel(lab)
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

#%% Average over time
itav = 4 # number of time steps to average over -> len(plttime) MUST BE MULTIPLE OF THIS

spec_qt_mn = block_reduce(spec_qt,(itav,1,1),func=np.mean)
spec_thl_mn = block_reduce(spec_thl,(itav,1,1),func=np.mean)
spec_thlv_mn = block_reduce(spec_thlv,(itav,1,1),func=np.mean)
spec_w_mn = block_reduce(spec_w,(itav,1,1),func=np.mean)
spec_ql_mn = block_reduce(spec_ql,(itav,1,1),func=np.mean)

spec_wqt_mn = block_reduce(spec_wqt,(itav,1,1),func=np.mean)
spec_wthl_mn = block_reduce(spec_wthl,(itav,1,1),func=np.mean)
spec_wthlv_mn = block_reduce(spec_wthlv,(itav,1,1),func=np.mean)
spec_wql_mn = block_reduce(spec_wql,(itav,1,1),func=np.mean)

plttime_mn = plttime[::itav]

#%% Plot
izpl = 9

# Variances
plot_spectrum(k1d, spec_qt_mn, r"$k\widehat{q}_t'^2$", plttime_mn)
plot_spectrum(k1d, spec_thl_mn, r"$k\widehat{\theta}_l'^2$", plttime_mn)
plot_spectrum(k1d, spec_thlv_mn, r"$k\widehat{\theta}_{lv}'^2$", plttime_mn)
plot_spectrum(k1d, spec_w_mn, r"$k\widehat{w}'^2$", plttime_mn)
plot_spectrum(k1d, spec_ql_mn, r"$k\widehat{q_l}'^2$", plttime_mn)

# Fluxes
plot_spectrum(k1d, spec_wqt_mn, r"$k\widehat{wq}_t'$", plttime_mn)
plot_spectrum(k1d, spec_wthl_mn, r"$k\widehat{w\theta}_l'$", plttime_mn)
plot_spectrum(k1d, spec_wthlv_mn, r"$k\widehat{w\theta}_{lv}'$", plttime_mn)
plot_spectrum(k1d, spec_wql_mn, r"$k\widehat{wq}_{l}'$", plttime_mn)

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
