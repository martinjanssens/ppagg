#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:35:15 2021

@author: janssens
"""

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
plt.plot(np.mean(wff_moist_time[its:ite,:],axis=0),zflim,c='C0',label='iadv2, moist')
plt.plot(np.mean(wff_moist_time100[its:ite,:],axis=0),zflim,c='C0',linestyle='--',label='iadv5, moist')
plt.plot(np.mean(wff_dry_time[its:ite,:],axis=0),zflim,c='C1',label='iadv2, dry')
plt.plot(np.mean(wff_dry_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--',label='iadv5, dry')
plt.legend()
plt.show()

plt.plot(np.mean(wthlvpf_moist_time[its:ite,:]-wthlvp_av_time[its:ite,:],axis=0),zflim,c='C0',label='iadv2, moist')
plt.plot(np.mean(wthlvpf_moist_time100[its:ite,:]-wthlvp_av_time100[its:ite,:],axis=0),zflim,c='C0',linestyle='--',label='iadv5, moist')
plt.plot(np.mean(wthlvpf_dry_time[its:ite,:]-wthlvp_av_time[its:ite,:],axis=0),zflim,c='C1',label='iadv2, dry')
plt.plot(np.mean(wthlvpf_dry_time100[its:ite,:]-wthlvp_av_time100[its:ite,:],axis=0),zflim,c='C1',linestyle='--',label='iadv5, dry')
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

labbsl = r"a2"
labopt = r"a5 from bsl"

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