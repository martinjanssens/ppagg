#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:47:12 2022

@author: janssens
"""


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ppagg_io import load_ppagg
from functions import getRad, lowPass, vint, mean_mask
from dataloader import DataLoaderDALES, DataLoaderMicroHH

lps = ['/scratch-shared/janssens/bomex200_e12/ppagg',
       '/scratch-shared/janssens/bomex200a5',
       '/scratch-shared/janssens/bomex200_fiso_from100_12hr/ppagg_merged',
       '/scratch-shared/janssens/bomex100_e12/ppagg',
       '/scratch-shared/janssens/bomex50',
       '/scratch-shared/janssens/tmp.bomex/bomex_200m/ppagg',
       '/scratch-shared/janssens/tmp.bomex/bomex_100m/ppagg',
       # '/scratch-shared/janssens/tmp.bomex/bomex_50m/ppagg',
       ]
sp = '/scratch-shared/janssens/bomex_comparisons'

labs = [
        r'D1: $\Delta x = 200m$',
        r'D2: $\Delta x = 200m$, a5',
        r'D3: $\Delta x = 200m$, fiso',
        r'D4: $\Delta x = 100m$',
        r'D5: $\Delta x = 50m$',
        r'M1: $\Delta x = 200m$',
        r'M2: $\Delta x = 100m$',
        #r'M3 - $\Delta x = 50m$',
        ]
mods = [
        'dales',
        'dales',
        'dales',
        'dales',
        'dales',
        'microhh',
        'microhh',
        # 'microhh',
        ]

src = ['ppagg',
       'cape',
       'ppagg',
       'ppagg',
       'cape',
       'ppagg',
       'ppagg',
       # 'ppagg',
       ]

ls = ['-',
      '-.',
      (0, (3, 2, 1, 2, 1, 2)),
      '--',
      ':',
      '-',
      '--',
      # ':'
      ]

tmin = 6.
tmax = [
        18.,
        36.,
        24.,
        36.,
        36.,
        12.,
        36.,
        # 36.
        ]

klp = 4
alpha=0.9
lw=2

col_moist_mhh = plt.cm.RdYlBu(0.7)
col_dry_mhh = plt.cm.RdYlBu(0.3)
col_moist_dal = plt.cm.RdYlBu(0.99)
col_dry_dal = plt.cm.RdYlBu(0.01)

f1 = plt.figure(figsize=(5,10/3)); axs1 = plt.gca()
for i in range(len(lps)):
    lp = lps[i]
 
    if src[i] == 'cape':
        # FIXME assumes it's DALES output
        ds = nc.Dataset(lp+'/cape2d.001.nc')
        time  = np.ma.getdata(ds.variables['time'][:]) / 3600

        itpltmin = np.where(time>=tmin)[0][0]
        itpltmax = np.where(time<tmax[i])[0][-1]+1
        plttime_var = np.arange(itpltmin,itpltmax,1)

        sz = ds.dimensions['xt'].size
        circ_mask = np.zeros((sz,sz))
        rad = getRad(circ_mask)
        circ_mask[rad<=klp] = 1

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
        
        [exp_moist,fac_moist], cov = curve_fit(lambda x, a, b: b * np.exp(x*a), 
                                               time[plttime_var]*3600, 
                                               twppf_moist,
                                               p0=[1e-5,0])
        tau = 1/exp_moist
        
    elif src[i] == 'ppagg':
        if mods[i] == 'dales':
            dl = DataLoaderDALES(lp+'/..')
        elif mods[i] == 'microhh':
            dl = DataLoaderMicroHH(lp+'/..')
        ld = load_ppagg(dl, lp)
        time = ld['time']

        itpltmin = np.where(time>=tmin)[0][0]
        itpltmax = np.where(time<tmax[i])[0][-1]+1
        plttime_var = np.arange(itpltmin,itpltmax,1)

        zflim = ld['zflim']
        rhobfi = ld['rhobf'][0,ld['izmin']:ld['izmax']]
        qtpfmi = ld['qtpf_moist_time']
        qtpfdi = ld['qtpf_dry_time']

        wthlvpf_moist_anom = ld['wthlvpf_anom_moist_time']
        wthlvpf_dry_anom = ld['wthlvpf_anom_dry_time']
        Gamrat = ld['Gamrat_av_time']

        twppf_moist = vint(qtpfmi,rhobfi,zflim,plttime_var)
        twppf_dry = vint(qtpfdi,rhobfi,zflim,plttime_var)
        
        # Time scale estimate
        Gamratz = (Gamrat[:,1:] - Gamrat[:,:-1])/(zflim[1] - zflim[0])
        wthlvpf_anomi_moist = -vint(wthlvpf_moist_anom[:,2:-1]*Gamratz,rhobfi[2:-1], zflim[2:-1], plttime=plttime_var)
        wthlvpf_anomi_dry = -vint(wthlvpf_dry_anom[:,2:-1]*Gamratz,rhobfi[2:-1], zflim[2:-1], plttime=plttime_var)

        coef = np.polyfit(np.concatenate((twppf_dry, twppf_moist)),
                          np.concatenate((wthlvpf_anomi_dry, wthlvpf_anomi_moist)), 1)
        tau = 1./coef[0]

        # Plot for debugging
        # qtpfi_mod = np.linspace(twppf_dry.min(),twppf_moist.max(),10)
        # wthlvpf_anomi_mod = qtpfi_mod  / tau
        # fs=14
        # fig = plt.figure(); ax = plt.gca()
        # ax.scatter(np.concatenate((twppf_dry, twppf_moist)),
        #            np.concatenate((wthlvpf_anomi_dry, wthlvpf_anomi_moist)),c='k',s=0.5)
        # ax.plot(qtpfi_mod,wthlvpf_anomi_mod,'k')
        # ax.set_ylabel(r"$-\left\langle F_{{\theta_{lv}}_m'} \frac{\partial}{\partial z}\left(\frac{\Gamma_{q_t}}{\Gamma_{\theta_{lv}}}\right) \right\rangle$ [kg/$m^2$/s]", fontsize=fs)
        # ax.set_xlabel(r"$\left\langle q_{t_m}'\right\rangle$ [kg/m$^2$]", fontsize=fs)
        # ax.annotate(r"$\tau_{q_{t_m}'} =$ %.2f hr"%(tau/3600), (0.55,0.06), xycoords='axes fraction', fontsize=14)
        # plt.show()
        
    if mods[i] == 'dales':
        col_moist = col_moist_dal
        col_dry = col_dry_dal
    elif mods[i] == 'microhh':
        col_moist = col_moist_mhh
        col_dry = col_dry_mhh
    axs1.plot(time[plttime_var],twppf_moist,c=col_moist,linestyle=ls[i],lw=lw,alpha=alpha)
    axs1.plot(time[plttime_var],twppf_dry,c=col_dry,linestyle=ls[i],lw=lw,alpha=alpha)
    labs[i]=labs[i]+r", $\tau_{q_{t_m}'} =$ %.1f hr"%(tau/3600)
axs1.set_xlabel('Time [hr]')
axs1.set_ylabel(r"$TWP_m'$ [kg/m$^2$]")
lines = axs1.get_lines()
dalind = [i for i in range(len(mods)) if mods[i] == 'dales']
mhhind = [i for i in range(len(mods)) if mods[i] == 'microhh']
leg1 = plt.legend([lines[2*i] for i in dalind],[labs[i] for i in dalind],title='DALES',bbox_to_anchor=(1,1.05),loc='upper left')
leg2 = plt.legend([lines[2*i] for i in mhhind],[labs[i] for i in mhhind],title='MicroHH',bbox_to_anchor=(1,0.3),loc='upper left')
axs1.add_artist(leg1)
axs1.add_artist(leg2)
# axs1.legend(loc='upper left',bbox_to_anchor=(1,1),ncol=2)
plt.savefig(sp+'/twp_evo_num.pdf',bbox_inches='tight')


