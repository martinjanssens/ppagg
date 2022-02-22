#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:47:12 2022

@author: janssens
"""


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
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
        r'$\Delta x = 200m$',
        r'$\Delta x = 200m$, a5',
        r'$\Delta x = 200m$, fiso',
        r'$\Delta x = 100m$',
        r'$\Delta x = 50m$',
        r'$\Delta x = 200m$',
        r'$\Delta x = 100m$',
        # r'$\Delta x = 50m$',
        ]
mods = ['dales',
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
tmax = [18.,
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

col_moist_mhh = plt.cm.RdYlBu(0.8)
col_dry_mhh = plt.cm.RdYlBu(0.2)
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

        twppf_moist = vint(qtpfmi,rhobfi,zflim,plttime_var)
        twppf_dry = vint(qtpfdi,rhobfi,zflim,plttime_var)

    if mods[i] == 'dales':
        col_moist = col_moist_dal
        col_dry = col_dry_dal
    elif mods[i] == 'microhh':
        col_moist = col_moist_mhh
        col_dry = col_dry_mhh
    axs1.plot(time[plttime_var],twppf_moist,c=col_moist,linestyle=ls[i],lw=lw,alpha=alpha,label=labs[i])
    axs1.plot(time[plttime_var],twppf_dry,c=col_dry,linestyle=ls[i],lw=lw,alpha=alpha)
axs1.set_xlabel('Time [hr]')
axs1.set_ylabel(r"$TWP_m'$ [kg/m$^2$]")
lines = axs1.get_lines()
dalind = [i for i in range(len(mods)) if mods[i] == 'dales']
mhhind = [i for i in range(len(mods)) if mods[i] == 'microhh']
leg1 = plt.legend([lines[2*i] for i in dalind],[labs[i] for i in dalind],title='DALES',bbox_to_anchor=(1,1),loc='upper left')
leg2 = plt.legend([lines[2*i] for i in mhhind],[labs[i] for i in mhhind],title='MicroHH',bbox_to_anchor=(1,0.3),loc='upper left')
axs1.add_artist(leg1)
axs1.add_artist(leg2)
# axs1.legend(loc='upper left',bbox_to_anchor=(1,1),ncol=2)
plt.savefig(sp+'/twp_evo_num.pdf',bbox_inches='tight')


