#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:47:12 2022

@author: janssens
"""


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from 
lps = ['/scratch-shared/janssens/bomex200_e12/ppagg',
       '/scratch-shared/janssens/bomex200a5',
       # '/scratch-shared/janssens/bomex200_fiso_from100_12hr',
       '/scratch-shared/janssens/bomex100_e12/ppagg',
       '/scratch-shared/janssens/bomex50',
       '/scratch-shared/janssens/tmp.bomex/bomex_200m/ppagg',
       '/scratch-shared/janssens/tmp.bomex/bomex_100m/ppagg',
       # '/scratch-shared/janssens/tmp.bomex/bomex_50m/ppagg',
       ]
labs = [
        r'$\Delta x = 200m$',
        r'$\Delta x = 100m$, a5',
        # r'$\Delta x = 200m$, fiso',
        r'$\Delta x = 100m$',
        r'$\Delta x = 50m$',
        r'$\Delta x = 200m$',
        r'$\Delta x = 100m$',
        r'$\Delta x = 50m$',
        ]
mods = ['dales',
        'dales',
        # 'dales',
        'dales',
        'dales',
        'microhh',
        'microhh',
        # 'microhh',
        ]

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
    lp = lps[i]
    if 'ppagg' not in lp.split('/')[-1]:
        if mods[i] == 'dales':
            dl = DataLoaderDALES(lp)
        elif mods[i] == 'microhh':
            dl = DataLoaderMicroHH(lp) 
    
    else:
        time = ld[i]['time']


    itpltmin = np.where(time>=tmin)[0][0]
    itpltmax = np.where(time<tmax[i])[0][-1]+1
    plttime_var = np.arange(itpltmin,itpltmax,1)

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


