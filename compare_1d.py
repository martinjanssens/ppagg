#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 10:41:51 2021

@author: janssens
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from dataloader import DataLoaderDALES, DataLoaderMicroHH

rp = '/home/hp200321/data'

lps = [rp+'/botany-6-768/runs/Run_40',
       rp+'/botany-6-1536-50/runs/Run_40'
      ]
labs = ['dx = 100m',
        'dx = 200m'
        ]
sp = '/home/hp200321/data/pp/figs'
mods = ['dales',
        'dales']

#rp = '/scratch-shared/janssens'

#lps = [rp+'/bomex200_from100_12hr',
#       rp+'/bomex100a5_from100_12hr',
#       rp+'/bomex200_fiso_from100_12hr',
#       rp+'/bomex100_e12',
#       rp+'/bomex50',
#       rp+'/tmp.bomex/bomex_200m',
#       rp+'/tmp.bomex/bomex_100m',
       # '/scratch-shared/janssens/tmp.bomex/bomex_50m/ppagg',
#       ]
#sp = '/scratch-shared/janssens/bomex_comparisons'

#labs = [
#        r'D1: $\Delta x = 200m$',
#        r'D2: $\Delta x = 200m$, a5',
#        r'D3: $\Delta x = 200m$, fiso',
#        r'D4: $\Delta x = 100m$',
#        r'D5: $\Delta x = 50m$',
#        r'M1: $\Delta x = 200m$',
#        r'M2: $\Delta x = 100m$',
        #r'M3 - $\Delta x = 50m$',
#        ]
#mods = [
#        'dales',
#        'dales',
#        'dales',
#        'dales',
#        'dales',
#        'microhh',
#        'microhh',
        # 'microhh',
#        ]

ls = ['-',
      '-.',
      (0, (3, 2, 1, 2, 1, 2)),
      '--',
      ':',
      '-',
      '--',
      # ':'
      ]

def mav(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w
#%% Plot vertical profiles at a single time

var  = 'qt'
xlab = r'$q_t$]'
tplt = 2
tav = 2
zmin = 100
zmax = 3000

fig=plt.figure(); ax=plt.gca()
for i in range(len(lps)):
    lp = lps[i]
    
    if mods[i] == 'dales':
        dl = DataLoaderDALES(lp)
        col =  plt.cm.RdYlBu(0.99)
    elif mods[i] == 'microhh':
        dl = DataLoaderMicroHH(lp)
        col  =  plt.cm.RdYlBu(0.7)
    
    time = dl.time1d/3600
    zt   = dl.zf
    
    itmin  = np.argmin(abs(tplt-time))
    itmax =  np.argmin(abs(tplt+tav-time))
        
    izs = np.argmin(abs(zmin-zt))
    ize = np.argmin(abs(zmax-zt))
    
    load_func = getattr(dl, 'load_'+var+'av')
    pltvar = np.mean(load_func(izs,ize)[itmin:itmax,:],axis=0)

    ax.plot(pltvar,zt[izs:ize],label=labs[i],color=col,linestyle=ls[i])
            
ax.set_ylabel('z [m]')
ax.set_xlabel(xlab if len(xlab)>0 else var.split('/')[-1])
# ax.set_title('t = %.2f'%(tplt) +' hr')
ax.legend(loc='best',bbox_to_anchor=(1,1))
ax.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
plt.savefig(sp+'/prof_'+var+'_t'+str(tplt)+'.pdf',bbox_inches='tight')
plt.show()

