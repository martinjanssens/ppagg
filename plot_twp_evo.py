#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 09:47:12 2022

@author: janssens
"""


import numpy as np
import netCDF4 as nc

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




