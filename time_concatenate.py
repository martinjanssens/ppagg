#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 16:54:04 2021

@author: janssens
"""

import numpy as np

lps = ['/scratch-shared/janssens/bomex100_e12',
       '/scratch-shared/janssens/bomex200_from100_12hr']
savePath = '/scratch-shared/janssens/bomex200_from100_12hr/ppagg_merged'

fnames = ['time',
          'plttime',
	  'zf',
          'qtpf_moist_time',
          'qtpf_dry_time',
          'qtpf_prod_moist_time',
          'qtpf_prod_dry_time',
          'qtpf_prod_moist_wex_time',
          'qtpf_prod_dry_wex_time',
          'qtpf_vdiv_moist_time',
          'qtpf_vdiv_dry_time',
          'qtpf_hdiv_moist_time',
          'qtpf_hdiv_dry_time',
          'qtpf_subs_moist_time',
          'qtpf_subs_dry_time',
          'qtpf_diff_moist_time',
          'qtpf_diff_dry_time',
          'thlvpf_moist_time',
          'thlvpf_dry_time',
          'thlvpf_prod_moist_time',
          'thlvpf_prod_dry_time',
          'thlvpf_vdiv_moist_time',
          'thlvpf_vdiv_dry_time',
          'thlvpf_hdiv_moist_time',
          'thlvpf_hdiv_dry_time',
          'thlvpf_subs_moist_time',
          'thlvpf_subs_dry_time',
          'thlvpf_diff_moist_time',
          'thlvpf_diff_dry_time',
          'thlvpp_moist_time',
          'thlvpp_dry_time',
          'thlvpp_prod_moist_time',
          'thlvpp_prod_dry_time',
          'thlvpp_vdiv_moist_time',
          'thlvpp_vdiv_dry_time',
          'thlvpp_hdiv_moist_time',
          'thlvpp_hdiv_dry_time',
          'thlvpp_subs_moist_time',
          'thlvpp_subs_dry_time',
          'thlvpp_diff_moist_time',
          'thlvpp_diff_dry_time',    
          'thl_av_time',
          'thlv_av_time',
          'qt_av_time',    
          'thlpf_moist_time',
          'thlpf_dry_time',
          'wff_moist_time',
          'wff_dry_time',
          'qlpf_moist_time',
          'qlpf_dry_time',
          'thlpp_moist_time',
          'thlpp_dry_time',
          'wfp_moist_time',
          'wfp_dry_time',
          'qlpp_moist_time',
          'qlpp_dry_time',
          'wthlpf_moist_time',
          'wthlpf_dry_time',
          'wqtpf_moist_time',
          'wqtpf_dry_time',
          'wqlpf_moist_time',
          'wqlpf_dry_time',
          'wthlvp_av_time',
          'wthlvpf_moist_time',
          'wthlvpf_dry_time',
          'wthlvpf_l_moist_time',
          'wthlvpf_l_dry_time',
          'wthlvpf_c_moist_time',
          'wthlvpf_c_dry_time',
          'wthlvpf_r_moist_time',
          'wthlvpf_r_dry_time',
          'wthlvpp_moist_time',
          'wthlvpp_dry_time',
          ]

def process(fname, lps, savePath):
    arrs_in = []
    for i in range(len(lps)):
        arrs_in.append(np.load(lps[i]+'/'+fname+'.npy'))
        if fname == 'plttime' and i>0:
            arrs_in[i] += arrs_in[i-1][-1]+1
    arr_out = np.concatenate(arrs_in,axis=0) # along time dimension
    if fname == 'zf':
        arr_out = arrs_in[0] # Assumes vertical grid stays the same over restart
    np.save(savePath+'/'+fname+'.npy', arr_out)

for f in fnames:
    process(f, lps, savePath)
