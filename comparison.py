#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:35:15 2021

@author: janssens
"""

## Just lines to copy into a console with preloaded 100m res run for now
runcell(0, '/nfs/home4/janssens/scripts/pp3d/stats3d_eco_load.py')
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


plt.plot(np.mean(wff_moist_time[35:39,:],axis=0),zflim)
plt.plot(np.mean(wff_moist_time100[23:25,:],axis=0),zflim)
plt.plot(np.mean(wff_dry_time[35:39,:],axis=0),zflim)
plt.plot(np.mean(wff_dry_time100[23:25,:],axis=0),zflim)