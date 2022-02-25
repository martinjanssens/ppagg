#!/usr/bin/env python3

import numpy as np
from functions import zderivef

def load_ppagg(dl, lp):
    ld = {}

    ld['time'] = np.load(lp+'/time.npy')
    ld['zf']   = dl.zf_inp

    ld['time1d'] = dl.time1d
    ld['rhobf'] = dl.rhobf

    dzh = np.diff(ld['zf'])[0] # FIXME only valid if spacing is constant. Only used internally
    
    # Larger-scale subsidence
    ld['wfls'] = dl.wfls

    ld['plttime'] = np.load(lp+'/plttime.npy')
    ld['zflim'] = np.load(lp+'/zf.npy')

    ld['izmin'] = np.where(ld['zflim'][0] == ld['zf'])[0][0]
    ld['izmax'] = np.where(ld['zflim'][-1] == ld['zf'])[0][0]+1

    ld['qtpf_moist_time'] = np.load(lp+'/qtpf_moist_time.npy')
    ld['qtpf_dry_time'] = np.load(lp+'/qtpf_dry_time.npy')
    ld['qtpf_prod_moist_time'] = np.load(lp+'/qtpf_prod_moist_time.npy')
    ld['qtpf_prod_dry_time'] = np.load(lp+'/qtpf_prod_dry_time.npy')
    ld['qtpf_prod_wex_moist_time'] = np.load(lp+'/qtpf_prod_moist_wex_time.npy')
    ld['qtpf_prod_wex_dry_time'] = np.load(lp+'/qtpf_prod_dry_wex_time.npy')
    ld['qtpf_vdiv_moist_time'] = np.load(lp+'/qtpf_vdiv_moist_time.npy')
    ld['qtpf_vdiv_dry_time'] = np.load(lp+'/qtpf_vdiv_dry_time.npy')
    ld['qtpf_hdiv_moist_time'] = np.load(lp+'/qtpf_hdiv_moist_time.npy')
    ld['qtpf_hdiv_dry_time'] = np.load(lp+'/qtpf_hdiv_dry_time.npy')
    ld['qtpf_subs_moist_time'] = np.load(lp+'/qtpf_subs_moist_time.npy')
    ld['qtpf_subs_dry_time'] = np.load(lp+'/qtpf_subs_dry_time.npy')

    ld['thlvpf_moist_time'] = np.load(lp+'/thlvpf_moist_time.npy')
    ld['thlvpf_dry_time'] = np.load(lp+'/thlvpf_dry_time.npy')
    ld['thlvpf_prod_moist_time'] = np.load(lp+'/thlvpf_prod_moist_time.npy')
    ld['thlvpf_prod_dry_time'] = np.load(lp+'/thlvpf_prod_dry_time.npy')
    ld['thlvpf_vdiv_moist_time'] = np.load(lp+'/thlvpf_vdiv_moist_time.npy')
    ld['thlvpf_vdiv_dry_time'] = np.load(lp+'/thlvpf_vdiv_dry_time.npy')
    ld['thlvpf_hdiv_moist_time'] = np.load(lp+'/thlvpf_hdiv_moist_time.npy')
    ld['thlvpf_hdiv_dry_time'] = np.load(lp+'/thlvpf_hdiv_dry_time.npy')
    ld['thlvpf_subs_moist_time'] = np.load(lp+'/thlvpf_subs_moist_time.npy')
    ld['thlvpf_subs_dry_time'] = np.load(lp+'/thlvpf_subs_dry_time.npy')

    ld['thlvpp_moist_time'] = np.load(lp+'/thlvpp_moist_time.npy')
    ld['thlvpp_dry_time'] = np.load(lp+'/thlvpp_dry_time.npy')
    ld['thlvpp_prod_moist_time'] = np.load(lp+'/thlvpp_prod_moist_time.npy')
    ld['thlvpp_prod_dry_time'] = np.load(lp+'/thlvpp_prod_dry_time.npy')
    ld['thlvpp_vdiv_moist_time'] = np.load(lp+'/thlvpp_vdiv_moist_time.npy')
    ld['thlvpp_vdiv_dry_time'] = np.load(lp+'/thlvpp_vdiv_dry_time.npy')
    ld['thlvpp_hdiv_moist_time'] = np.load(lp+'/thlvpp_hdiv_moist_time.npy')
    ld['thlvpp_hdiv_dry_time'] = np.load(lp+'/thlvpp_hdiv_dry_time.npy')
    ld['thlvpp_subs_moist_time'] = np.load(lp+'/thlvpp_subs_moist_time.npy')
    ld['thlvpp_subs_dry_time'] = np.load(lp+'/thlvpp_subs_dry_time.npy')

    ld['wthlvpf_prod_moist_time'] = np.load(lp+'/wthlvpf_prod_moist_time.npy')
    ld['wthlvpf_prod_dry_time'] =  np.load(lp+'/wthlvpf_prod_dry_time.npy')
    ld['wthlvpf_vdiv_moist_time'] =  np.load(lp+'/wthlvpf_vdiv_moist_time.npy')
    ld['wthlvpf_vdiv_dry_time'] = np.load(lp+'/wthlvpf_vdiv_dry_time.npy')
    ld['wthlvpf_hdiv_moist_time'] = np.load(lp+'/wthlvpf_hdiv_moist_time.npy')
    ld['wthlvpf_hdiv_dry_time'] = np.load(lp+'/wthlvpf_hdiv_dry_time.npy')
    ld['wthlvpf_buoy_moist_time'] = np.load(lp+'/wthlvpf_buoy_moist_time.npy')
    ld['wthlvpf_buoy_dry_time'] = np.load(lp+'/wthlvpf_buoy_dry_time.npy')
    ld['wthlvpf_pres_moist_time'] = np.load(lp+'/wthlvpf_pres_moist_time.npy')
    ld['wthlvpf_pres_dry_time'] = np.load(lp+'/wthlvpf_pres_dry_time.npy')
    ld['wthlvpf_subs_moist_time'] = np.load(lp+'/wthlvpf_subs_moist_time.npy')
    ld['wthlvpf_subs_dry_time'] = np.load(lp+'/wthlvpf_subs_dry_time.npy')
    ld['wthlvpf_diff_moist_time'] = np.load(lp+'/wthlvpf_diff_moist_time.npy')
    ld['wthlvpf_diff_dry_time'] = np.load(lp+'/wthlvpf_diff_dry_time.npy')

    ld['thl_av_time'] = np.load(lp+'/thl_av_time.npy')
    ld['thlv_av_time'] = np.load(lp+'/thlv_av_time.npy')
    ld['qt_av_time'] = np.load(lp+'/qt_av_time.npy')

    ld['thlpf_moist_time'] = np.load(lp+'/thlpf_moist_time.npy')
    ld['thlpf_dry_time'] = np.load(lp+'/thlpf_dry_time.npy')
    ld['wff_moist_time'] = np.load(lp+'/wff_moist_time.npy')
    ld['wff_dry_time'] = np.load(lp+'/wff_dry_time.npy')
    ld['qlpf_moist_time'] = np.load(lp+'/qlpf_moist_time.npy')
    ld['qlpf_dry_time'] = np.load(lp+'/qlpf_dry_time.npy')

    ld['thlpp_moist_time'] = np.load(lp+'/thlpp_moist_time.npy')
    ld['thlpp_dry_time'] = np.load(lp+'/thlpp_dry_time.npy')
    ld['wfp_moist_time'] = np.load(lp+'/wfp_moist_time.npy')
    ld['wfp_dry_time'] = np.load(lp+'/wfp_dry_time.npy')
    ld['qlpp_moist_time'] = np.load(lp+'/qlpp_moist_time.npy')
    ld['qlpp_dry_time'] = np.load(lp+'/qlpp_dry_time.npy')

    ld['wthlpf_moist_time'] = np.load(lp+'/wthlpf_moist_time.npy')
    ld['wthlpf_dry_time'] = np.load(lp+'/wthlpf_dry_time.npy')

    ld['wqtpf_moist_time'] = np.load(lp+'/wqtpf_moist_time.npy')
    ld['wqtpf_dry_time'] = np.load(lp+'/wqtpf_dry_time.npy')

    ld['wqlpf_moist_time'] = np.load(lp+'/wqlpf_moist_time.npy')
    ld['wqlpf_dry_time'] = np.load(lp+'/wqlpf_dry_time.npy')
    ld['wthlvp_av_time'] = np.load(lp+'/wthlvp_av_time.npy')
    ld['wthlvpf_moist_time'] = np.load(lp+'/wthlvpf_moist_time.npy')
    ld['wthlvpf_dry_time'] = np.load(lp+'/wthlvpf_dry_time.npy')
    ld['wthlvpf_l_moist_time'] = np.load(lp+'/wthlvpf_l_moist_time.npy')
    ld['wthlvpf_l_dry_time'] = np.load(lp+'/wthlvpf_l_dry_time.npy')
    ld['wthlvpf_c_moist_time'] = np.load(lp+'/wthlvpf_c_moist_time.npy')
    ld['wthlvpf_c_dry_time'] = np.load(lp+'/wthlvpf_c_dry_time.npy')
    ld['wthlvpf_r_moist_time'] = np.load(lp+'/wthlvpf_r_moist_time.npy')
    ld['wthlvpf_r_dry_time'] = np.load(lp+'/wthlvpf_r_dry_time.npy')
    ld['wthlvpp_moist_time'] = np.load(lp+'/wthlvpp_moist_time.npy')
    ld['wthlvpp_dry_time'] = np.load(lp+'/wthlvpp_dry_time.npy')

    ld['wthlvpf_anom_moist_time'] = ld['wthlvpf_moist_time'] - ld['wthlvp_av_time']
    ld['wthlvpf_anom_dry_time'] = ld['wthlvpf_dry_time'] - ld['wthlvp_av_time']

    ld['Gamma_qt_av_time'] = zderivef(ld['qt_av_time'],dzh)
    ld['Gamma_thlv_av_time'] = zderivef(ld['thlv_av_time'],dzh)
    ld['Gamrat_av_time'] = ld['Gamma_qt_av_time']/ld['Gamma_thlv_av_time']
    ld['Gamrat_av_time'][np.abs(ld['Gamrat_av_time'])>0.03] = np.nan
    
    ## Reconstruct slab-mean budget terms
    ## FIXME Not working yet, would need support for time1d vs time and handling different time dimension sizes in restart and original
    # thl_av_1d = ds1['thl'][:,ld[i]['izmin']:ld[i]['izmax']]
    # qt_av_1d = ds1['qt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # thlv_av_1d = thl_av_1d*(1 + 0.608*qt_av_1d)

    # # Tendencies
    # ld[i]['ddt_thlv_av_time'] = tderive(thlv_av_1d, ld[i]['time1d']/3600)
    # ld[i]['ddt_qt_av_time'] = tderive(qt_av_1d, ld[i]['time1d']/3600)

    # # Flux divergence (approximately, i.e. ignoring rho)
    # ld[i]['wthl_av_time'] = ds1['wthlt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # ld[i]['wqt_av_time'] = ds1['wqtt'][:,ld[i]['izmin']:ld[i]['izmax']]
    # ld[i]['wthlv_av_time'] = ld[i]['wthl_av_time'] + 0.608*thl_av_1d*ld[i]['wqt_av_time']

    # ld[i]['ddz_wthlv_av_time'] = zderivef(ld[i]['wthlv_av_time'],dzh)
    # ld[i]['ddz_wqt_av_time'] = zderivef(ld[i]['wqt_av_time'],dzh)

    return ld
