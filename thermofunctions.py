#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copied from https://github.com/fjansson/cloudbotany/blob/cloud_botany/thermo.py

import numpy as np
import netCDF4 as nc
import scipy

# Python version of the thermodynamic calculations in DALES

# # Inputs
# itmin = 23
# itmax = 24
# di    = 1
# izmin = 0
# izmax = 80
# lp = '/scratch-shared/janssens/bomex200_e12'

# ds = nc.Dataset(lp+'/fielddump.001.nc')
# ds1= nc.Dataset(lp+'/profiles.001.nc')
# ds0= nc.Dataset(lp+'/tmser.001.nc')

# time  = np.ma.getdata(ds.variables['time'][:]) / 3600
# zf    = np.ma.getdata(ds.variables['zt'][:]) # Cell centres (f in mhh)
# zh    = np.ma.getdata(ds.variables['zm'][:]) # Cell edges (h in mhh)
# xf    = np.ma.getdata(ds.variables['xt'][:]) # Cell centres (f in mhh)
# yf    = np.ma.getdata(ds.variables['yt'][:]) # Cell centres (f in mhh)

# time1d = np.ma.getdata(ds1.variables['time'][:])
# rhobf = np.ma.getdata(ds1.variables['rhobf'][:])
# rhobh = np.ma.getdata(ds1.variables['rhobh'][:])

# plttime = np.arange(itmin, itmax, di)
# zflim = zf[izmin:izmax]

# constants as in DALES
tup     = 268.            # Temperature range over which mixed phase occurs (high)
tdn     = 253.            # Temperature range over which mixed phase occurs (low)
rd      = 287.04          # gas constant for dry air.
rv      = 461.5           # gas constant for water vapor.
cp      = 1004.           # specific heat at constant pressure (dry air).
rlv     = 2.53e6          # latent heat for vaporisation
pref0   = 1.e5            # standard pressure used in exner function.
grav    = 9.81            # gravity acceleration.

# exner function of pressure
def exnf(pres):
    return (pres/pref0)**(rd/cp)


# Saturation water pressure over liquid and ice
# as in DALES and 
#     D. M. Murphy and T. Koop 2005, "Review of the vapour
#     pressures of ice and supercooled water for atmospheric
#     applications."  Q. J. R. Meteorol. Soc. 131:1539.

def esatl(T):
     return np.exp(54.842763-6763.22/T-4.21*np.log(T)+0.000367*T+
         np.tanh(0.0415*(T-218.8))*(53.878-1331.22/T-9.44523*np.log(T)+ 0.014025*T))

def esati(T):
    return np.exp(9.550426-5723.265/T+3.53068*np.log(T)-0.00728332*T)

# ice - liquid ratio as function of temperature
# 1 for all liquid, 0 for all ice
def ilratio(T):
    return np.maximum(0.,np.minimum(1.,(T-tdn)/(tup-tdn)))

# saturation humidity as in DALES
# linear interpolation between qsatur over liquid and qsatur over ice
#
# note IFS makes a different choice - linear interpolation vapor pressure
def qsatur(T, pres):
    esl1 = esatl(T)
    esi1 = esati(T)

    # this breaks if vapor pressure > real pressure
    # i.e. above boiling point of water at current pressure
    # original:
#    qsatur = (ilratio(T)     *(rd/rv)*esl1 / (pres - (1.-rd/rv)*esl1) +
#              (1.-ilratio(T))*(rd/rv)*esi1 / (pres - (1.-rd/rv)*esi1) )
    # with clamp
    qsat  = (ilratio(T)     *(rd/rv)*esl1 / (pres - np.minimum( (1.-rd/rv)*esl1, pres*0.8) ) +
              (1.-ilratio(T))*(rd/rv)*esi1 / (pres - np.minimum( (1.-rd/rv)*esi1, pres*0.8) ) )

    
    return np.clip(qsat, 0, 0.9)

# get (T, ql) from (thl, qt) at pressure pres
# accepts scalars or numpy arrays
def T_and_ql(thl, qt, pres):
    
    T = exnf(pres) * thl # first guess
    qsat = qsatur(T, pres)

#    if (qt < qsat):  # early return if below saturation
#        return T, 0  # doesn't work with vectors (could use .all for speed)
#                     # not needed for correctness
#    else:            # above saturation

# find T, ql such that ...
#   ql = max(qt - qsatur, 0.)
#   thl = T/exnf(pres) - (rlv/(cp*exnf(pres))) * ql

    def thl_err(t):
        ql = np.maximum(qt - qsatur(t, pres), 0.)
        thl1 = t/exnf(pres) - (rlv/(cp*exnf(pres))) * ql
        return thl - thl1

    def thl_err_scalar(t, qt, pres, thl):
        ql = np.maximum(qt - qsatur(t, pres), 0.)
        thl1 = t/exnf(pres) - (rlv/(cp*exnf(pres))) * ql
        return thl - thl1

    # T = scipy.optimize.broyden1(thl_err, T, f_tol=1e-5)
    #print(thl_err(200), thl_err(330))
    if type(thl)==np.ndarray:
        T = np.zeros_like(thl)
        for i in range(len(thl)):
            try:
                T[i] = scipy.optimize.brentq(thl_err_scalar, 200, 330, xtol=1e-3, args = (qt[i], pres[i], thl[i]))
            except:
                print('T, ql calculation failed, assigning nan')
                T[i] = np.nan
    else:
        print('scalar version')
        T = scipy.optimize.brentq(thl_err_scalar, 200, 330, xtol=1e-3, args = (qt, pres))
        
    ql = np.maximum(qt - qsatur(T, pres), 0.)
    return T, ql


# base pressure profile, for the option  ibas_prf=3
# "use standard atmospheric lapse rate with surface temperature offset"
# for now this version is valid only below 11 km
def pressure(zf, ps=101300, thls=300):
    # zmat=(/11000.,20000.,32000.,47000./)           # heights of lapse rate table
    lapserate=[-6.5/1000., 0., 1./1000, 2.8/1000 ]   # lapse rate table
        
    tsurf=thls*(ps/pref0)**(rd/cp) # surface temperature
    zsurf = 0
    #pmat = np.exp((log(ps)*lapserate(1)*rd+np.log(tsurf+zsurf*lapserate(1))*grav-
    #               np.log(tsurf+zmat(1)*lapserate(1))*grav)/(lapserate(1)*rd))

    pb = np.exp((np.log(ps)*lapserate[0]*rd + np.log(tsurf+zsurf*lapserate[0])*grav-
                 np.log(tsurf+zf*lapserate[0])*grav)/(lapserate[0]*rd))
    
    tb = tsurf+lapserate[0]*(zf - zsurf)

    rhobf = pb / (rd*tb) # dry estimate

    return pb