#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 20:24:33 2020

@author: martinjanssens
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

class DataLoaderDALES:

    def __init__(self, load_path):
        print('Initialising dataloader...')
        self.lp = load_path
        self.ds = nc.Dataset(self.lp+'/fielddump.001.nc')
        self.ds1= nc.Dataset(self.lp+'/profiles.001.nc')
        self.ilp = np.loadtxt(self.lp+'/lscale.inp.001')

        self.time  = np.ma.getdata(self.ds.variables['time'][:]) / 3600
        self.zf    = np.ma.getdata(self.ds.variables['zt'][:])
        self.zh    = np.ma.getdata(self.ds.variables['zm'][:])
        self.xf    = np.ma.getdata(self.ds.variables['xt'][:]) 
        self.xh    = np.ma.getdata(self.ds.variables['xm'][:])
        self.yf    = np.ma.getdata(self.ds.variables['yt'][:])
        self.yh    = np.ma.getdata(self.ds.variables['ym'][:])

        self.time1d = np.ma.getdata(self.ds1.variables['time'][:])
        self.rhobf = np.ma.getdata(self.ds1.variables['rhobf'][:])
        self.rhobh = np.ma.getdata(self.ds1.variables['rhobh'][:])

        self.zf_inp = self.ilp[:,0]
        self.wfls = self.ilp[:,3]
        self.dqdt_ls = self.ilp[:,6]
        self.dthldt_ls = self.ilp[:,7]
        print('Set paths to all datasets and extracted dimensions')

    def load_qt(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['qt'][it,izmin:izmax,:,:])

    def load_wh(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['w'][it,izmin:izmax+1,:,:])

    def load_thl(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['thl'][it,izmin:izmax,:,:])

    def load_ql(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['ql'][it,izmin:izmax,:,:])

    def load_u(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['u'][it,izmin:izmax,:,:])

    def load_v(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['v'][it,izmin:izmax,:,:])

    def load_e12(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['e12'][it,izmin:izmax,:,:])

    def load_p(self, it, izmin, izmax):
        return np.ma.getdata(self.ds.variables['p'][it,izmin:izmax,:,:])

    def load_presh(self, it, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['presh'][it,izmin:izmax])

    def load_qtav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['qt'][:,izmin:izmax])

    def load_thlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['thl'][:,izmin:izmax])

    def load_qlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['ql'][:,izmin:izmax])

    def load_wqtav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wqtt'][:,izmin:izmax])

    def load_wthlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wthlt'][:,izmin:izmax])

    def load_wqtrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wqtr'][:,izmin:izmax])

    def load_wthlrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wthlr'][:,izmin:izmax])

    def load_wqlrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wqlr'][:,izmin:izmax])

    def load_wthvrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1.variables['wthvr'][:,izmin:izmax])

class DataLoaderMicroHH:

    def __init__(self, load_path, case='bomex'):
        print('Initialising dataloader...')
        self.lp = load_path
        self.ds1 = nc.Dataset(self.lp+'/'+case+'.default.0000000.nc')
        self.ilp = nc.Dataset(self.lp+'/'+case+'_input.nc')
        self.dsqt = nc.Dataset(self.lp+'/qt.nc')
        self.dsthl = nc.Dataset(self.lp+'/thl.nc')
        self.dsql = nc.Dataset(self.lp+'/ql.nc')
        self.dsw = nc.Dataset(self.lp+'/w.nc')
        self.dsu = nc.Dataset(self.lp+'/u.nc')
        self.dsv = nc.Dataset(self.lp+'/v.nc')
        self.dsp = nc.Dataset(self.lp+'/p.nc')

        self.time  = np.ma.getdata(self.dsqt.variables['time'][:]) / 3600
        self.zf    = np.ma.getdata(self.dsqt.variables['z'][:])
        self.zh    = np.ma.getdata(self.dsw.variables['zh'][:])
        self.xf    = np.ma.getdata(self.dsqt.variables['x'][:])
        self.xh    = np.ma.getdata(self.dsu.variables['xh'][:])
        self.yf    = np.ma.getdata(self.dsqt.variables['y'][:])
        self.yh    = np.ma.getdata(self.dsv.variables['yh'][:])

        self.time1d = np.ma.getdata(self.ds1.variables['time'][:])
        self.rhobf = np.ma.getdata(self.ds1['thermo']['rho'][:])
        self.rhobh = np.ma.getdata(self.ds1['thermo']['rhoh'][:])

        self.zf_inp = np.ma.getdata(self.ilp['z'][:])
        self.wfls = np.ma.getdata(self.ilp['init']['w_ls'][:])
        self.dqdt_ls = np.ma.getdata(self.ilp['init']['qt_ls'][:])
        self.dthldt_ls = np.ma.getdata(self.ilp['init']['thl_ls'][:])
        print('Set paths to all datasets and extracted dimensions')

    def load_qt(self, it, izmin, izmax):
        return np.ma.getdata(self.dsqt.variables['qt'][it,izmin:izmax,:,:])

    def load_wh(self, it, izmin, izmax):
        return np.ma.getdata(self.dsw.variables['w'][it,izmin:izmax+1,:,:])

    def load_thl(self, it, izmin, izmax):
        return np.ma.getdata(self.dsthl.variables['thl'][it,izmin:izmax,:,:])

    def load_ql(self, it, izmin, izmax):
        return np.ma.getdata(self.dsql.variables['ql'][it,izmin:izmax,:,:])

    def load_u(self, it, izmin, izmax):
        return np.ma.getdata(self.dsu.variables['u'][it,izmin:izmax,:,:])

    def load_v(self, it, izmin, izmax):
        return np.ma.getdata(self.dsv.variables['v'][it,izmin:izmax,:,:])

    def load_e12(self, it, izmin, izmax):
        raise NotImplementedError('MicroHH does not solve for e12')

    def load_p(self, it, izmin, izmax):
        return np.ma.getdata(self.dsp.variables['p'][it,izmin:izmax,:,:])

    def load_presh(self, it, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['phydroh'][it,izmin:izmax])

    def load_qtav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['qt'][:,izmin:izmax])

    def load_thlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['thl'][:,izmin:izmax])

    def load_qlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['ql'][:,izmin:izmax])

    def load_wqtav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['qt_flux'][:,izmin:izmax])

    def load_wthlav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['thl_flux'][:,izmin:izmax])

    def load_wqtrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['qt_w'][:,izmin:izmax])

    def load_wthlrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['thl_w'][:,izmin:izmax])

    def load_wqlrav(self, izmin, izmax):
        print('Warning: MicroHH does not output wql. Returning zeros...')
        return np.zeros((len(self.time1d),izmax-izmin))

    def load_wthvrav(self, izmin, izmax):
        return np.ma.getdata(self.ds1['thermo']['thv_w'][:,izmin:izmax])


