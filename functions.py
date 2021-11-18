#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:10:43 2021

@author: janssens
"""
import numpy as np
from scipy import fftpack, ndimage
import matplotlib.pyplot as plt

# Constants
cp = 1004.
Rd = 287.04
Rv = 461.5
Lv = 2.53e6
grav = 9.81


def getRad(data):
    h  = data.shape[-2];  hc = h//2
    w  = data.shape[-1];  wc = w//2
   
    # create an array of integer radial distances from the center
    Y, X = np.ogrid[0:h, 0:w]
    r    = np.hypot(X - wc, Y - hc).astype(np.int)
    
    return r

def lowPass(field, mask): 
    # FFT
    F = fftpack.fft2(field)       # 2D FFT (no prefactor)
    F = fftpack.fftshift(F,axes=(-1,-2))       # Shift so k0 is centred
    if len(F.shape) == 2:
        F = F*mask
    elif len(F.shape) == 3:
        F = F*mask[np.newaxis,:,:]
    F = fftpack.fftshift(F,axes=(-1,-2))
    
    return np.real(fftpack.ifft2(F))

def scaleDecomposeFlux(xt,xp,yt,yp,mask):
    # *t -> Filtered quantity
    # *p -> Remainder
    
    Le = lowPass(xt*yt, mask)
    Cr = lowPass(xt*yp + xp*yt, mask)
    Re = lowPass(xp*yp, mask)
    
    return Le, Cr, Re

def ddzwx_2nd(wh,x,dzh,rhobf=None):
    # Assumes constant dz
    if type(rhobf) == type(None):
        rhobf = np.ones(x.shape[0])
    return ((wh[2:,:,:]*(rhobf[2:,np.newaxis,np.newaxis]*x[2:,:,:] + 
                         rhobf[1:-1,np.newaxis,np.newaxis]*x[1:-1,:,:]) - 
             wh[1:-1,:,:]*(rhobf[:-2,np.newaxis,np.newaxis]*x[:-2,:,:]+
                           rhobf[1:-1,np.newaxis,np.newaxis]*x[1:-1,:,:]))/
             (2*dzh*rhobf[1:-1,np.newaxis,np.newaxis]))

def kddza_2nd(k,x,zf,rhobf=None):
    # Assumes constant dz
    if type(rhobf) == type(None):
        rhobf = np.ones(k.shape)
    if k.shape != x.shape:
        x = x[1:-1,:,:]
        rhobf = rhobf[1:-1]
    dzh = np.diff(zf)[0]
    return (k[1:-1,:,:]*((rhobf[2:,np.newaxis,np.newaxis]*x[2:,:,:]+
                           rhobf[1:-1,np.newaxis,np.newaxis]*x[1:-1,:,:]) - 
                         (rhobf[:-2,np.newaxis,np.newaxis]*x[:-2,:,:] + 
                          rhobf[1:-1,np.newaxis,np.newaxis]*x[1:-1,:,:]))/
             (2*dzh*rhobf[1:-1,np.newaxis,np.newaxis]))

def ddxhuha_2nd(u,v,a,dx,dy):
    # Assumes constant dx, dy, shape (z,y,x)
    return ((np.roll(u,-1,axis=2)*(np.roll(a,-1,axis=2) + a) -
             u                  *(np.roll(a,1,axis=2) + a))/(2*dx) +
            (np.roll(v,-1,axis=1)*(np.roll(a,-1,axis=1) + a) -
             v                  *(np.roll(a,1,axis=1) + a))/(2*dy))

def wsubdxdz(wfls,x,dzh):
    # Assumes constant dz, wfls < 0, returns values starting at full level 0
    # wsub defined at cell centres -> Interpolate to edges 
    whls = (wfls[1:] +  wfls[:-1])*0.5
    # upwinding
    return whls[:,np.newaxis,np.newaxis]*(x[1:,:,:] - x[:-1,:,:])/dzh

def compute_dthvdz(thl,qt,ql,exnf,dzh):
    # Using 2nd order scheme, ignores first model level
    
    # Bodges to have an upper boundary (which will be screwed up)
    exnf = np.append(exnf,exnf[-1]) 
    
    thv = (thl + Lv*ql/(cp*exnf[:,np.newaxis,np.newaxis]))* (1+(Rv/Rd-1)*qt-Rv/Rd*ql)
    dthvdz = (thv[2:,:,:] - thv[:-2])/(2*dzh)
    return dthvdz

def compute_ek(e12,dthvdz,thv_av,delta,cf=2.5,cn=0.76,Pri=3,ekmin=1e-6):
    # Assumes TKE-eq scheme, lmason=lanisotrop=ldelta=False
    alpha_kolm=1.5
    cm = cf / (2. * np.pi) * (1.5*alpha_kolm)**(-1.5)

    # Bodge to get shape right    
    e12 = e12[1:-1,:,:] # Remove first, last model level where you don't have dthvdz
    delta = delta[1:]
    thv_av = thv_av[1:]
        
    # Stable regime correction
    zlt = np.ones(e12.shape)*delta[:,np.newaxis,np.newaxis]
    mask = ( grav*np.abs(dthvdz) * zlt**2 > (cn*e12**2 * thv_av[:,np.newaxis,np.newaxis] ))
    zlt[mask] = cn*e12[mask]/(np.sqrt(grav/thv_av[:,np.newaxis,np.newaxis]*np.abs(dthvdz)))[mask]
    
    ekm = cm*zlt*e12
    ekh = (1 + (Pri-1)*zlt/delta[:,np.newaxis,np.newaxis])*ekm
    
    ekm[ekm<ekmin] = ekmin
    ekh[ekh<ekmin] = ekmin
    
    return ekm,ekh

def diffzeka(ekh,a,dzh,rhobf,rhobh):
    # Set shapes to ekh shape
    rhobf = rhobf[1:-1]
    rhobh = rhobh[1:-1]
    if ekh.shape != a.shape:
        a = a[1:-1,:,:]
    
    return 0.5*(rhobh[2:,np.newaxis,np.newaxis]/rhobf[1:-1,np.newaxis,np.newaxis]
                *(ekh[1:-1,:,:] + ekh[2:,:,:])*(a[2:,:,:] - a[1:-1,:,:]) / dzh**2
              -
               rhobh[1:-1,np.newaxis,np.newaxis]/rhobf[1:-1,np.newaxis,np.newaxis]
               *(ekh[1:-1,:,:] + ekh[:-2,:,:])*(a[1:-1,:,:] - a[:-2,:,:]) / dzh**2)

def diffeka(ekh,a,dx,dy,zf,rhobf=None,rhobh=None):
    # Assumes constant vertical spacing
    if type(rhobf) == type(None):
        rhobf = np.ones(zf.shape)
    if type(rhobh) == type(None):
        rhobh = np.ones(zf.shape)
    
    # bodge
    dzh = np.diff(zf)[0]
    
    # set to ekh shape
    a = a[1:-1,:,:]
    
    diffxh = 0.5*(((np.roll(ekh,-1,axis=2) + ekh)*(np.roll(a,-1,axis=2) - a)
                -(ekh + np.roll(ekh,1,axis=2))*(a - np.roll(a,1,axis=2)))/(dx**2)
                +
                ((np.roll(ekh,-1,axis=1) + ekh)*(np.roll(a,-1,axis=1) - a)
                -(ekh + np.roll(ekh,1,axis=1))*(a - np.roll(a,1,axis=1)))/(dy**2))
    diffz = diffzeka(ekh, a, dzh, rhobf, rhobh)
    
    return diffxh[1:-1,:,:] + diffz

def mean_mask(field,mask):
    masked = np.ma.masked_equal(field*mask,0)
    if len(field.shape) == 2:
        return masked.mean()
    elif len(field.shape) == 3:
        return masked.mean(axis=(1,2))
    else:
        print('Input field has wrong shape')
        return

def plot_2d(pltinp,xf,vmin=-1,vmax=1,fluct=True,cmap='RdYlBu'):
    if fluct:
        pltvar = pltinp - np.mean(pltinp)
    else:
        pltvar = pltinp
    plt.figure(figsize=(5.5,5))
    plt.imshow(pltvar,extent=np.array([xf.min(), xf.max(), xf.min(), xf.max()])/1000,
               vmin=vmin,vmax=vmax,cmap=cmap)
    plt.colorbar()
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.show()
    
def get_rad(data):
    # From https://medium.com/tangibit-studios/2d-spectrum-characterization-e288f255cc59
    h = data.shape[0]
    hc = h // 2
    w = data.shape[1]
    wc = w // 2

    # create an array of integer radial distances from the center
    Y, X = np.ogrid[0:h, 0:w]
    r = np.hypot(X - wc, Y - hc)

    return r


def get_psd_1d_radial(psd_2d, dx):

    # TODO: Assumes even number of points in psd_2d and square domain

    N = np.min(psd_2d.shape)
    L = int(N * dx)

    # Index radii corresponding to horizontal wavenumber 2*pi/L*r
    r = get_rad(psd_2d)
    r_int = np.round(r).astype(int)
    rp = np.arange(1, N // 2 + 1)

    # SUM all psd_2d pixels with label 'kh' for 0<=kh<=N/2 * 2*pi*L
    # Will miss power contributions in 'corners' kh>N/2 * 2*pi*L
    # Ignores centre (mean, wavenumber 0)
    # This is still a variance quantity.
    psd_1d = ndimage.sum(psd_2d, r_int, index=rp)

    # Compute prefactor that converts to spectral density and corrects for
    # annulus discreteness
    Ns = ndimage.sum(np.ones(psd_2d.shape), r_int, index=rp)

    kp = 2 * np.pi / L * ndimage.sum(r, r_int, index=rp) / Ns

    psd_1d *= L ** 2 * kp / (2 * np.pi * N ** 2 * Ns)

    return psd_1d

def compute_spectrum(cloud_scalar,dx,cloud_scalar_2=None,sqrt=False):
    # FFT
    F = fftpack.fft2(cloud_scalar)  # 2D FFT (no prefactor)
    F = fftpack.fftshift(F)  # Shift so k0 is centred
    if type(cloud_scalar_2) == type(None):
        # psd
        psd_2d = np.abs(np.conjugate(F) * F)
    else:
        # coherence magnitude spectrum
        F2 = fftpack.fft2(cloud_scalar_2)
        F2 = fftpack.fftshift(F2)
        psd_2d = np.abs(np.conjugate(F) * F2)
    if sqrt:
        psd_2d = np.sqrt(psd_2d)
    psd_2d /= np.prod(cloud_scalar.shape)  # Energy-preserving 2D PSD
    psd_1d_rad = get_psd_1d_radial(psd_2d, dx)  # Azimuthal integral-> 1D radial PSD

    # Wavenumbers
    N = np.min(cloud_scalar.shape)
    L = dx * N
    k1d = 2 * np.pi / L * np.arange(1, N // 2 + 1)
    
    return k1d, psd_1d_rad  
