#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:09:54 2019

@author: joshua
"""

import numpy as np
import math as m
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import astropy.stats as astat
from scipy.special import betainc 

def update_progress(progress):
    '''
    progress bar function
    '''    
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Finishing...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def bpt(o3, hb, n2, ha, unred=False):
    '''
    0 = Ambiguous
    1 = Star forming
    2 = Composite
    3 = LINER
    4 = Seyfert
    -1 = Invalid
    '''
    if unred == True:
        ebv = EBV(ha, hb)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv = ebv)
        ha = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
    
    o3hb = log10(divide(o3, hb))
    n2ha = log10(divide(n2, ha))
    
    clas = []
    for i in range(len(o3hb)):
        if np.isnan(o3hb[i])==True or np.isnan(n2ha[i])==True:
            '''Invalid'''
            clas.append(-1)
        else:
            if (o3hb[i]>=0.61/(n2ha[i]-0.47) + 1.19 or n2ha[i]>=0.47) and \
            (o3hb[i]>=1.05*n2ha[i]+0.45):
                '''Seyfert'''
                clas.append(4)
            elif (((o3hb[i]>=0.61/(n2ha[i]-0.47) + 1.19)or(n2ha[i]>=0.47)) and \
            (o3hb[i]<1.05*n2ha[i]+0.45)):
                '''LINER'''
                clas.append(3)
            elif ((o3hb[i]<0.61/(n2ha[i]-0.05)+1.30)and(n2ha[i]<0.05)):
                '''Star Forming'''
                clas.append(1)
            elif (((o3hb[i]>=0.61/(n2ha[i]-0.05)+1.30)or(n2ha[i]>0.05)) and \
            ((o3hb[i]<0.61/(n2ha[i]-0.47)+1.19)and(n2ha[i]<0.47))):
                '''Composite'''
                clas.append(2)
            else:
                '''Ambiguous'''
                clas.append(0)
    return np.asarray(clas)

def whan(n2, ha, ewha):
    n2ha = log10(divide(n2, ha))
    
    clas = []
    for i in range(len(n2ha)):
        if np.isnan(n2ha[i])==True or np.isnan(ewha[i])==True:
            clas.append(-1)
        else:
            if n2ha[i]>=-0.4 and ewha[i]>=6:
                clas.append(4)
            elif n2ha[i]>=-0.4 and ewha[i]<6 and ewha[i]>=3:
                clas.append(3)
            elif n2ha[i]<-0.4 and ewha[i]>3:
                clas.append(1)
            elif ewha[i]<3:
                clas.append(0)
    return np.asarray(clas)

def bpt_whan(o3, hb, n2, ha, ewha, unred=False):
    '''
    0 = Retired
    1 = Star forming
    2 = Composite
    3 = LINER
    4 = Seyfert
    -1 = Invalid
    -2 = Ambiguous
    '''
    if unred == True:
        ebv = EBV(ha, hb)
        o3 = ccm_unred([5008.]*len(o3), o3, ebv = ebv)
        n2 = ccm_unred([6585.]*len(n2), n2, ebv = ebv)
        ha = ccm_unred([6565.]*len(ha), ha, ebv = ebv)
        hb = ccm_unred([4863.]*len(hb), hb, ebv = ebv)
    
    o3hb = np.log10(np.divide(o3, hb))
    n2ha = np.log10(np.divide(n2, ha))

    bpt = np.ones(o3hb.shape)*-2
    bpt =  4 * (((o3hb>=0.61/(n2ha-0.47) + 1.19)|(n2ha>=0.47))&(o3hb>=1.05*n2ha+0.45)) \
         + 3 * (((o3hb>=0.61/(n2ha-0.47) + 1.19)|(n2ha>=0.47))&(o3hb<1.05*n2ha+0.45)) \
         + 2 * (((o3hb>=0.61/(n2ha-0.05)+1.30)|(n2ha>0.05))&((o3hb<0.61/(n2ha-0.47)+1.19)&(n2ha<0.47))) \
         + 1 * ((o3hb<0.61/(n2ha-0.05)+1.30)&(n2ha<0.05))
    
    bpt = -1 * (np.isnan(o3hb)|np.isnan(n2ha)|~np.isfinite(o3hb)|~np.isfinite(n2ha)) \
          + bpt * ~(np.isnan(o3hb)|np.isnan(n2ha)|~np.isfinite(o3hb)|~np.isfinite(n2ha))
    
    whan = 1 * (ewha > 3) + 0 * (ewha <= 3)

    clas = bpt * (whan == 1) + 0 * (whan == 0)
    
    return np.asarray(clas)

def ifu_rin(ifu):
    '''
    convert ifudesignsize to ifu_rin
    '''
    ifu_dia = 12.5 * (ifu==19) + 17.5 * (ifu==37) + 22.5 * (ifu==61) + 27.5 * (ifu==91) + 32.5 * (ifu==127)
    ifu_rad = (ifu_dia+3.5-2.31)/2.*(np.sqrt(3.)/2.0)
    return ifu_rad

def in_IFU(ifudesignsize, ifura, ifudec, objra, objdec):
    '''
    Check if object is within MaNGA IFU. Takes the ifudesignsize and converts it
    to its equivalent arcsecond size. It then builds a artificial WCS for the 
    IFU and converts the objects RA and Dec to a pixel coordinate. The function
    then uses shapely to deterine if the object's pixel coordinates are within 
    the IFU.
    
    ifudesignsize    fiber size of the IFU from MaNGA 19, 37, 61, 91, or 127
    ifura            RA of the IFU
    ifudec           Dec of the IFU
    objra            Array of object RAs to check 
    objdec           Array of object Decs to check
    
    returns          Array of True/False for if the object is in the IFU
    '''
    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon
    from astropy import wcs
    
    IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
    size = IFU_arcsec[1, np.where(ifudesignsize==IFU_arcsec[0])[0][0]]/60.0 #in arcminutes
    
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0.0, 0.0]
    w.wcs.cdelt = np.array([-size/60.0/2, size/60.0/2])
    w.wcs.crval = [ifura, ifudec]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    
    vert = np.empty((0,2))
    for j in range(6):
        vert = np.append(vert, np.array([[m.cos((j)*m.pi/3.0), m.sin((j)*m.pi/3.0)]]), axis=0)
    polygon = Polygon(vert)
    
    inside = []
    for j in range(len(objra)):
        point = Point(w.wcs_world2pix([[objra[j], objdec[j]]], 1)[0])
        inside.append(polygon.contains(point))
    inside = np.asarray(inside)
    return inside

def ccm_unred(wave, flux, a_v=None, ebv=None, r_v=3.1):
    from numpy.lib.polynomial import poly1d
    import numpy
        
    x = 10000. / numpy.array(wave)                # Convert to inverse microns 
    npts = x.size
    a = numpy.zeros(npts, dtype=numpy.float)
    b = numpy.zeros(npts, dtype=numpy.float)    
    good = numpy.where( (x >= 0.3) & (x < 1.1) )
    if len(good[0]) > 0:    
        a[good] = 0.574 * x[good] ** (1.61)
        b[good] = -0.527 * x[good] ** (1.61)
        
    good = numpy.where( (x >= 1.1) & (x < 3.3) )
    if len(good[0]) > 0:
        y = x[good] - 1.82
        #Use new constants from O'Donnell (1994)
        c1 = numpy.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])
        c2 = numpy.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])
        a[good] = poly1d(c1[::-1])(y)
        b[good] = poly1d(c2[::-1])(y)

    good = numpy.where( (x >=3.3) & (x < 8))
    if len(good[0]) > 0:        
        y = x[good]
        f_a = numpy.zeros([len(good[0])], dtype=numpy.float)    
        f_b = numpy.zeros([len(good[0])], dtype=numpy.float)
        good1 = numpy.where(np.ravel((y > 5.9)))[0]
        if len(good1[0]) > 0:    
            y1 = y[good1] - 5.9
            f_a[good1] = -0.04473 * y1 ** 2 - 0.009779 * y1 ** 3
            f_b[good1] = 0.2130 * y1 ** 2 + 0.1207 * y1 ** 3
        
        a[good] = 1.752 - 0.316 * y - divide(0.104 , ((y - 4.67) ** 2 + 0.341)) + f_a
        b[good] = -3.090 + 1.825 * y + divide(1.206 , ((y - 4.62) ** 2 + 0.263)) + f_b
    good = numpy.where( (x >= 8) & (x <= 11) )
    if len(good[0]) > 0:    
        y = x[good] - 8.
        c1 = numpy.array([-1.073, -0.628, 0.137, -0.070])
        c2 = numpy.array([13.670, 4.257, -0.420, 0.374])
        a[good] = poly1d(c1[::-1])(y)
        b[good] = poly1d(c2[::-1])(y)
    if a_v is None:
        a_v = r_v * ebv
    a_lambda = a_v * (a + b / r_v)
    funred = flux * 10. ** (0.4 * a_lambda)       #Derive unreddened flux
    return funred

# reddening curve
def a(x):
    from numpy.lib.polynomial import poly1d
    c = np.array([1., 0.104, -0.609, 0.701, 1.137, -1.718, -0.827, 1.647, -0.505])
    return  poly1d(c[::-1])(x)

def b(x):
    from numpy.lib.polynomial import poly1d
    c = np.array([0., 1.952, 2.908, -3.989, -7.985, 11.102, 5.491, -10.805, 3.347])
    return  poly1d(c[::-1])(x)

def alam_anu(x, Rv=3.1):
    x = x - 1.82
    return a(x) + b(x)/Rv

def klam(x, Rv=3.1):
    x = 10000./x
    return alam_anu(x, Rv) * Rv

def EBV(ha, hb):
    return (np.log(2.85/divide(ha, hb)) / (-0.339783))/3.1

def log10(x):
    return np.log10(x, out=np.zeros(x.shape), where=x>0)
def divide(x, y):
    return np.divide(x, y, out=np.zeros(y.shape), where=y!=0)

def relav(z1, z2):
    c = 299792.458
    B = ( (1+z1)**2 - (1+z2)**2 )/( (1+z1)**2 + (1+z2)**2 ) 
    return c*B

def binfit(bins, x, y, order):
    """
    Take scatter plot data, bin them and fit a polynomial function to it.
    """  
    bin_val = []
    bin_std = []
    for i in range(len(bins)-1):
        z = y[((x>=bins[i])&(x<bins[i+1]))]
        if len(z) <= 4:
            bin_val.append(np.nan)
            bin_std.append(np.nan)
        else:
            bin_val.append(np.median(reject_outliers(z, m=1)))
            bin_std.append(np.std(reject_outliers(z)))
    
    real = (~np.isnan(bin_val))
    bins = bins[0:-1][real]
    bin_val = np.asarray(bin_val)[real]
    bin_std = np.asarray(bin_std)[real]
#    function = np.poly1d(np.polyfit(bins, bin_val, order))
    function = [0]
    return bins, bin_val, bin_std, function

def reject_outliers(data, m = 3., index=False):
    d = np.abs(data - np.nanmedian(data))
    mdev = np.nanmedian(d)
    s = d/mdev if mdev else 0.
    if index is False:
        return data[s<m]
    elif index is True:
        return (s<m)
#______________________________________________________________________________
def SFR_K98(ha, hb, z, dha=None, dhb=None):
    """Kennicutt 1998 SFR Law"""
    c = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.27)
    d = c.luminosity_distance(z).value * 10**6 * 3.086*10**16# *100
    #kpc_arcsec = []
    #for i in range(len(z)):
    #    kpc_arcsec.append(cosmo(z[i]))
    #kpc_arcsec = np.asarray(kpc_arcsec)
    #ha = ccm_unred([6565.]*len(ha), ha, ebv = EBV(ha, hb))
    #ha_lum = ha * (1.0 * kpc_arcsec *3.086**21)**2 / 10**-17
    ha_lum = ha * 4*np.pi * d**2 * 10**-7
    #sfr = log10(ha_lum / (1.26 * 10**41))
    sfr = log10(ha_lum / (10**41.1))
    if dha is not None and dhb is not None:
        dsfr = sfr
    else:
        dsfr = np.zeros_like(ha)
    return sfr, dsfr

def SFR_M11(Lha):
    '''
    Input L_Ha in erg/s to get SFR in Mo/yr
    '''
    return 5.37*10**-42 * Lha

def line_lum(flux, z): 
    cosmo = FlatLambdaCDM(H0=70*u.km/u.s/u.Mpc, Om0=0.3)
    r = cosmo.luminosity_distance(z).value * 10**8 * (3.09*10**16)
    lum = flux * 10**-17 * (4 * np.pi * r**2)
    lum = lum.astype('float')
    return lum

def stack(x, z=False, zbin=False, error='std', trim=None):
    '''
    Function for stacking radial profiles. Can break up profiles to be binned
    along some 3rd dimension variable (like mass). The 3rd dim is z, the value
    for that dimension, and zbin, the grid to bin on.
    
    x should be a 2-D array where axis=1 is the radius dimension and axis=0 is 
    the radially binned value.
    
    The error can be set to std to return standard deviation or sem to return
    standard error of the mean.
    
    Trim is an optional argument which will, if set, only keep the data in a 
    bin if the bin contains at least the given percent of the total data. Input
    trim as a decimal.
    '''
    if z is False and zbin is False:
        val = np.empty(len(x[0]))*np.nan    
        err = np.empty(len(x[0]))
        for j in range(len(x[0])):
            xr = x[:, j]
            v = np.nanmedian(xr[~np.isnan(xr)&np.isfinite(xr)])
            val[j] = v
            if error=='std':
                err[j] = astat.biweight_scale(xr[~np.isnan(xr)&np.isfinite(xr)], c=6)
            elif error=='sem':
                err[j] = astat.biweight_scale(xr[~np.isnan(xr)&np.isfinite(xr)], c=6)/np.sqrt(len(xr))
        if error=='std':
            err = astat.biweight_scale(x[~np.isnan(x)&np.isfinite(x)], c=6)
        elif error=='sem':
            err = astat.biweight_scale(x[~np.isnan(x)&np.isfinite(x)], c=6)/np.sqrt(len(x))
    
    else:
        step = zbin[1]-zbin[0]
        zhig = zbin + step/2.
        zlow = zbin - step/2.
        
        if trim is None:
            lim = 0.0
        else:
            lim = trim
        
        val = np.empty((len(zbin), len(x[0])))*np.nan    
        err = np.empty((len(zbin), len(x[0])))
        for i in range(len(zbin)):        
            mix = (z >= zlow[i])&(z < zhig[i])
            xx = x[mix]
            
            for j in range(len(x[0])):
                if len(xx) > 1:
                    xr = xx[:, j]
                    if len(xr[~np.isnan(xr)&np.isfinite(xr)]) > len(xx) * lim:
                        v = np.nanmedian(xr[~np.isnan(xr)&np.isfinite(xr)])
                        val[i,j] = v
                        if error=='std':
                            err[i,j] = astat.biweight_scale(xr[~np.isnan(xr)&np.isfinite(xr)], c=6)
                        elif error=='sem':
                            err[i,j] = astat.biweight_scale(xr[~np.isnan(xr)&np.isfinite(xr)], c=6)/np.sqrt(len(xr))
            
            # Cutoff the profile once a radius bin becomes NaN. Pass over first 
            # bins where val is NaN.
            """
            start = 0
            end = 0
            for j in range(len(val[i])):
                if np.isnan(val[i,j]) and start == 1 and end == 0:
                    end += 1
                if ~np.isnan(val[i,j]) and start == 0:
                    start += 1
                
                if start == 1 and end == 1:
                    val[i,j] = np.nan
                    err[i,j] = np.nan
            """
            # remove any single points
            for j in range(len(val[i])-2):
                j+=1
                if np.isnan(val[i,j-1]) and np.isnan(val[i,j+1]):
                    val[i,i] = np.nan
                    err[i,j] = np.nan
            
    return val, err

def Rotate2D(pts,cnt,ang):
    '''pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian'''
    return np.dot(pts-cnt,np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]]))+cnt

def dist_ellipse(x, y, x0, y0, ratio, pos_ang, q):
    ang = pos_ang * np.pi/180.0

    cosang = np.cos(ang)
    sinang = np.sin(ang)
    
    ratio = np.sqrt(((ratio)**2 - q**2)/(1 - q**2))
    
    x = x - x0
    y = y - y0
    
    xcosang = x*cosang
    xsinang = x*sinang

    xtemp = xcosang + y*sinang
    ytemp = -1*xsinang + y*cosang  
    r = np.sqrt( (xtemp/ratio)**2 + ytemp**2)
    return r

def azimuth_profiles(x, y, xbins):
    bin_val = [] 
    bin_std = []
    n = []
    for j in range(len(xbins)-1):
        ix = (x >= xbins[j])&(x < xbins[j+1])
        if len( y[ix] ) > 0:
            #print y[ix]
            bin_val.append(np.nanmedian(y[ix]))
            bin_std.append(np.nanstd(y[ix]))
            n.append(len(y[ix]))
        else:
            #print [np.nan]
            bin_val.append(np.nan)
            bin_std.append(np.nan)
            n.append(0)
            
        #print bin_val
    return np.asarray(bin_val), np.asarray(bin_std), np.asarray(n)
    
def binormial_er(a, b, cl=0.683):
    '''
    error between 
    '''
    frac = np.empty((a.shape))
    lower = np.empty((a.shape))
    upper = np.empty((a.shape))
    for i in range(len(a)):
        k = np.float64(a[i])
        n = np.float64(b[i])
        kn = np.divide(k, n, out=np.zeros(n.shape), where=n!=0)
        z = np.arange(10**4)*10**-4
        bet = betainc(k+1, n-k+1, z)
        il = np.where(abs(bet - (1.-cl)/2.) == min(abs(bet - (1.-cl)/2.)))[0]
        iu = np.where(abs(bet - (1. - (1-cl)/2.)) == min(abs(bet - (1. - (1.-cl)/2.))))[0]
        p_lower=z[il]
        p_upper=z[iu]
        el = kn - p_lower
        eu = p_upper - kn
        frac[i] = kn
        lower[i] = el
        upper[i] = eu
    return frac, lower, upper