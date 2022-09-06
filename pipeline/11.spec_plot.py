#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:08:23 2021

Print the extracted and modelled spectra.

@author: joshua
"""

print('\nInitiallizing...\n')

from astropy.io import fits
import numpy as np
import sys, os
from glob import glob
import functions as f
import matplotlib.pyplot as plt
from PIL import Image
from astropy import wcs
from astropy import units as u
import matplotlib.patches as patches
from astropy.coordinates import SkyCoord
import math as m
import plotstyle as ps
from copy import copy

import warnings
warnings.filterwarnings("ignore")

def red(wave, z):
    'basic redshift function'
    return wave/(1+z)

# Filepath for this script
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
# Filepath for the SDSS pseudocolor images
datapath = os.path.abspath(os.path.join(filepath, '..', 'data'))

drpall = fits.getdata(filepath + '/0.drpall/drpall_clean_weight.fits')
size = 300
rad = 4000.0/size
my_dpi = 96.0
scale = 40.0/size

ras = drpall['ifura']
decs = drpall['ifudec']
tarra = drpall['objra']
tardec = drpall['objdec']
ifusize = drpall['ifudesignsize']
redshift = drpall['nsa_z']

plt.style.use('dark_background')

outpath = filepath + '/11.spec_plot/' 
if not os.path.exists(outpath):
    os.makedirs(outpath)
    os.makedirs(outpath + '/blagn')
    os.makedirs(outpath + '/border')
    os.makedirs(outpath + '/good')
    os.makedirs(outpath + '/low_snr')
    os.makedirs(outpath + '/others')
    os.makedirs(outpath + '/star')
    os.makedirs(outpath + '/z_off')
    os.makedirs(outpath + '/z_off_qso')
#______________________________________________________________________________
# data from emlines.dat
emlines = open(datapath + '/emlines.dat').readlines()
emission_wave = np.empty((0))
emission_name = np.empty((0))
for i in emlines:
    emission_wave = np.append(emission_wave, [float(i.split(' ')[0])], axis=0)
    emission_name = np.append(emission_name, [i.split(' ')[1]], axis=0)

# data from ellipticalabsorptionlines.dat
ellipticalabsorption = open(datapath + '/ellipticalabsorptionlines.dat').readlines()
absorption_wave = np.empty((0))
absorption_name = np.empty((0))
for i in ellipticalabsorption:
    if '-' not in i:
        absorption_wave = np.append(absorption_wave, [float(i.split(' ')[0])], axis=0)
        absorption_name = np.append(absorption_name, [i.split(' ')[1]], axis=0)
    else:
        pass

# select a subsection of the absorption lines
abselect = [2,3,4,5,10,17,40,42,43,45,46]
absorption_wave = absorption_wave[abselect] 

# Collect spfit files
files = glob(filepath + '/10.spfit/*_spec.fits')

# Make txt file listing all empty spectra
txtpath = filepath + '/11.empty.txt'
txtpath2 = filepath + '/11.printed.txt'
if not os.path.exists(txtpath):
    empty = np.array([])
else:
    empty = np.loadtxt(txtpath, dtype='str')
if not os.path.exists(txtpath2):
    printed = np.array([])
else:
    printed = np.loadtxt(txtpath2, dtype='str')

plt.ioff()
#______________________________________________________________________________
for i in range(len(files)):
    plateifu = os.path.basename(files[i]).strip('_spec.fits')
    
    with Image.open(datapath +'/png/'+str(plateifu)+'.png') as im:
        file = fits.getdata(filepath + '/10.spfit/'+plateifu+'_spec.fits')
        cat = fits.getdata(filepath + '/9.final_table/'+plateifu+'.fits')
        
        cra = cat['ra']
        cdec = cat['dec']
        ctype = cat['type']
#______________________________________________________________________________
        # Select IFU parameters for the observation        
        ind = np.where(drpall['plateifu']==plateifu)[0][0]
        ra = ras[ind]
        dec = decs[ind]
        IFU_size = ifusize[ind]
        z = redshift[ind]
        
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
#______________________________________________________________________________
        # Make a Astropy WCS for Skyserver image
        cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
        imw = wcs.WCS(naxis=2)
        imw.wcs.crpix = [cent_pix[0], cent_pix[1]]
        imw.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
        imw.wcs.crval = [ra, dec]
        imw.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        imw.wcs.cunit = ['deg', 'deg']
        imw.wcs.set_pv([(2, 1, 45.0)])
#______________________________________________________________________________
        # Find which Object is the target object
        c1 = SkyCoord(cra * u.deg, cdec * u.deg, frame = 'icrs')
        c2 = SkyCoord(tarra[ind] * u.deg, tardec[ind] * u.deg, frame = 'icrs')
        sep = c1.separation(c2).arcsecond
        target = np.where(sep == sep.min())[0][0]
        
        tarcard = np.array([[cra[target], cdec[target]]], np.float_)
        im_target = imw.wcs_world2pix(tarcard, 1)[0]
#______________________________________________________________________________
        for j in range(len(cat)):
            if plateifu+'1-'+str(j) not in empty and plateifu+'1-'+str(j) not in printed:            
                wave, flux, model, emlines = 10**file['LOG10LAM'][j], \
                    file['GALAXY'][j], file['BEST'][j], file['EMIS'][j]
                twave, tflux, tmodel, temlines = 10**file['LOG10LAM'][0], \
                    file['GALAXY'][0], file['BEST'][0], file['EMIS'][0]

                
                if not all(np.isnan(flux)):
                    
                    wave = red(wave, drpall['nsa_z'][drpall['plateifu']==plateifu][0])
                    twave = red(twave, drpall['nsa_z'][drpall['plateifu']==plateifu][0])
                    
                    fontsize = 20
                    fig = plt.figure(figsize=(20,8))
                    gs = fig.add_gridspec(2, 2, height_ratios=(1,1), width_ratios=(10, 3),
                                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                                          wspace=0.0, hspace=0.0)
                    
                    ax1 = fig.add_subplot(gs[0, 0])
                    ax2 = fig.add_subplot(gs[0, 1])
                    ax3 = fig.add_subplot(gs[1, 1])
                    ax4 = fig.add_subplot(gs[1, 0])

                    ax1.plot(twave, tflux, lw=0.5, c='w')
                    ax1.plot(twave, tmodel, lw=0.75, c='r', linestyle='--')
                    ax1.plot(twave, temlines, lw=0.5, c='cyan')
                    
                    ax4.plot(wave, flux, lw=0.5, c='w')
                    ax4.plot(wave, model, lw=0.75, c='r', linestyle='--')
                    ax4.plot(wave, emlines, lw=0.5, c='cyan')
                    
                    for k in range(len(emission_wave)):
                        ax1.axvline(emission_wave[k], linestyle=':', color='lime', lw=0.75)
                        ax4.axvline(emission_wave[k], linestyle=':', color='lime', lw=0.75)
                    for k in range(len(absorption_wave)):
                        ax1.axvline(absorption_wave[k], linestyle=':', color='cyan', lw=0.75)
                        ax4.axvline(absorption_wave[k], linestyle=':', color='cyan', lw=0.75)
                    
                    # Highlight regions with significant sky lines
                    ax1.axvspan(6775, 6875, color = 'yellow', alpha=0.25)
                    ax1.axvspan(6950, 7050, color = 'yellow', alpha=0.25)
                    ax1.axvspan(7575, 7700, color = 'yellow', alpha=0.25)
                    ax4.axvspan(6775, 6875, color = 'yellow', alpha=0.25)
                    ax4.axvspan(6950, 7050, color = 'yellow', alpha=0.25)
                    ax4.axvspan(7575, 7700, color = 'yellow', alpha=0.25)
                    
                    ax1.set_xlim(3500, 7500)
                    ax4.set_xlim(3500, 7500)
                    
                    
                    if not all(np.isnan(tflux)):
                        ax1.set_ylim(-0.01, np.nanmax(tflux)*1.1)
                        ax1.set_ylim(-0.01, np.nanmean(tflux)*2.5)
                    if not all (np.isnan(flux)):
                        ax4.set_ylim(-0.01, np.nanmax(flux)*1.1)
                    
                    ax4.set_ylim(-0.01, np.nanmean(flux)*2.5)
                    
                    ax1.set_ylabel('Target Flux', fontsize=fontsize)
                    ax4.set_ylabel('Object Flux', fontsize=fontsize)
                    ax4.set_xlabel(r'Wavelength [$\AA$]', fontsize=fontsize)
#______________________________________________________________________________
                    # Label all object in Frame
                    ax2.imshow(im)
                    ax3.imshow(im)
                    hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6, radius = size/2.0 * ifu_scale/40.0, ec = 'w', fc = 'none', orientation=m.pi/2)
                    ax2.add_patch(copy(hexagon))
                    ax3.add_patch(copy(hexagon))
                    
                    for k in range(len(ctype)):
                        card = np.array([[cra[k], cdec[k]]], np.float_)
                        obj_x, obj_y = imw.wcs_world2pix(card, 1)[0]
                        if ctype[k] == 3: # galaxy
                            ax2.scatter(obj_x, obj_y, marker = 'o', s=250, facecolors='none', edgecolors='lime')
                        elif ctype[k] == 6: # Star
                            ax2.scatter(obj_x, obj_y, marker = '*', s=250, facecolors='none', edgecolors='lime')
                        else: # Added object with no photometry
                            ax2.scatter(obj_x, obj_y, marker = 'P', s=250, facecolors='none', edgecolors='lime')
                        ax2.text(obj_x-5, obj_y - 17, str(k), fontsize=9)
                    
                    card = np.array([[cra[j], cdec[j]]], np.float_)
                    obj_x, obj_y = imw.wcs_world2pix(card, 1)[0]
                    ax2.scatter(obj_x, obj_y, marker='x', s=250, edgecolors='cyan')
#______________________________________________________________________________            
                    plt.suptitle(plateifu+'-'+str(j), fontsize=fontsize*1.5)
                    
                    ps.style(ax1, fontsize, labelbottom=False)
                    ps.style(ax4, fontsize)
                    ps.style(ax2, fontsize, labelbottom=False, labelleft=False)
                    ps.ticks(ax2, 80, 80, xminor=16, yminor=16)
                    ps.style(ax3, fontsize, labelbottom=False, labelleft=False)
                    ps.ticks(ax3, 80, 80, xminor=16, yminor=16)

                    
                    plt.savefig(outpath+'/'+plateifu+'-'+str(j)+'.png', bbox_inches='tight')
                    plt.close(fig)
                    
                    printed = np.append(printed, plateifu+'-'+str(j))
                    np.savetxt(txtpath2, printed, fmt='%s')
                else:
                    empty = np.append(empty, plateifu+'-'+str(j))
                    np.savetxt(txtpath, empty, fmt='%s')
    f.update_progress((i+1.0)/np.float64(len(files)))