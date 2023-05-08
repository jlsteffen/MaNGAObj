"""
Created on Thu Aug 30 15:06:30 2018

Compare the full SDSS photoObj catalog versus the processed photometric 
catalog.

@author: joshua
"""
print('\nInitiallizing...\n')

from astropy.io import fits
import numpy as np
import sys, os
import math as m
from PIL import Image
from astropy import wcs
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from copy import copy
from glob import glob
import functions as f
import plotstyle as ps

# Filepath for this script
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
# Filepath to SDSS's pseudocolor pngs
datapath = filepath + '/0.data'
drpall = fits.getdata(filepath + '/0.drpall/drpall_clean_weight.fits')

files = glob(filepath + '/4.remov_overdebl/*.fits')

outpath = filepath + '/5.catalog_check'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
# pixel Size of the image and radius size of targeting circles
size = 300
r = 4000.0/size
scale = 40.0/size
#______________________________________________________________________________
# Information from the drpall
ras = drpall['ifura']
decs = drpall['ifudec']
ifus = drpall['ifudesignsize']
z = drpall['nsa_z']
#______________________________________________________________________________
for i in range(len(files)):
    plateifu = str(os.path.basename(files[i])).strip('.fits')
    ind = np.where(drpall['plateifu']==plateifu)[0]
    if not os.path.exists(outpath + '/' + plateifu + '1.png') and \
        os.path.exists(datapath +'/png/'+str(plateifu)+'.png') and len(ind)>0:
    
        ind = ind[0]        
    
        red_table = fits.getdata(files[i])
        obj_ra_red = red_table['ra']
        obj_dec_red = red_table['dec']
        Type_red = red_table['type']
#______________________________________________________________________________
        # Download images from skyserver    
        ra = ras[ind]
        dec = decs[ind]
        IFU_size = int(drpall['ifudsgn'][ind][0:-2])
        
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
#______________________________________________________________________________
        # Open the downloaded png for editting
        im = Image.open(datapath +'/png/'+str(plateifu)+'.png')
        im_mark = Image.open(datapath +'/png_mark/'+str(plateifu)+'.png')
        cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
        im_wide = Image.open(datapath + '/png_wide/'+str(plateifu)+'.png')
#______________________________________________________________________________
        # Make a Astropy WCS
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [cent_pix[0], cent_pix[1]]
        w.wcs.cdelt = np.array([-scale/3600.0, -scale/3600.0])
        w.wcs.crval = [ra, dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cunit = ['deg', 'deg']
        w.wcs.set_pv([(2, 1, 45.0)])
        my_dpi = 96.0
#______________________________________________________________________________
        # Set up subplots
        plt.style.use('dark_background')
        fig, ax = plt.subplots(1, 4, sharey=False, sharex=False, 
             figsize=((size)*4.0/my_dpi, size/my_dpi), dpi=my_dpi)
        ax[0].imshow(im, aspect='auto')
        ax[1].imshow(im_mark, aspect='auto')
        ax[2].imshow(im, aspect='auto')
        ax[3].imshow(im_wide, aspect='auto')
#______________________________________________________________________________
        # IFU hexagons    
        hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6, 
                  radius = size/2.0 * ifu_scale/40.0, ec = 'w', fc = 'none', orientation=m.pi/2)
        ax[1].add_patch(copy(hexagon))
        ax[2].add_patch(copy(hexagon))
        
        # Elliptical Petrosian Apertures
        crd = np.array([[obj_ra_red[0], obj_dec_red[0]]], np.float_)        
        pix = w.wcs_world2pix(crd, 1)
        obj_x = pix[0][0]
        obj_y = pix[0][1]
        rad = drpall['nsa_elpetro_th50_r'][ind] * 300/40.
        ba = drpall['nsa_elpetro_ba'][ind]
        angle = drpall['nsa_elpetro_phi'][ind]
        ellipse = patches.Ellipse((obj_x, obj_y), rad/ba, rad, angle=90-angle, ec = 'lime', fc='none')
        ax[3].add_patch(ellipse)
        
        rad = drpall['nsa_elpetro_th50_r'][ind] * 300/160.
        ba = drpall['nsa_elpetro_ba'][ind]
        ellipse = patches.Ellipse((obj_x, obj_y), rad/ba, rad, angle=90-angle, ec = 'lime', fc='none')
        ax[3].add_patch(ellipse)
#______________________________________________________________________________
        for j in range(len(obj_ra_red)):
            # Convert object ra and dec to by in pixel coordinates
            crd = np.array([[obj_ra_red[j], obj_dec_red[j]]], np.float_)        
            pix = w.wcs_world2pix(crd, 1)
            obj_x = pix[0][0]
            obj_y = pix[0][1]
            if Type_red[j] == 3: # galaxy
                ax[2].scatter(obj_x, obj_y, marker = 'o', s=250, facecolors='none', edgecolors='red')
            elif Type_red[j] == 6: # Star
                ax[2].scatter(obj_x, obj_y, marker = '*', s=250, facecolors='none', edgecolors='red')
            else: # Added object with no photometry
                ax[2].scatter(obj_x, obj_y, marker = 'P', s=250, facecolors='none', edgecolors='red')
            ax[2].text(obj_x-5, obj_y - 17, str(j), fontsize=9)
#______________________________________________________________________________
        ax[0].text(15, 30, 'RA: %.5f' % (obj_ra_red[0]), fontsize=12)
        ax[0].text(15, 50, 'Dec: %.5f'% (obj_dec_red[0]), fontsize=12)
        
        for j in range(len(ax)):
            ax[j].set_xlim(0,300)
            ax[j].set_ylim(300,0)
        
        fontsize=20
        ps.style(ax[0], fontsize, labelbottom=False, labelleft=False)
        ps.ticks(ax[0], 75, 75, xminor=15, yminor=15)
        ps.style(ax[1], fontsize, labelbottom=False, labelleft=False)
        ps.ticks(ax[1], 75, 75, xminor=15, yminor=15)
        ps.style(ax[2], fontsize, labelbottom=False, labelleft=False)
        ps.ticks(ax[2], 75, 75, xminor=15, yminor=15)
        ps.style(ax[3], fontsize, labelbottom=False, labelleft=False)
        ps.ticks(ax[3], 75, 75, xminor=15, yminor=15)
        
        ax[0].set_title(plateifu+ ' (40")', color='w')
        ax[1].set_title('PhotObj', color='w')
        ax[2].set_title('Cleaned Photometric Catalog', color='w')
        ax[3].set_title('Wide Field (160")', color='w')
        plt.subplots_adjust(wspace=0)
        plt.savefig(outpath + '/'+str(plateifu)+'.png', dpi=my_dpi, bbox_inches='tight')
        plt.close('all')
        
        im.close()
        im_mark.close()
#______________________________________________________________________________
    f.update_progress((i+1.0)/np.float64(len(files)))
print('\nComplete\n')