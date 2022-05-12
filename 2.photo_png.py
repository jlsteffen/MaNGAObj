"""
Created on Wed May 16 11:27:14 2018

@author: jlsteffen

Prints plots for the downloaded photometric catalogs. The code will print a 
figure for every identified object, so each MaNGA field may have several pngs
printed. In the figures, the left panel is a SDSS pseudocolor image and the 
right panel is the same image with MaNGA's hexagonal field of view overlaid in
white and with the collected photObjs highlighted. The squares represnt 
galaxies and the circles represent stars as classified by SDSS's photObj table.
The white object is the MaNGA target, red object is the photObj associated with
the image, and the remaining blue objects are the other photObjs in the frame.

These images will be used to sort between real sources and over-deblended 
objects in the photObj catalog.

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
import functions as f
import plotstyle as ps
from glob import glob

plt.style.use('dark_background')

# Filepath of this script
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
# Filepath to the SDSS pseudocolor images
datapath = os.path.abspath(os.path.join(filepath, '..', 'data'))
# Filepath to the photObj tables
photopath = os.path.abspath(os.path.join(filepath, '1.photo_table'))

drpall = fits.getdata(filepath + '/0.drpall-v3_1_1.fits')

# Information from the drpall
ras = drpall['ifura']
decs = drpall['ifudec']
ifus = drpall['ifudesignsize']
#______________________________________________________________________________
outpath1 = filepath + '/2.photo_png'
if not os.path.exists(outpath1):
    os.makedirs(outpath1)
# create a identical png folder for the subsequent image sorting
outpath2 = filepath + '/3.image_sort'
if not os.path.exists(outpath2):
    os.makedirs(outpath2)
#______________________________________________________________________________
# pixel Size of the image and radius size of targeting circles
size = 300
r = 4000.0/size
#______________________________________________________________________________
def color(ind):
    '''
    Function for determining how to color the object, based on relative 
    magnitude as calculated in photo_load.py.
    '''
    if ind == 0:
        color = 'white'
    elif ind == k+1:
        color = 'red'
    else:
        color = 'blue'
    return color
#______________________________________________________________________________
# Collect all photObj table
files = glob(photopath + '/*.fits')
for i in range(len(files)):
    plateifu = str(os.path.basename(files[i])).strip('.fits')
    table = fits.getdata(files[i])
    ind = np.where(drpall['plateifu']==plateifu)[0][0]
    
    if ifus[ind] != -9999:
        obj_ra = table['ra']
        obj_dec = table['dec']
        Type = table['type']
#______________________________________________________________________________
        # Get object positions  
        ra = ras[ind]
        dec = decs[ind]
        # Get ifu position and size
        IFU_size = int(drpall['ifudsgn'][ind][0:-2])
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        ifu_scale = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]
        # arcsec to pix image scale
        scale = 40.0/size
#______________________________________________________________________________
        # Open the downloaded png for editting
        im = Image.open(datapath +'/png/'+str(plateifu)+'.png')
        cent_pix = np.array([im.size[0]/2.0, im.size[1]/2.0])
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
        for k in range(len(obj_ra)-1):
            if not os.path.exists(outpath1 + '/'+str(plateifu)+'-'+str(k+1)+'.png'):
                fig, ax = plt.subplots(1, 2, sharey=True, sharex=True, 
                     figsize=((size)*2.0/my_dpi, size/my_dpi), dpi=my_dpi)
                ax[0].imshow(im, aspect='auto')
                ax[1].imshow(im, aspect='auto')
#______________________________________________________________________________
                # IFU hexagons    
                hexagon = patches.RegularPolygon([size/2.0, size/2.0], 6, 
                          radius = size/2.0 * ifu_scale/40.0, ec = 'w', fc = 'none', orientation=m.pi/2)
                ax[1].add_patch(copy(hexagon))
#______________________________________________________________________________
                # Highlight each photObj
                for j in range(len(obj_ra)):
                    # Convert object ra and dec to by in pixel coordinates
                    crd = np.array([[obj_ra[j], obj_dec[j]]], np.float_)        
                    pix = w.wcs_world2pix(crd, 1)
                    obj_x = pix[0][0]
                    obj_y = pix[0][1]
#______________________________________________________________________________
                    # Label objects as stars, circle, or galaxies, square
                    if Type[j] == 3:
                        psf_square = patches.Rectangle([obj_x-r, obj_y-r], 2*r, 2*r, ec=color(j), fc='none')
                        ax[1].add_patch(copy(psf_square))
                    elif Type[j] == 6:
                        psf_circle = patches.Circle([obj_x, obj_y], radius=r, ec=color(j), fc='none')
                        ax[1].add_patch(copy(psf_circle))
#______________________________________________________________________________
                loc1 = 300/4.
                loc2 = 300/4./5               
                ps.style(ax, fontsize=10, labelbottom=False, labelleft=False)
                ps.ticks(ax, xmajor=loc1, ymajor=loc1, xminor=loc2, yminor=loc2)
                
                ax[0].set_title(plateifu +'-'+str(k+1), color='w')
                ax[1].set_title('Photometric Catalog', color='w')
                
                plt.subplots_adjust(wspace=0)
                plt.savefig(outpath1 + '/'+str(plateifu)+'-'+str(k+1)+'.png', dpi=my_dpi, bbox_inches='tight')
                plt.savefig(outpath2 + '/'+str(plateifu)+'-'+str(k+1)+'.png', dpi=my_dpi, bbox_inches='tight')
                
                plt.close('all')
            else:
                pass
#______________________________________________________________________________
    f.update_progress((i+1.0)/np.float64(len(files)))
    
print('\nComplete\n')