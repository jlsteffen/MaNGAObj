"""
Created on Fri May 18 10:36:07 2018

The scipt is to add objects that are on the border of an IFU. Search around a 
nearby circle outside of the IFU and add objects from selected plateifus.

@author: Joshua Steffen
"""
print('\nInitializing...\n')

import numpy as np
from SciServer import SkyServer
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
import os, sys
import pandas as pd
import math as m
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import functions as f

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

drpall = fits.getdata(filepath + '/0.drpall-v3_1_1.fits')
#______________________________________________________________________________
# Creates the directory for the photometry tables
outpath = filepath + '/8.partial_overlap'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
# Make list of ra and dec for the target galaxies
ra_list = drpall['objra']
dec_list = drpall['objdec']
# Find IFU sizes to set size of the search region
ifus = drpall['ifudesignsize']
ifu_ra_list = drpall['ifura']
ifu_dec_list = drpall['ifudec']

# Text file of plateifus that have overlapping objects
plateifus = np.loadtxt(filepath + '/8.partial_overlap.txt', dtype=str)
pd.set_option('precision', 10)
#______________________________________________________________________________
for i in range(len(plateifus)):
    if not os.path.exists(outpath +'/'+ plateifus[i]+'.fits'):
        ix = np.where(drpall['plateifu']==plateifus[i])[0][0]
        
        ra_target = ra_list[ix]
        dec_target = dec_list[ix]
        ifu_ra = ifu_ra_list[ix]
        ifu_dec = ifu_dec_list[ix]
#______________________________________________________________________________
        file = fits.open(filepath + '/9.final_cat/'+plateifus[i]+'.fits')
        data = file[1].data
#______________________________________________________________________________
        # Conversion from IFU size to size on the sky in arcminutes

        IFU_size = int(drpall['ifudsgn'][ix][0:-2])
        
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        size = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]/60.0 #in arcminutes
        
        # Make an inscribed circle within the IFU
        circle_radius = size * m.cos(m.pi/6.0)/2.0    
        
        c1 = SkyCoord(ifu_ra * u.deg, ifu_dec * u.deg, frame='icrs')
#______________________________________________________________________________
        # Build artificial wcs to map point to 
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0.0, 0.0]
        w.wcs.cdelt = np.array([-size/60.0/2, size/60.0/2])
        w.wcs.crval = [ifu_ra, ifu_dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cunit = ['deg', 'deg']
        
        # Build hexagonal IFU area
        vert = np.empty((0,2))
        for j in range(6):
            vert = np.append(vert, np.array([[m.cos((j)*m.pi/3.0), m.sin((j)*m.pi/3.0)]]), axis=0)
        polygon = Polygon(vert)
#______________________________________________________________________________
        # The SQL query script
        sql1 = "SELECT TOP 100 P.ra, P.dec, P.type, P.psfMag_r, P.psfMagErr_r, P.petroMag_r, P.petroMagErr_r, P.modelMag_r, P.modelMagErr_r FROM PhotoPrimary as P JOIN dbo.fGetNearbyObjEq("+str(ra_target)+","+str(dec_target)+", "+str(size)+") AS PN ON P.objID = PN.objID ORDER BY distance"
        photo_table = SkyServer.sqlSearch(sql=sql1, dataRelease='DR14')    
#______________________________________________________________________________
        # Collect just the objects in the IFU
        plateifu = [plateifus[i]] * len(photo_table)        
        ra = photo_table.values[:,0]
        dec = photo_table.values[:,1]
    
        inside = np.empty((0), dtype = 'bool')
        for j in range(len(photo_table)):
            point = Point(w.wcs_world2pix([[ra[j], dec[j]]], 1)[0])
            if polygon.contains(point) == True:
                inside = np.append(inside, [True], axis=0)
            elif polygon.contains(point) == False:
                inside = np.append(inside, [False], axis=0)
        
        ind = []
        again = 'y'
        while again == 'y':
            print(plateifus[i])
            print(str(ra[0])+' '+str(dec[0]))
            
            targra = float(input('RA '))
            targdec = float(input('Dec '))
            
            c1 = SkyCoord(ra[~inside] * u.deg, dec[~inside] * u.deg, frame='icrs')
            c2 = SkyCoord(targra * u.deg, targdec * u.deg, frame = 'icrs')
            distance = c1.separation(c2).arcminute
            ind.append(np.where(distance == distance.min())[0][0])
            
            more = input('More Objects? ')
            if more == 'y':
                again = 'y'
            else:
                again = 'n'
        
        ind = np.array(ind)
        ra = np.append(data['ra'], ra[~inside][ind])
        dec = np.append(data['dec'], dec[~inside][ind])
        types = np.append(data['type'], photo_table.values[:,2][~inside][ind])
        psf_mag = np.append(data['psfMag_r'], photo_table.values[:,3][~inside][ind])
        psf_mag_err = np.append(data['psfMagErr_r'], photo_table.values[:,4][~inside][ind])
        petro_mag = np.append(data['petroMag_r'], photo_table.values[:,5][~inside][ind])
        petro_mag_err = np.append(data['petroMagErr_r'], photo_table.values[:,6][~inside][ind])
        model_mag = np.append(data['modelMag_r'], photo_table.values[:,7][~inside][ind])
        model_mag_err = np.append(data['modelMagErr_r'], photo_table.values[:,8][~inside][ind])
#______________________________________________________________________________    
        # Save photometric table to a .fits file
        col1 = fits.Column(name='ra', format='1D', array=ra)
        col2 = fits.Column(name='dec', format='1D', array=dec)
        col3 = fits.Column(name='type', format='1I', array=types)
        col4 = fits.Column(name='psfMag_r', format='1E', array=psf_mag)
        col5 = fits.Column(name='psfMagErr_r', format='1E', array=psf_mag_err)
        col6 = fits.Column(name='petroMag_r', format='1E', array=petro_mag)
        col7 = fits.Column(name='petroMagErr_r', format='1E', array=petro_mag_err)
        col8 = fits.Column(name='modelMag_r', format='1E', array=model_mag)
        col9 = fits.Column(name='modelMagErr_r', format='1E', array=model_mag_err)
        cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outpath +'/'+ plateifus[i]+'.fits', overwrite=True)
    
    f.update_progress((i+1.0)/np.float64(len(plateifus)))
#______________________________________________________________________________
print('\nComplete\n')