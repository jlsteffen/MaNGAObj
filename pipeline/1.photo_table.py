"""
Created on Fri May 18 10:36:07 2018

@author: Joshua Steffen

Downloads the photometric catalog from SDSS for the objects in MaNGA's data
reduction pipeline catalog. SDSS's database can be accessed using SciServer's
Skyserver function. Object's are queried using SQL. The script collects all 
photometric objects within twice the radius of the IFU. The photometric table 
build is organized based on distance from the plate-IFU's target object with 
the target object being the leading object in the table.
"""
print('\nInitializing...\n')

import numpy as np
from SciServer import SkyServer
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
import os, sys
import math as m
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import functions as f

# The filepath for this script
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
drpall = fits.getdata(filepath + '/0.drpall/drpall_clean_weight.fits')
#______________________________________________________________________________
# Creates the directory for the photometry tables
outpath = filepath + '/1.photo_table'
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
plateifu = drpall['plateifu']
#______________________________________________________________________________
for i in range(len(plateifu)):
    #if not os.path.exists(outpath +'/'+ str(plateifu_list[i])+'_photo.fits'):
    if ifus[i] != -9999:
        ix = np.where(drpall['plateifu']==plateifu[i])[0][0]
        
        ra_target = ra_list[ix]
        dec_target = dec_list[ix]
        ifu_ra = ifu_ra_list[ix]
        ifu_dec = ifu_dec_list[ix]
#______________________________________________________________________________
        # Conversion from IFU size to size on the sky in arcminutes
        IFU_size = int(drpall['ifudsgn'][ix][0:-2])
        
        IFU_arcsec = np.array(([19, 37, 61, 91, 127], [16.0, 21.0, 26.0, 31.0, 36.0]))
        #Size in arcminutes
        size = IFU_arcsec[1, np.where(IFU_size==IFU_arcsec[0])[0][0]]/60.0
        
        # Make an inscribed circle within the IFU
        circle_radius = size * m.cos(m.pi/6.0)/2.0    
        
        # Build coordinates for ifu center
        c1 = SkyCoord(ifu_ra * u.deg, ifu_dec * u.deg, frame='icrs')
#______________________________________________________________________________
        # Build artificial wcs to map point to 
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0.0, 0.0]
        w.wcs.cdelt = np.array([-size/60.0/2, size/60.0/2])
        w.wcs.crval = [ifu_ra, ifu_dec]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cunit = ['deg', 'deg']
        
        # Build IFU hexagon
        vert = np.empty((0,2))
        for j in range(6):
            vert = np.append(vert, np.array([[m.cos((j)*m.pi/3.0), m.sin((j)*m.pi/3.0)]]), axis=0)
        polygon = Polygon(vert)
#______________________________________________________________________________
        # The SQL query script
        sql1 = "SELECT TOP 100 P.ra, P.dec, P.type, P.psfMag_r, P.psfMagErr_r,"+\
        "P.petroMag_r, P.petroMagErr_r, P.modelMag_r, P.modelMagErr_r FROM "+\
        "PhotoPrimary as P JOIN dbo.fGetNearbyObjEq("+str(ra_target)+","+\
        str(dec_target)+", "+str(size)+") AS PN ON P.objID = PN.objID ORDER BY"+\
        " distance"
        
        photo_table = SkyServer.sqlSearch(sql=sql1, dataRelease='DR14')    
#______________________________________________________________________________
        # Collect just the objects in the IFU
        ra = photo_table.values[:,0]
        dec = photo_table.values[:,1]
    
        inside = np.empty((0), dtype = 'int')
        for j in range(len(photo_table)):
            point = Point(w.wcs_world2pix([[ra[j], dec[j]]], 1)[0])
            if polygon.contains(point) == True:
                inside = np.append(inside, [j], axis=0)
            elif polygon.contains(point) == False:
                pass
        ra = ra[inside]
        dec = dec[inside]
        types = photo_table.values[:,2][inside]
        psf_mag = photo_table.values[:,3][inside]
        psf_mag_err = photo_table.values[:,4][inside]
        petro_mag = photo_table.values[:,5][inside]
        petro_mag_err = photo_table.values[:,6][inside]
        model_mag = photo_table.values[:,7][inside]
        model_mag_err = photo_table.values[:,8][inside]
#______________________________________________________________________________    
        # Save photometric table to a .fits file
        if len(inside) > 0:    
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
            tbhdu.writeto(outpath +'/'+ str(plateifu[i])+'.fits', overwrite=True)
        else:
            # If for some reason no photometric objects are found
            print("empty - " + str(plateifu[i]))
            nan = np.array([np.nan])
            col1 = fits.Column(name='ra', format='1D', array=[ra_target])
            col2 = fits.Column(name='dec', format='1D', array=[dec_target])
            col3 = fits.Column(name='type', format='1I', array=nan)
            col4 = fits.Column(name='psfMag_r', format='1E', array=nan)
            col5 = fits.Column(name='psfMagErr_r', format='1E', array=nan)
            col6 = fits.Column(name='petroMag_r', format='1E', array=nan)
            col7 = fits.Column(name='petroMagErr_r', format='1E', array=nan)
            col8 = fits.Column(name='modelMag_r', format='1E', array=nan)
            col9 = fits.Column(name='modelMagErr_r', format='1E', array=nan)
            cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.writeto(outpath +'/'+ str(plateifu[i])+'.fits', overwrite=True)

    f.update_progress((i+1.0)/np.float64(len(plateifu)))

#______________________________________________________________________________
print('\nComplete\n')