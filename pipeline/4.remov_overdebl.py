"""
Created on Tue Aug 21 11:01:35 2018

This script creates updated photometric tables for MaNGA targets where 
over-deblended sources have been removed.

@author: joshua
"""
print('\nInitializing ...\n')

import os, sys
from glob import glob
import numpy as np
from astropy.io import fits
import functions as f
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
photopath = os.path.abspath(os.path.join(filepath,  '1.photo_table'))
overpath = os.path.abspath(os.path.join(filepath,  '3.image_sort'))

drpall = fits.getdata(filepath + '/0.drpall/drpall_clean_weight.fits')

files = glob(overpath + '/true/*.png')
photo_files = glob(photopath + '/*.fits')
allplateifus = np.empty((0))
for i in photo_files:
    basename = os.path.basename(i).strip('.fits').split('-')
    allplateifus = np.append(allplateifus, [str(basename[0]) + '-' + str(basename[1])], axis=0)

outpath = filepath + '/4.remov_overdebl/'
if not os.path.exists(outpath):
        os.makedirs(outpath)

plateifus = np.empty((0))
objnums = np.empty((0), dtype='int')
for i in files:
    basename = os.path.basename(i).strip('.png').split('-')
    plateifus = np.append(plateifus, [str(basename[0]) + '-' + str(basename[1])], axis=0)
    objnums = np.append(objnums, [int(basename[2])], axis=0)
#______________________________________________________________________________
for i in range(len(allplateifus)):
    ind = np.where(plateifus == allplateifus[i])[0]
    objnum = np.append(np.array([0], dtype='int'), objnums[ind], axis=0)
    
    photometry = fits.getdata(photopath + '/'+str(allplateifus[i])+'.fits')
    
    if np.isnan(photometry['ra'][0]) == False:
        # Table is there is data        
        new_table = photometry[objnum]
        
        col1 = fits.Column(name='ra', format='1D', array=new_table['ra'])
        col2 = fits.Column(name='dec', format='1D', array=new_table['dec'])
        col3 = fits.Column(name='type', format='1I', array=new_table['type'])
        col4 = fits.Column(name='psfMag_r', format='1E', array=new_table['psfMag_r'])
        col5 = fits.Column(name='psfMagErr_r', format='1E', array=new_table['psfMagErr_r'])
        col6 = fits.Column(name='petroMag_r', format='1E', array=new_table['petroMag_r'])
        col7 = fits.Column(name='petroMagErr_r', format='1E', array=new_table['petroMagErr_r'])
        col8 = fits.Column(name='modelMag_r', format='1E', array=new_table['modelMag_r'])
        col9 = fits.Column(name='modelMagErr_r', format='1E', array=new_table['modelMagErr_r'])
        cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outpath+str(allplateifus[i])+'.fits', overwrite=True)
    else:
        # Table for handling empty files
        ind = np.where(allplateifus[i] == drpall['plateifu'])[0]
        col1 = fits.Column(name='ra', format='1D', array=[drpall['objra'][ind]])
        col2 = fits.Column(name='dec', format='1D', array=[drpall['objdec'][ind]])
        col3 = fits.Column(name='type', format='1I', array=[-99])
        col4 = fits.Column(name='psfMag_r', format='1E', array=[-99])
        col5 = fits.Column(name='psfMagErr_r', format='1E', array=[-99])
        col6 = fits.Column(name='petroMag_r', format='1E', array=[-99])
        col7 = fits.Column(name='petroMagErr_r', format='1E', array=[-99])
        col8 = fits.Column(name='modelMag_r', format='1E', array=[-99])
        col9 = fits.Column(name='modelMagErr_r', format='1E', array=[-99])
        cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(outpath+str(allplateifus[i])+'.fits', overwrite=True)
#______________________________________________________________________________
    f.update_progress((i+1.0)/np.float64(len(allplateifus)))
    
print('\nComplete\n')