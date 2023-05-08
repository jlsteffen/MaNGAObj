"""
Created on Wed May 16 11:27:14 2018

This script is used to manually add missing objects back into the catalogs. 
The user will be prompted to manually fill in right ascensions and declinations
for the missing objects. After filing in one missing object, the user will be
asked if there are more missing objects in the field. If yes is chosen, the 
previous query will be repeated, if no, the script will continue to the next
field.

Fields with missing objects are written in the 6.add_object.txt file.

@author: jlsteffen
"""
from astropy.io import fits
import numpy as np
import os, sys

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

# A text file list of plateifus that have missing objects
plateifus = np.loadtxt(filepath + '/6.add_object.txt', dtype=str)

outpath = filepath + '/6.add_object'
if not os.path.exists(outpath):
    os.makedirs(outpath)
#______________________________________________________________________________
def add(obj_ra, obj_dec, types, psf, petro, model, psferr, petroerr, modelerr):
    ra = input('RA [deg] = ')
    dec = input('Dec [deg] = ')
    
    obj_ra = np.append(obj_ra, [ra], axis=0)
    obj_dec = np.append(obj_dec, [dec], axis=0)
    types = np.append(types, [-99], axis=0)
    psf = np.append(psf, [-99], axis=0)
    petro = np.append(petro, [-99], axis=0)
    model = np.append(model, [-99], axis=0)
    psferr = np.append(psferr, [-99], axis=0)
    petroerr = np.append(petroerr, [-99], axis=0)
    modelerr = np.append(modelerr, [-99], axis=0) 
    return obj_ra, obj_dec, types, psf, petro, model, psferr, petroerr, modelerr
#______________________________________________________________________________
def coord(plateifu):
    """
    Add more objects to a plateifu that already had added objects.
    """    
    photo = fits.getdata(filepath + '/4.remov_overdebl/'+str(plateifu)+'.fits')
    obj_ra = photo['ra']
    obj_dec = photo['dec']
    types = photo['type']
    psf = photo['psfMag_r']
    petro = photo['petroMag_r']
    model = photo['modelMag_r']
    psferr = photo['psfMagErr_r']
    petroerr = photo['petroMagErr_r']
    modelerr = photo['modelMagErr_r']
    
    print('Target RA Dec: ' +str(photo['ra'][0])+' '+str(photo['dec'][0]))
    
    again = 'y'
    while again == 'y':
        obj_ra, obj_dec, types, psf, petro, model, psferr, petroerr, modelerr = add(obj_ra, obj_dec, types, psf, petro, model, psferr, petroerr, modelerr)
        
        more = input('More Objects? ')
        if more == 'y':
            again = 'y'
        else:
            again = 'n'
    
    # Write the FITs file
    col1 = fits.Column(name='ra', format='1D', array=obj_ra)
    col2 = fits.Column(name='dec', format='1D', array=obj_dec)
    col3 = fits.Column(name='type', format='1I', array=types)
    col4 = fits.Column(name='psfMag_r', format='1E', array=psf)
    col5 = fits.Column(name='psfMagErr_r', format='1E', array=psferr)
    col6 = fits.Column(name='petroMag_r', format='1E', array=petro)
    col7 = fits.Column(name='petroMagErr_r', format='1E', array=petroerr)
    col8 = fits.Column(name='modelMag_r', format='1E', array=model)
    col9 = fits.Column(name='modelMagErr_r', format='1E', array=modelerr)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outpath + '/'+str(plateifu)+'.fits', overwrite=True)
#______________________________________________________________________________
for i in range(len(plateifus)):
    if not os.path.exists(outpath + '/' + str(plateifus[i]) + '.fits'):
        print('\n'+plateifus[i])
        coord(plateifus[i])
