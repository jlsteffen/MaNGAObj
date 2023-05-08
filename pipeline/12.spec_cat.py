#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 16:34:51 2021

Make table summarizing the spectral classifications of the identified objects.

	* blagn/4   - broad line AGN w/ correct redshift
	* good/1    - good fit (correct z and enough S/N), 
	              regardless of index
	* border/0  - low S/N spectra, cannot decide if the fit is good
	* star/-1   - prominent stellar features at z=0
	* zoff/-2   - galaxy at z different from that of the target
	* zoff_qso/-3 - BLAGN at z different from that of the target
	* lowSN/-4  - poor S/N spectra, unable to measure redshift
	* others/-9 - the photoObj is part of the main galaxy (it should
	  	      have been removed when cleanning PhotoObjs)
		    - no SPFIT spectrum (aperture lies outside
		      of IFU or bin SNR < 3*0.85, so got skipped)
		    - in either case above, the object should be
		      considered as invalid. 

@author: joshua
"""

print('\nInitiallizing...\n')

from astropy.io import fits
import numpy as np
import sys, os
from glob import glob
import functions as f


filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

proc = fits.open(filepath + '/12.photo_proc.fits')
data = proc[1].data

blagn = glob(filepath + '/12.spec_plot/blagn/*.png')
border = glob(filepath + '/12.spec_plot/border/*.png')
good = glob(filepath + '/12.spec_plot/good/*.png')
lowsnr = glob(filepath + '/12.spec_plot/low_snr/*.png')
others = glob(filepath + '/12.spec_plot/others/*.png')
star = glob(filepath + '/12.spec_plot/star/*.png')
zoff = glob(filepath + '/12.spec_plot/z_off/*.png')
zoffqso = glob(filepath + '/12.spec_plot/z_off_qso/*.png')

empty = np.loadtxt(filepath + '/12.empty.txt', dtype='str')

sptype = np.array([-99]*len(data))
for i in range(len(data)):
    obj = data['plateifu'][i]+'-'+str(data['index'][i])
    
    if os.path.exists(filepath + '/12.spec_plot/blagn/'+obj+'.png'):
        sptype[i] = 4
    elif os.path.exists(filepath + '/12.spec_plot/border/'+obj+'.png'):
        sptype[i] = 0
    elif os.path.exists(filepath + '/12.spec_plot/good/'+obj+'.png'):
        sptype[i] = 1
    elif os.path.exists(filepath + '/12.spec_plot/low_snr/'+obj+'.png'):
        sptype[i] = -4
    elif os.path.exists(filepath + '/12.spec_plot/star/'+obj+'.png'):
        sptype[i] = -1
    elif os.path.exists(filepath + '/12.spec_plot/z_off/'+obj+'.png'):
        sptype[i] = -2
    elif os.path.exists(filepath + '/12.spec_plot/z_off_qso/'+obj+'.png'):
        sptype[i] = -3
    elif os.path.exists(filepath + '/12.spec_plot/others/'+obj+'.png'):
        sptype[i] = -9
    elif obj in empty:
        sptype[i] = -9
    else:
        print(obj)
        
    f.update_progress((i+1.0)/np.float64(len(data)))
    
data['sptype'] = sptype
proc.writeto(filepath + '/12.photo_proc_sp.fits', overwrite=True)

print('\nComplete\n')