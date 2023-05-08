"""
Created on Wed Feb 24 16:41:00 2021

For plate-IFUs with offset markers for the target galaxy, go and replace the 
RA and DEC with the objra and objdec from the MaNGA DRPALL

@author: joshua
"""

from astropy.io import fits
import numpy as np
import os, sys

filepath = os.path.abspath(os.path.dirname(sys.argv[0]))
datapath = os.path.abspath(os.path.join(filepath,  '..', 'data'))

# Text file with a list of plateifus that need to be recentered
plateifus = np.loadtxt(filepath + '/7.recenter.txt', dtype=str)

drpall = fits.getdata(filepath + '/0.drpall/drpall_clean_weight.fits')

outpath = filepath + '/7.recenter'
if not os.path.exists(outpath):
    os.makedirs(outpath)
    
# these are plate where the orig center is another object
# these will keep the object as ancillary object istead of replacing the data
special = ['8565-12704', '9187-6102', '12081-9102']
    
for i in range(len(plateifus)):
    if not os.path.exists(outpath + '/' + plateifus[i] + '.fits'):
        ix = np.where(drpall['plateifu']==plateifus[i])[0]

        photo = fits.open(filepath + '/4.remov_overdebl/' + plateifus[i] +'.fits')
        data = photo[1].data
        
        if plateifus[i] not in special:
                
            data['ra'][0] = drpall['objra'][ix]
            data['dec'][0] = drpall['objdec'][ix]
            
            for j in range(len(data.names) - 2):
                data[data.names[j+2]][0] = -99
            photo.writeto(outpath + '/' + plateifus[i] + '.fits', overwrite=True)
       
        else:
            dataframe = np.asarray(data.tolist())
            
            newdata = np.zeros((1,9))
            
            newdata[:,0] = drpall['objra'][ix][0]
            newdata[:,1] = drpall['objdec'][ix][0]
            
            for j in range(len(data.names) - 2):
                newdata[:,j+2] = -99
            
            newframe = np.concatenate((newdata, dataframe), axis=0)
            
            col1 = fits.Column(name='ra', format='1D', array=newframe[:,0])
            col2 = fits.Column(name='dec', format='1D', array=newframe[:,1])
            col3 = fits.Column(name='type', format='1I', array=newframe[:,2])
            col4 = fits.Column(name='psfMag_r', format='1E', array=newframe[:,3])
            col5 = fits.Column(name='psfMagErr_r', format='1E', array=newframe[:,4])
            col6 = fits.Column(name='petroMag_r', format='1E', array=newframe[:,5])
            col7 = fits.Column(name='petroMagErr_r', format='1E', array=newframe[:,6])
            col8 = fits.Column(name='modelMag_r', format='1E', array=newframe[:,7])
            col9 = fits.Column(name='modelMagErr_r', format='1E', array=newframe[:,8])
            cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            tbhdu.writeto(outpath +'/'+ plateifus[i] +'.fits', overwrite=True)