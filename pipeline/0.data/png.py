#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 11:10:47 2018

Retrieve SDSS pseudocolor images from SDSS's databases. The code uses input 
for command line promps. 

1. It will have you give the directory to save the pngs. If the directory 
does not exit, the script will create it. 

2. It will have you request a plate-ifu to print out. You can give a single
plate-ifu, a list, or if you input "ALL", it will print all available 
plate-ifu.

3. It will ask if you want to determine the png cutout based on the
arcsecond or kiloparsec size of the region.

4. It will ask you to give the arcsecond/kpc size of the cutout. The cutout 
will be a square and the size will define the length of each of the sides.

5. It will ask if you want any optional parameters. "S" will displacy SDSS 
spectroscopic objects, "P" will display sdss photometric objects, "SP" will 
display both, and no input will give you a field without marked targets.

For my pngs I use 40" for my default thumbnails and 140 kpc for my wide-field
cutouts.

@author: joshua
"""
import os, sys, urllib
from astropy.io import fits
import numpy as np
import functions as f
from urllib.error import HTTPError
#______________________________________________________________________________
filepath = os.path.abspath(os.path.dirname(sys.argv[0]))

directory = input("Choose Directory: ")

outpath = filepath + '/' + directory
if not os.path.exists(outpath):
    os.makedirs(outpath)

drppath = os.path.abspath(os.path.join(filepath, '..', \
                          '0.drpall/originals/drpall-v3_1_1.fits'))
drpall = fits.getdata(drppath)
#______________________________________________________________________________
mngtarg1 = drpall['mngtarg1']
mngtarg3 = drpall['mngtarg3']
#______________________________________________________________________________
plateifus = drpall['plateifu']

plifus = str(input("Choose Plate-IFU(s), ALL for all Plate-IFUs: "))
if plifus == "ALL":
    pass
else:
    plateifus = [plifus]

ifuras = drpall['ifura']
ifudecs = drpall['ifudec']
z = drpall['nsa_z']
#______________________________________________________________________________

baseurl = 'https://data.sdss.org/sas/mangawork/manga/spectro/redux/MPL-7'
username = 'sdss'
password = '2.5-meters'
size = 300

a = input("Fixed arsec or kpc? [arc/kpc] ")
if a == 'arc':
    b = float(input("Arcsecond size: ")) # I use 40" for the regular cuttouts
elif a == 'kpc':
    b = float(input("Kiloparsec size: ")) # I use 140 kpc for wide frames
else:
    print("Input ERROR, choose arc or kpc")
    quit
c = input("Optional parameters? e.g. S (specObj) or P (photoObj) ")
#______________________________________________________________________________
for i in range(len(plateifus)):
    ix = np.where(drpall['plateifu']==plateifus[i])[0][0]
    plateifu = drpall['plateifu'][ix]
    if not os.path.exists(outpath + '/'+str(plateifu)+'.png'):
        try:
            ra = ifuras[ix]
            dec = ifudecs[ix]
        
            if a == 'arc':
                scale = np.float64(b/size)
            elif a == 'kpc':
                outerad = b / f.cosmo(z[ix])
                scale = np.float64(outerad/size)
            baseurl = 'https://skyserver.sdss.org/dr14/SkyServerWS/ImgCutout/getjpeg?'
            urllib.request.urlretrieve(baseurl + 'ra='+str(ra)+'&dec='+str(dec)+'&scale='+str(scale)+'&width='+str(size)+'&height='+str(size)+'&opt='+str(c),
                                       outpath +'/'+str(plateifu)+'.png')
        except HTTPError:
            pass
    else:
        pass
    f.update_progress((i+1.0)/np.float64(len(plateifus)))

#______________________________________________________________________________