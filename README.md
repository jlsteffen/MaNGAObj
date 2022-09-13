# MaNGAObj
This pipeline is used for the contstruction of the MaNGA object catalog (MaNGAObj) used in [Fu et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F/abstract), [Steffen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..120S/abstract), and Steffen et al. 2022a,b (in preparation). The pipeline has the following main steps;
- Cross-match MaNGA fields with SDSS's photometric catalog (PhotObj)
- Clean overdeblended sources from PhotObj
- Add in any missing sources
- Extract spectra from the MaNGA data-cubes for each PhotObj
- Build models for the extracted spectra
- Classify the extracted spectra
- Build a summary catalog

# Tables
MaNGA_SpecObj.fits contains the classifications for the identified objects. The table is organized into the following columns;

| Column | Shape | Description |
| ------ | ----- | ----------- |
| PLATEIFU | 15592, | MaNGA Observation ID |
| INDEX | 15592, | Obj Index in IFU |
| RA | 15592, | Right Ascension J2000 (deg) |
| DEC | 15592, | Declination J2000 (deg) |
| TYPE | 15592, | SDSS PhotObj Type |
| PSFMAG_R | 15592, | SDSS PSF r mag |
| PSFMAGERR_R | 15592, | mag error |
| PETROMAG_R | 15592, | SDSS Petro r mag |
| PETROMAGERR_R | 15592, | mag error |
| MODELMAG_R | 15592, | SDSS Model r mag |
| MODELMAGERR_R | 15592, | mag error |
| MAINTARG | 15592, | MaNGA Target Galaxy |
| DIS2OBJ | 15592, | Sep w/ Primary (") |
| ZCLASS | 15592, | Z classification (1=z_corr, 0=z_off, -1=unclass) |
| SCLASS | 15592, | Spectral classification (1=galaxy, 2=BLAGN, 3=star, -1=lowSN, -2=defect) |

MaNGAObj_aper2arc_v1.fits and MaNGAObj_aper2kpc_v1.fits contain positions and derived parameters for all of our identified objects. The two files are identical except that the derived parameters in MaNGAObj_aper2arc_v1.fits are calculated from a 2 arcsecond diameter circular aperture while the derived parameters in MaNGAObj_aper2kpc_v1.fits are calculated from a 2 kpc diameter circular aperture (where the aperture size is calculated from the object's redshift). 

Both data tables have the following columns;

| Column | Shape | Description |
| ------ | ----- | ----------- |
| RA | 15592, | Right Ascension J2000 (deg) |
| DEC | 15592, | Declination J2000 (deg) |
| SNR | 15592, | Signal-to-noise of the stellar continuum |
| REDSHIFT | 15592, | redshift |
| FLUX | 15592, 2, 17 | Emission line flux (value, Error), (lines) |
| EW | 15592, 17 | Equivalent width (lines) |
| VEL | 15592, 2, 17 | Gas kinematics (value, Error), (lines) |
| SIGMA | 15592, 2, 17 | Gas velocity dispersion (value, Error), (lines) |
| H3 | 15592, 2, 17 | 3rd moment of the Gauss-Hermite series (value, Error), (lines) |
| H4 | 15592, 2, 17 | 4th moment of the Gauss-Hermite series (value, Error), (lines) |
| AON | 15592, 17 | Amplitude over noise for 17 emission lines (lines) |
| SIGMA_OBS | 15592, 17 | noise calculated from residual spectrum (lines) |
| EBMV | 15592, 2 | intrinsic reddening and err (value, error) |
| KINSTAR | 15592, 2, 4 | stellar kinematics (Value, error), (velocity, velocity dispersion, ?, ?) |
| WEIGHTS | 15592, 2, 78 | weights of SSP templates in units of 1e-30 (value, error), (weights) |
| M_STAR | 15592, 2, 78 | stellar mass (value, error), (weights) |
| CHI2PPXF | 15592, | Chi^2/DOF from PPXF |
| CHI2NU | 15592, | Chi^2/DOF from SPFIT |

The extracted emission lines are given below, the index refers to its position in the above tables.

| Index | Linename | Wavelength (A) |
| ----- | -------- | -------------- |
| 0 | OII | 3730 |
| 1 | NeIII | 3870 |
| 2 | NeIII | 3969 |
| 3 | Hg | 4342 |
| 4 | OIII | 4364 |
| 5 | HeII | 4687 |
| 6 | Hb | 4863 |
| 7 | OIII | 5008 |
| 8 | NI | 5199 |
| 9 | NaI | 5892 |
| 10 | NaI | 5898 |
| 11 | OI | 6302 |
| 12 | Ha | 6565 |
| 13 | NII | 6585 |
| 14 | SII | 6718 |
| 15 | SII | 6733 |
| 16 | ArIII | 7138 |

The files can be accessed in Python using the package, [AstroPy](https://www.astropy.org/). The table may be opened and columns may be called with the following example code block.
```
from astropy.io import fits

mangaobj = fits.open('/path_to_file/MaNGA_SpecObj_v1.fits')
obj = mangaobj[1].data

mangaspfit = fits.open('/path_to_file/MaNGAObj_aper2kpc_v1.fits')
dat = mangaspfit[1].data

# Right ascensions and declinations
ra = obj['RA']
Dec = obj['DEC']

# ZCLASS and SCLASS
zclass = obj['zclass']
sclass = obj['sclass']

# plateifus with BLAGN and zclass of z_corr
ix = (zclass == 1)&(sclass == 2)
plateifus = obj['plateifu'][ix]

# H-alpha flux and error
flux = dat['flux']

halpha = flux[:,0,12]
halpha_err = flux[:,1,12]

# H-alpha flux of galaxies with zclass equal to z_corr
iy = (zclass == 1)&(sclass == 1)
halpha = flux[iy,0,12]
```

