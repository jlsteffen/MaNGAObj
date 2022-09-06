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
15.r1arc.fits and 15.r1kpc.fits contain positions, classifications, and derived parameters for all of our identified objects. The two file are identical except that the derived parameters in 15.r1arc.fits are calculated from a 1 arcsecond circular aperture while the derived parameters in 15.r1kpc.fits are calculated from a 1 kpc circular aperture (where the aperture size is calculated from the object's redshift). 
