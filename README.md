# MaNGAObj
This pipeline is used for the contstruction of the MaNGA object catalog (MaNGAObj) used in [Fu et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJ...856...93F/abstract), [Steffen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..120S/abstract), and Steffen et al. 2022a,b (in preparation). The pipeline has the following main steps;
- Cross-match MaNGA fields with SDSS's photometric catalog (PhotObj)
- Clean overdeblended sources from PhotObj
- Add in any missing sources
- Extract spectra from the MaNGA data-cubes for each PhotObj
- Build models for the extracted spectra
- Classify the extracted spectra
- Build a summary catalog
