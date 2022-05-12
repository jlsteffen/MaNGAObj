PRO  SPFIT, galaxyin, errorin, log10lamin, objz, templatesin, log10lam_tpl, $
	EMISSION_SETUP_IN=emission_setup_in, BITMASK=bitmask, SRES_DATA=sres_data, $
   	DEGREE=degree, MDEGREE=mdegree, SMOMENTS=smoments, GMOMENTS=gmoments, BIAS=bias, $
	EBV_GAL=ebv_gal, REDDENING=reddening, PLOT=plot, QUIET=quiet, $ 
	FIT_RESULTS=fit_results, EMISSION_SETUP_out=emission_setup_out
;+
; NAME
;	SPFIT
;
; PURPOSE
; 	Calls pPXF twice to get BVLS-solved weights and MPFIT-solved
;	LOSVDs for the stellar continuum and the continuum-subtracted
;	(emission-line) spectrum. Then it supplies the pPXF-fit
;	parameters as initial guess to use MPFIT to solve
;	non-linearly and jointly all parameters describing the model
;	spectrum. The stellar continuum is matched with a linear 
;	combination of stellar population models, whereas emission-lines are 
; 	represented by Gauss-Hermite series, with interdependencies regulated by 
;	the input emission-setup file.
;
; INPUT
;	GALAXY 	- [N_pix] array, observed flux in f_lambda unit (1d-17 erg/s/cm2/A)
;	ERROR 	- [N_pix] array, 1-sigma error of the observed flux (1d-17 erg/s/cm2/A)
;	LOG10LAM - [N_pix] array, log10 of the observed wavelengths in observed frmae (A)
;	objz 	- constant, systemic redshift of the galaxy
;	TEMPLATES - [N_pix_tpl, N_tpl] array, L_lambda unit (erg/s/A/Msun)
;	LOG10LAM_TPL - [N_pix_tpl] array, log10 of template wavelengths
;		(rest-frame, A), must have the same spacing as LOG10LAM
;	EMISSION_SETUP_IN - structure, defined from the emission line setup file
;	BITMASK - [N_pix] array, bad pixel mask, 0: good, non-zero:bad
;	SRES_DATA - [N_pix] array, spectral resolution array (R = lambda/FWHM) of 
;		the instrument. This is used to calculate INT_DISP
;		(in km/s) array for the emission line templates (km/s) 
;
; OPTIONAL INPUT
;	DEGREE - degree of the Legendre polynomial used for correcting
;       	the template continuum shape during the fit (default: -1).
;       	This correction is ADDITIVE. Controls both pPXF and GANDALF.
;     		- IMPORTANT (GANDALF only): Additive polynomials cannot be used when
;       	the REDDENING keyword is set.
;	MDEGREE - degree of the *multiplicative* Legendre polynomial (with mean of 1)
;       	used to correct the continuum shape during the fit (default: 0). The
;       	zero degree multiplicative polynomial is always included in the fit as
;       	it corresponds to the weights assigned to the templates.
;       	Note that the computation time is longer with multiplicative
;       	polynomials than with the same number of additive polynomials.
;     		- IMPORTANT (PPXF & GANDALF): Multiplicative polynomials cannot be 
;		used when the REDDENING keyword is set.
;		- Even if not set, PPXF/GANDALF will assign MDEGREE=0 if
;		REDDENING is set
;	SMOMENTS: Order of the Gauss-Hermite moments to fit for SSPs.
;		Default is set to 4 to fit [V,sigma,h3,h4], set it to 2
;		for only [V, sigma]. 
;	GMOMENTS: Order of the Gauss-Hermite moments to fit for emission lines.
;		Default: 4
;	BIAS: This parameter biases the (h3,h4,...) measurements towards zero
;		(Gaussian LOSVD) unless their inclusion significantly decreses the
;		error in the fit. Set this to BIAS=0.0 not to bias the fit: the
;		solution (including [V,sigma]) will be noisier in that case. The
;		default BIAS should provide acceptable results in most cases, but it
;		would be safe to test it with Monte Carlo simulations. This keyword
;		precisely corresponds to the parameter \lambda in the Cappellari &
;		Emsellem (2004) paper. Note that the penalty depends on the *relative*
;		change of the fit residuals, so it is insensitive to proper scaling
;		of the NOISE vector. A nonzero BIAS can be safely used even without a
;		reliable NOISE spectrum, or with equal weighting for all pixels.
;	EBV_GAL - E(B-V) of foreground Galactic extinction to be applied before any fitting.
;  	REDDENING - Allows to include in the fit the effect of reddening by
;       	dust, by specifying a single E(B-V) guess for extinction, or a
;       	two-element array of E(B-V) guesses. A single guess will
;       	trigger the use of a single-screen dust model, affecting the
;       	entire spectrum, whereas passing a two-elements array of
;       	guesses will add a second dust component affecting only the
;       	nebular fluxes. This second option is recommended only when
;       	recombination lines are clearly detected, so that a
;       	temperature-dependent prior (e.g. the decrement of the Balmer
;       	lines) on their relative intensity can be use to constrain
;       	such a second, internal component. 
;	/PLOT - set to plot the best fitting solution at the end of pPXF fit.
;	/QUIET - set this keyword to mute the output
;
; OUTPUT
;	FIT_RESULTS - a structure contains the initial and best-fit
;		parameters of gas and stellar continuum
;	EMISSION_SETUP_OUT - structure, updated with best-fit parameters
;
;	fit_results = create_struct($
;		'i',emission_setup.i[i_l],'name',emission_setup.name[i_l],$
;		'lambda',emission_setup.lambda[i_l],$
;	        'Flux',Fgas,'Vel',Vgas,'Sigma',Sgas,'H3',H3gas,'H4',H4gas,$
;	        'AoN',AoN,'Sigma_obs',sigma_obs,$
;		'EBmV',EBmV,$		; intrinsic reddening and err
;		'apoly',apoly,'apolypars',apolypars,$ ; additive polynomial and parameters
;		'mpoly',mpoly,'mpolypars',mpolypars,$ ; multiplicative polynomial
;		'kinstar',kinstar,$	; stellar kinematics
;		'Weights',w_ssp,$       ; weights of SSP templates in units of 1e-30
;		'M_star',m_star,$       ; stellar mass
;		'chi2ppxf',chi2ppxf,$   ; Chi^2/DOF from PPXF
;		'chi2nu',chi2spfit,$    ; Chi^2/DOF from SPFIT
;		'redshift',objz,$     	; redshift
;		'EBmV_gal',ebv_gal,$    ; foreground Galactic extinction
;		'log10lam',log10lam,$   ; wavelength in rest-frame in A
;		'galaxy_in',galaxyin,$  ; input flux in Flam, observed-frame
;		'galaxy',galaxy,$ 	; Flam, rest-frame & Gal dust de-reddened 
;		'err', error,$ 		; 1-sigma errors of galaxy
;		'good',good,$ 		; mask: good pixels are 1, bad are 0.
;		'best',bestfit,$ 	; best-fit model (stellar + emission)	
;		'emis',emission,$ 	; emission-line only spectrum
;		'neat',spec_neat,$ 	; galaxy spectrum w/ emission lines removed
;		'embr',emis_br)  	; broad line only model
;
; HISTORY
;	2015/8/7 HF - Written
;	2015/8/12 - fixed bug in creating Gaussian templates for pPXF
;		emission line fit.
;		- redefined chi2ppxf using the bestfit evaluated with
;		the start_pars from PPXF
;		- decided to use FTOL = 1d-3 for MPFIT through trying
;		different values (1d-2 to 1d-6) on 8326-12701. 
;		chi^2/nu doesn't improve when FTOL is below 1d-3
;	2015/8/16 - broad lines' start_pars now takes V and S values in
;		emission_setup_in, instead of at fixed (0, 1000 km/s) 
;		values; this solves the inconsistency w/ the /BLR option. 
;		H3 & H4 are fixed at 0 for broad lines, i.e., we do not
;		fit H3 & H4 for broad lines.
;		- If degree != -1, then add the best-fit additive
;		polynomial to the templates array, and return the
;		altered array as templatesin
;	2015/8/20 - change Galactic unreddening law from Calzetti to 
;		Cardelli+89
;	2015/8/27 - added MDEGREE option to modify input templates with
;		a multiplicative polynomial after first pPXF fit to
;		stellar continuum. Save both APOLY and MPOLY to
;		fit_results structure
;		- MDEGREE = 6 significantly improves the quality of the
;		fit, but it cannot be used together with REDDENING. If
;		REDDENING is set, default MDEGREE to 0
;-

c = 299792.4580d ; Speed of light in km/s

; default parameters
IF NOT KEYWORD_SET(DEGREE) THEN degree=-1
IF KEYWORD_SET(REDDENING) OR ~KEYWORD_SET(MDEGREE) THEN mdegree=0
IF NOT KEYWORD_SET(ebv_gal) THEN ebv_gal=0.0
if ~keyword_set(smoments) then smoments = 4
if ~keyword_set(gmoments) then gmoments = 4

; rename to new variables to preserve the input data; otherwise the
; input data will be alterred upon return
galaxy = galaxyin
error = errorin
log10lam = log10lamin
templates = templatesin*1e-30 ; unit change to 1e30 erg/s/Mo/AA to
			      ; increase values of weights

; De-redden the spectra for Galactic extinction
IF ebv_gal gt 0 THEN BEGIN
    ; Cardelli et al. 1989 law
    ccm_unred,10.^log10lam,galaxy,ebv_gal,R_V=3.1
    ccm_unred,10.^log10lam,error,ebv_gal,R_V=3.1
ENDIF

; Deredshift observed galaxy spectrum 
; lambda_rf = lambda_obs/(1+z)
; flam_rf = flam_obs * (1+z) so that lambda*flam is conserved 
log10lam = log10lam - alog10(1+objz)
galaxy = galaxy * (1+objz)
error  = error * (1+objz)

; estimate stellar mass from weights
; flam_obsved (1d-17 erg/s/A/cm) = templates (erg/s/A/Msun) * Mass / fl2ll = templates * weights
; Mass = fl2ll * weight
; flam (1d-17 erg/s/A/cm2) <-> Llam (erg/s/A)
pc = 3.0856776d18 ; pc in cm
fl2ll = 1d-17 * 4 * !dpi * (lumdist(objz,/silent)*1d6*pc)^2

; Computing the velocity offset between the starting wavelengths
; of the galaxy and template spectra. You need to include a alog(10)
; factor since we are using log10 rebinned data instead of ln rebinned.
;
; The galaxy and the template spectra do not have the same starting wavelength.
; For this reason an extra velocity shift DV has to be applied to the template
; to fit the galaxy spectrum. We remove this artificial shift by using the
; keyword VSYST in the call to PPXF below, so that all velocities are
; measured with respect to DV. 
l0_gal   = log10lam[0]     ; log10(wavelength in A) 
l0_templ = log10lam_tpl[0] ; log10(wavelength in A)
tpl_offset = -(l0_gal-l0_templ)*c*alog(10.0d) 
; dV = z*c ~= ln(lambda1/lambda0)*c = log10(lambda1/lambda0)*ln(10)*c

; Preamble to PPXF
; this is the velocity scale spanned by each log_lambda interval
; You need to include a alog(10) factor as well.
lstep_gal = mean(log10lam[1:*]-log10lam) 
velscale  = c*lstep_gal*alog(10.0d)
; Initial V and sigma guesses, in km/s. For V we use the
; previously derived artificial velocity offset between 
; templates and data. For stellar velocity dispersion we 
; assume 150 km/s
start  = [tpl_offset,150.0d]

; load original input emission setup structure, to be used for masking the
; emission-line contaminated regions and fitting the emission lines in
; GANDALF. Modify emission_setup instead of emission_setup_in so that
; the latter is preserved upon return
emission_setup = emission_setup_in 

; set action='i' (ignore) for lines outside of wavelength range
npix = n_elements(galaxy)
meml_cpix = (alog10(emission_setup.lambda)-l0_gal)/lstep_gal
ind = where(meml_cpix le 0 or meml_cpix ge npix-1,ct)
if ct gt 0 then emission_setup.action[ind] = 'i'
; output emission line setup structure
; Note: do this after changing lines outside of range to action='i' 
emission_setup_out = emission_setup

; redefine emission_setup to exclude from the input structure the lines that
; are being masked/sky or ignored (action = 'm'/'s' or 'i').
i_f = where(emission_setup.action eq 'f') 
dummy = emission_setup 
emission_setup = create_struct('i',dummy.i[i_f],'name',dummy.name[i_f],$
                               'lambda',dummy.lambda[i_f],'action',dummy.action[i_f],$
                               'kind',dummy.kind[i_f],'a',dummy.a[i_f],$
                               'v',dummy.v[i_f],'s',dummy.s[i_f],$
                               'fit',dummy.fit[i_f])

; For pPXF, mask all lines with action = 'f' or 'm'
i_fit = where(emission_setup.action eq 'f') 
emission_setup.action[i_fit] = 'm' ; toggle 'f' to 'm' for pPXF
l_rf_range = 10.^minmax(log10lam)
goodpixels = mask_emission_lines(npix,0.0,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)
; mask additional bad pixels using input bitmask (0 - good, 1 - bad pixel)
if n_elements(bitmask) ne 0 then begin
	ind = where(bitmask[goodpixels] eq 0,ct)
	if ct gt 0 then goodpixels = goodpixels[ind]
endif

; PPXF fit of stellar continuum with SSPs! 
; Note: After the fit the input REDDENING is replaced with the best fitting E(B-V) value.
; so we define a new variable ebv_ppxf to avoid overwriting
nssp = (size(templates))[2]
if keyword_set(reddening) then ebv_ppxf = reddening[0] > 0.001 
PPXF, templates, galaxy, error, velscale, start, sol, error=sol_star_err,$
    	goodpixels=goodpixels, moments=smoments, degree=degree, mdegree=mdegree,$
	reddening=ebv_ppxf, lambda=10.^log10lam, bestfit=bestfit_ssp,$
	quiet=quiet, plot=plot,weights=weights_ssp,polyweights=polyweights
kinstars = sol ; [Vel* + tpl_offset,Sigma,h3,h4,h5,h6,Chi^2/DOF]
;     - When MDEGREE > 0 then SOL contains in output the 7+MDEGREE elements 
; compute best-fit Chi^2
if ~keyword_set(quiet) then print,' PPXF SSP fit, Chi^2/DOF = ', kinstars[6]

; add bestfit additive polynomial to templates array
if degree gt -1 then begin
	apolypars = polyweights ; parameters to be saved
	x = cap_range(-1d,1d,(size(templates))[1])
	apoly = 0d ; Additive polynomial
	for j=0,DEGREE do apoly += legendre(x,j)*polyWeights[j]
	templates = [[templates],[apoly]]
	weights_ssp = [weights_ssp,1.0]
	nssp = (size(templates))[2]
	templates2 = templates ; save a copy before it's modified by mpoly
endif else begin
	apolypars = -1
	apoly = -1
endelse
; multiply the bestfit multiplicative poly to templates array
if mdegree gt 0 then begin
	mpolypars = sol[7:*]
	x = cap_range(-1d,1d,(size(templates))[1])
	mpoly = 1d ; Multiplicative polynomial
	for j=1,mdegree do mpoly += legendre(x,j)*sol[6+j]
	for j=0,nssp-1 do templates[*,j] *= mpoly
	; do not multiply the additive polynomial
	if degree gt -1 then templates[*,nssp-1] /= mpoly
endif else begin
	mpolypars = 0
	mpoly = 1
endelse

; compute instrumental broadening at the wavelength of emission lines
; R = lambda/FWHM = lambda/(2.355*sigma)
; V_disp = sigma/lambda * c = c/(2.355*R)
int_disp = interpol(c/sres_data/2.355,log10lam,alog10(emission_setup.lambda)) ; km/s
int_disp_pix = int_disp/velscale ; km/s -> pixels

; For the 2nd PPXF fit, unmask emission lines by re-assign the goodpixels array
; masking what are trully action='m' (e.g., Na absorption lines)
emission_setup.action[i_fit] = 'f' ; toggle 'm' back to 'f'
goodpixels = mask_emission_lines(n_elements(galaxy),kinstars[0]-tpl_offset,emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)
; use only the pixels around each emission line (+/- 1500 km/s range) 
; to increase the Chi^2's sensitivity 
; [proven not effective - so commented out]
;badpixels = mask_emission_lines(n_elements(galaxy),sol[0]-tpl_offset,emission_setup,velscale,$
;                                 l0_gal,lstep_gal,sigma=500.0,/log10,l_rf_range=l_rf_range)
;tmp = intarr(n_elements(galaxy))+1
;tmp[badpixels] = 0
;goodpixels = where(tmp)
; mask bad pixels using input bitmask (0 - good, 1 - bad pixel)
if n_elements(bitmask) ne 0 then begin
	ind = where(bitmask[goodpixels] eq 0,ct)
	if ct gt 0 then goodpixels = goodpixels[ind]
endif
; construct emission-line templates
; index of narrow lines only 
i_nal = where(emission_setup.i lt 200) 
dummy = emission_setup 
nal_setup = create_struct('i',dummy.i[i_nal],'name',dummy.name[i_nal],$
                          'lambda',dummy.lambda[i_nal],'action',dummy.action[i_nal],$
                          'kind',dummy.kind[i_nal],'a',dummy.a[i_nal],$
                          'v',dummy.v[i_nal],'s',dummy.s[i_nal],$
                          'fit',dummy.fit[i_nal])
; index of non-doublet, narrow lines to set up parameter array
i_nal = where(emission_setup.kind eq 'l' and emission_setup.i lt 200) 
nnal  = n_elements(i_nal)
eml_pars = fltarr(nnal*5)
for i=0,nnal-1 do begin
	eml_pars[i*5] = emission_setup.a[i_nal[i]] ; total area; b/c flux (not ampl) is preserved when
			    ; pPXF convolves the templates w/ LOSVD
	eml_pars[i*5+1] = $ ; centroid in pixel
	(alog10(emission_setup.lambda[i_nal[i]])-l0_gal)/lstep_gal 
endfor
eml_templates = spfit_ghtempl(EMISSION_SETUP=nal_setup,LSTEP_GAL=lstep_gal,$
			NPIX=npix, PARS=eml_pars, INT_DISP_PIX=int_disp_pix)
; run pPXF on the residual spectrum, force DEGREE = -1, MDEGREE = 0
start = [kinstars[0]-tpl_offset, kinstars[1] < 150] ; initial guess
resi = galaxy-float(bestfit_ssp)
PPXF, eml_templates, resi, error, velscale, start, sol, error=sol_gas_err,$
    	goodpixels=goodpixels, moments=gmoments, degree=-1, mdegree=0,$
	lambda=10.^log10lam, bestfit=bestfit_gas,$
	quiet=quiet, plot=plot, weights=weights_gas
kingas = sol ; [Vel,Sigma,h3,h4,h5,h6,Chi^2/DOF]

;---------
; Prelude to SPFIT
; combine all SSP and gas best-fit parameter into one parameter array
i_l = where(emission_setup.kind eq 'l') ; now include broad lines
nlines  = n_elements(i_l)
npars = 5*nlines + 4 + nssp + 1
start_pars = dblarr(npars)
for i=0,nnal-1 do begin
	; line fluxes (integrated over pixels instead of wavelength)
	start_pars[i*5] = eml_pars[i*5] * weights_gas[i]
	; convert velocity to pixels
	start_pars[i*5+1] = eml_pars[i*5+1] + kingas[0]/velscale
	start_pars[i*5+2] = kingas[1]/velscale
	; h3, h4
	start_pars[i*5+3:i*5+4] = kingas[2:3]	
endfor
; broad lines
for i=nnal,nlines-1 do begin
	; line fluxes
	start_pars[i*5] = 0.0 
	; r.f. centroid and line width inherited from input emission_setup
	start_pars[i*5+1] = $
		(alog10(emission_setup.lambda[i_l[i]])-l0_gal)/lstep_gal+$
		emission_setup.v[i_l[i]]/velscale
	start_pars[i*5+2] = emission_setup.s[i_l[i]]/velscale
	; h3, h4
	start_pars[i*5+3:i*5+4] = [0,0]	
endfor
; SSP pars
start_pars[nlines*5:nlines*5+3] = kinstars[0:3]
start_pars[nlines*5+4:nlines*5+4+nssp-1] = weights_ssp
if keyword_set(reddening) then start_pars[nlines*5+4+nssp] = ebv_ppxf

; set good pixels array for SPFIT
goodpixels = mask_emission_lines(n_elements(galaxy),kingas[0],emission_setup,velscale,$
                                 l0_gal,lstep_gal,sigma=200.0,/log10,l_rf_range=l_rf_range)
; mask bad pixels using input bitmask (0 - good, 1 - bad pixel)
if n_elements(bitmask) ne 0 then begin
	ind = where(bitmask[goodpixels] eq 0,ct)
	if ct gt 0 then goodpixels = goodpixels[ind]
endif

; set Bias to penalize h3, h4 ...
if n_elements(bias) eq 0 then bias = 0.7d*sqrt(500d/n_elements(goodpixels)) ; pPXF paper pg.144 left

; Set up parinfo that contains the limits and the appropriate inter-dependencies
; for the parameters to be fitted. 
SPFIT_SET_CONSTRAINTS, GALAXY=galaxy, NOISE=error, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,$
	templates=templates, EMISSION_SETUP=emission_setup, L0_TEMPL=l0_templ, $
        GOODPIXELS=goodpixels, INT_DISP=int_disp, velscale=velscale, $
	START_PARS=start_pars, PARINFO=parinfo, FUNCTARGS=functargs, $
	smoments=smoments, gmoments=gmoments, bias=bias
; if reddening is not set then keep E(B-V) at zero
if ~keyword_set(reddening) then parinfo[nlines*5+4+nssp].fixed = 1

; evaluate the fit residuals before MPFIT to see how good the intial
; pars are from PPXF
resi0= spfit_fitfunc(start_pars,GALAXY=galaxy, NOISE=error, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,$
	templates=templates, EMISSION_SETUP=emission_setup, L0_TEMPL=l0_templ, $
        GOODPIXELS=goodpixels, INT_DISP=int_disp, velscale=velscale, bias=bias,$
	BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates)
chi2ppxf = robust_sigma((galaxy[goodpixels]-bestfit[goodpixels])/error[goodpixels],/zero)^2
if ~keyword_set(quiet) then print,' Initial Chi^2/DOF = ', chi2ppxf
; run MPFIT starting from previous solution
pbest = mpfit('spfit_fitfunc',start_pars, FUNCTARGS=functargs, PARINFO=parinfo, $
                    FTOL=1d-3, NFEV=ncalls, ERRMSG=errmsg, PERROR=perror, $
		    STATUS=status, /QUIET)
; evaluate the fit residuals to re-assess the fit quality and to obtain
; bestfit emission+SSP model and emission-line-only model
resid = spfit_fitfunc(pbest,GALAXY=galaxy, NOISE=error, L0_GAL=l0_gal, LSTEP_GAL=lstep_gal,$
	templates=templates, EMISSION_SETUP=emission_setup, L0_TEMPL=l0_templ, $
        GOODPIXELS=goodpixels, INT_DISP=int_disp, velscale=velscale, bias=bias,$
	BESTFIT=bestfit, EMISSION_TEMPLATES=emission_templates)
; compute reduced chi^2
chi2spfit = robust_sigma((galaxy[goodpixels]-bestfit[goodpixels])/error[goodpixels],/zero)^2
if ~keyword_set(quiet) then print,' SPFIT Chi^2/DOF = ', chi2spfit
if ~keyword_set(quiet) then print,' Function evals = ', ncalls

;-------------
; rearrange bestfit results into an easy-to-read structure
; convert centroid and sigma_pix to Vgas and Sgas in km/s
rfcenter = (alog10(emission_setup.lambda[i_l])-l0_gal)/lstep_gal
Vgas = fltarr(nlines,2) ; radial velocity (km/s)
Sgas = fltarr(nlines,2) ; intrinsic vel dispersion (km/s)
Fgas = fltarr(nlines,2) ; flux (Unit: Flam * A)
H3gas = fltarr(nlines,2); h3
H4gas = fltarr(nlines,2); h4
Agas = fltarr(nlines) ; amplitude (Flam)
for i=0,nlines-1 do begin
	Vgas[i,*] = [pbest[5*i+1]-rfcenter[i],perror[5*i+1]] * velscale
	Sgas[i,*] = [pbest[5*i+2],perror[5*i+2]] * velscale
	Agas[i] = max(emission_templates[*,i])
	; best-fit parameter gives F' = Ampl. * Sigma_pix * sqrt(2PI)
	; while we want F = Ampl. * Sigma_A * sqrt(2PI)
	; so F = F' * Sigma_A/Sigma_pix = F' * lambda_obs * velscale/c
	; and velscale = c*lstep_gal*alog(10) 
	lam_obs = 10.^(pbest[5*i+1]*lstep_gal+l0_gal) ; pix -> A
	Fgas[i,*] = [pbest[5*i],perror[5*i]] * lam_obs * lstep_gal * alog(10.0)
	; correcting flux from Gauss to true Gauss-Hermite flux
	corr = ghscorr(pbest[5*i+3],pbest[5*i+4])
	Fgas[i,*] *= corr
	; h3, h4
	H3gas[i,*] = [pbest[5*i+3],perror[5*i+3]]
	H4gas[i,*] = [pbest[5*i+4],perror[5*i+4]]
endfor
; SSP parameters
kinstar = [[pbest[5*nlines:5*nlines+3]],[perror[5*nlines:5*nlines+3]]]
kinstar[0,0] -= tpl_offset ; remove wavelength offset between SSP templates and data
w_ssp = [[pbest[5*nlines+4:5*nlines+4+nssp-1]],[perror[5*nlines+4:5*nlines+4+nssp-1]]]
; E(B-V)
EBmV = [pbest[5*nlines+4+nssp],perror[5*nlines+4+nssp]]

; observed line width in pixels
sigma_obs = sqrt(Sgas[*,0]^2. + int_disp[i_l]^2.)/velscale ; km/s -> pix
; update emission_setup_out structure
emission_setup_out.v[i_l] = Vgas[*,0]
emission_setup_out.s[i_l] = Sgas[*,0]

; Compute A/N for each line using the median noise of residual spectrum
resi = galaxy-bestfit
Vsys = mean(Vgas[*,0])
sigma = mean(sigma_obs) * velscale
Aerr = noise_emission_lines(resi,error,Vsys,emission_setup,velscale,$
                            l0_gal,lstep_gal,sigma=sigma,/log10)
AoN = Agas/Aerr[i_l]
; remove detected emission lines, with a constant A/N cut of 4.
spec_neat = galaxy
for i = 0,n_elements(AoN)-1 do $
    if (AoN[i] ge 4.0) then spec_neat -= float(emission_templates[*,i])

; total emission line only model
if ((size(emission_templates))[0] eq 1) then emission = emission_templates          
if ((size(emission_templates))[0] eq 2) then emission = total(emission_templates,2)
; good pixel mask: 0 - bad, 1 - good, opposite to bitmask
good = intarr(n_elements(galaxy)) 
good[goodpixels] = 1
; broad lines only model
i_br = where(strpos(emission_setup.name[i_l],'_br') ge 0,ct)
if ct gt 0 then emis_br = total(emission_templates[*,i_br],2) else emis_br = galaxy*0.0 

; compute stellar mass for each SSP
m_star = w_ssp*fl2ll*1d-30 
if degree gt -1 then m_star[nssp-1,*] = 0 ; last "SSP" is the polynomial

; use a structure to store best-fit parameters and errors
fit_results = create_struct($
	; best-fit emission line parameters
	'i',emission_setup.i[i_l],'name',emission_setup.name[i_l],$
	'lambda',emission_setup.lambda[i_l],$
        'Flux',Fgas,'Vel',Vgas,'Sigma',Sgas,'H3',H3gas,'H4',H4gas,$
        'AoN',AoN,'Sigma_obs',sigma_obs,$ ; emission line Amplitude-to-Noise ratio
					  ; Note: noise calculated from residual spectrum
					  ; observed line width in pixels
	'EBmV',EBmV,$		; intrinsic reddening and err
	'apoly',apoly,'apolypars',apolypars,$ ; additive polynomial and parameters
	'mpoly',mpoly,'mpolypars',mpolypars,$ ; multiplicative polynomial
	'kinstar',kinstar,$	; stellar kinematics
	'Weights',w_ssp,$       ; weights of SSP templates in units of 1e-30
	'M_star',m_star,$       ; stellar mass
	'chi2ppxf',chi2ppxf,$   ; Chi^2/DOF from PPXF
	'chi2nu',chi2spfit,$    ; Chi^2/DOF from SPFIT
	'redshift',objz,$     	; redshift
	'EBmV_gal',ebv_gal,$    ; foreground Galactic extinction
	'log10lam',log10lam,$   ; wavelength in rest-frame in log(A)
	'galaxy_in',galaxyin,$  ; input flux in Flam, observed-frame (1d-17 erg/s/cm2/A)
	'galaxy',galaxy,$ 	; Flam, rest-frame & Gal dust de-reddened 
	'err', error,$ 		; 1-sigma errors of galaxy
	'good',good,$ 		; mask: good pixels are 1, bad are 0.
	'best',bestfit,$ 	; best-fit model (stellar + emission)	
	'emis',emission,$ 	; emission-line model spectrum
	'neat',spec_neat,$ 	; galaxy spectrum w/ emission-line model removed
	'embr',emis_br)  	; broad line only model

; update templates to return to upper level program
if degree gt -1 then templatesin = templates2*1e30

END
