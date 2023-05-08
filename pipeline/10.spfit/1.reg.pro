;+
; generate aperture list file based on cleaned PhotoObj catalog
; this is a preparation step to run spfit on these apertures
;-

; extraction parameters
a2b = 1.0 ; major/minor axis ratio
PA  = 0.0 ; PA of major axis

FOR jj = 0,1 do begin

if jj eq 0 then begin
   outdir = 'r1arc/'
   apradi = 1.0 
   apunit = 'arc'
endif else begin
   outdir = 'r1kpc/'
   apradi = 1.3 
   apunit = 'kpc'
endelse

if ~file_test(outdir) then spawn,'mkdir '+outdir

; drpall catalog
drpall = mrdfits('../0.drpall/drpall_clean.fits',1)
drpall.plateifu = strtrim(drpall.plateifu)

; specObj catalog that incl. RA/Dec of objects inside IFU
proc = mrdfits('../0.manga_specobj.fits',1)
proc.plateifu = strtrim(proc.plateifu)

fmt = '(" REG: ",i5," out of ",i5,$,%"\r")'
nfiles = n_elements(drpall)
;kk = where(drpall.plateifu eq '10511-12705')
;for i=kk[0],kk[0] do begin
for i=0,n_elements(drpall)-1 do begin
	; skip if exist
	if file_test(outdir+drpall[i].plateifu+'.reg') then continue
		
	; calculate radius in arcsec
	if apunit eq 'arc' then begin
		rad = apradi 
	endif
	if apunit eq 'kpc' then begin
		drp = drpall[i]
		objz = float(drp.nsa_z)
		if drp.nsa_z le 0 and drp.nsa_zdist gt 0 then objz = float(drp.nsa_zdist)
		if drp.nsa_z le 0 and drp.z gt 0 then objz = float(drp.z)
		if objz le 0 then $
			message,'Invalid redshift, z = '+strc(objz)+', skip '+drpall[i].plateifu
		pltscale = dangular(objz,/kpc)/(180.*3600./!pi) ; kpc/arcsec
		rad = apradi/pltscale ; kpc->arcsec
	endif
	
	; obtain RA & Dec for components
	comps = proc[where(proc.plateifu eq drpall[i].plateifu,ct)]
	if ct eq 0 then begin 
		print,'Clean catalog not found: ',drpall[i].plateifu
		continue
	endif
	forprint,f='(2(f0.7,1x),3(f0.2,1x))',comps.ra,comps.dec,$
		rad+fltarr(ct),a2b+fltarr(ct),PA+fltarr(ct),$
		text=outdir+drpall[i].plateifu+'.reg',/nocomment
	
 	print,f=fmt,i+1,nfiles
endfor
print, '' ;; don't overwrite the final line

ENDFOR

end

