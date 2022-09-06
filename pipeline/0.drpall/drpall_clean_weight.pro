;+
; Recompute the comoving volume weights for MPL-11 MaNGA data
; 	- update cosmology to (0.3, 0.7, 70)
; 	- use #s of observed IFUs to calculate Ps (allocation probabilities)
;
; Output FITS table w/ new weights, comoving volumes, and probabilities
; 
; PWEIGHT	The volume weight for the Primary sample. 
; SWEIGHT	The volume weight for the full Secondary sample. 
; SRWEIGHT	The volume weight for the down-sampled Secondary sample. 
; EWEIGHT	The volume weight for the Primary+ sample. 
; PSRWEIGHT	The volume weight for the combined Primary and down-sampled Secondary samples. 
; ESRWEIGHT	The volume weight for the combined Primary+ and down-sampled Secondary samples. 
; PSWEIGHT	The volume weight for the combined Primary and full Secondary samples. 
; ESWEIGHT	The volume weight for the combined Primary+ and full Secondary samples. 
;
; version for MPL-11 - Mar 2021
;-

; set default cosmology (0.3, 0.7, 70)
red,/default,/verbose
use_obsved = 1 ; use observed # of IFUs to compute allocation probs

; Input unique, non-overlapping datacubes in the Main Sample
drpfile = 'drpall_clean.fits'
; output drpfile
outdrp = repstr(drpfile,'.fits','_weight.fits')
; Input tiled catalog
tilfile = '../../../cats/nsa/MaNGA_targets_extNSA_tiled_ancillary.fits'

; Parent catalog - enhanced NSA catalog 
;	641409 galaxies at z < 0.15
;	583097 with IN_DR7_LSS = 129
; Target catalog - 
; 	41154 main galaxy sample - MANGA_TARGET1 != 0
; 	17682 / 17600 / 5910 = Primary/Secondary/Color Enhanced
;	 1407 ancillary targets - MANGA_TARGET3 != 0 and MANGA_TARGET1 = 0
tiled = mrdfits(tilfile,1,/silent)
nobj = n_elements(tiled)
; # of tiles: 1800, -999 is for unallocated objects 
help,rem_dup(tiled.manga_tileid) ; 1801
; copy to a new structure
newtags = ['obsved','Vp','Vs','Ve','Vem','Pp','Pe','Psd','Psnd','Ps']
str = mrd_struct(newtags,['0',replicate('0.0',9)],nobj)
new = struct_combine(str,tiled)
; null precalculated weights
new.pweight   = -99 	; primary
new.sweight   = -99 	; secondary (full)
new.srweight  = -99 	; secondary (downsampled)
new.eweight   = -99 	; primary+
new.psrweight = -99 	; primary & secondary (downsampled)
new.esrweight = -99 	; primary+ & secondary (downsampled)
new.psweight  = -99 	; primary & secondary (full)
new.esweight  = -99 	; primary+ & secondary (full)

; load observed DRP catalog
obsvd = mrdfits(drpfile,1,/silent) ; 9574 (MPL-11) 6142 (MPL-8) 
; # of plates completed
help,obsvd,rem_dup(obsvd.plate) ; 609 (MPL-11) 381 (MPL-8) 276 (DR15) 163 (DR14)
; identify observed galaxies
match2,obsvd.NSA_NSAID,tiled.NSA_NSAID,sa,sb
idx = where(sa ge 0) 
help,idx,where(tiled[sa[idx]].NSA_NSAID ne obsvd[idx].NSA_NSAID) ; 6142, -1
new[sa[idx]].obsved = 1 

; MaNGA survey area and fiducial volume (Wake17)
area = 7362.0/(360d0^2/!DPI) ; tiled fraction of the whole sky
VF = 1e6 ; Mpc^3 - fiducial volume

; Primary sample
iP = where((tiled.MANGA_TARGET1 and 2L^10) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0,NP) ; all tiled Primary sample
iPa = where((tiled.MANGA_TARGET1 and 2L^10) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.IFUDESIGNSIZEMAIN gt 0,NPa) ; those allocated w/ IFU
iPo = where((new.MANGA_TARGET1 and 2L^10) ne 0 and $
	new.obsved ne 0,NPo) ; observed Primary sample
print,NP,NPa,NPo,1.0*NPa/NP,1.0*NPo/NP,1.*NPo/NPa
; DR14: 17063       14050        1245     0.823419    0.0729649    0.0886121
; DR15: 17063       14050        2084     0.823419     0.122136     0.148327
; DR16: 17063       13669        2943     0.801090     0.172478     0.215305

; Primary+ = Primary + Color Enhanced (E)
iE = where((tiled.MANGA_TARGET1 and 2L^10+2L^12) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0,N_E) ; all tiled sample
iEa = where((tiled.MANGA_TARGET1 and 2L^10+2L^12) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.IFUDESIGNSIZEMAIN gt 0,NEa) ; those allocated w/ IFU
iEo = where((new.MANGA_TARGET1 and 2L^10+2L^12) ne 0 and $
	new.obsved ne 0,NEo) ; observed Primary+ sample
print,N_E,NEa,NEo,1.0*NEa/N_E,1.0*NEo/N_E,1.*NEo/NEa
; DR14: 22781       18739        1683     0.822571    0.0738774    0.0898127
; DR15: 22781       18739        2783     0.822571     0.122163     0.148514
; DR16: 22781       18230        3913     0.800228     0.171766     0.214646

; down-sampled Secondary sample
iSd = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.RANFLAG eq 1,NSd) ; tiled secondaries that pass the down-sampling
iSda = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.RANFLAG eq 1 and $
	tiled.IFUDESIGNSIZEMAIN gt 0,NSda) ; those allocated w/ IFU	
iSdo = where((new.MANGA_TARGET1 and 2L^11) ne 0 and $
	new.RANFLAG eq 1 and $
	new.obsved ne 0,NSdo) ; those observed
print,NSd,NSda,NSdo,1.*NSda/NSd,1.*NSdo/NSd,1.*NSdo/NSda
; DR14: 11365        9733         807     0.856401    0.0710075    0.0829138
; DR15: 11365        9733        1406     0.856401     0.123713     0.144457
; DR16: 11365        9449        1909     0.831412     0.167972     0.202032

; non-down-sampled Secondary sample
iSnd = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.RANFLAG eq 0,NSnd) ; tiled secondaries that don't pass the down-sampling
iSnda = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.RANFLAG eq 0 and $
	tiled.IFUDESIGNSIZEMAIN gt 0,NSnda) ; those allocated w/ IFU	
iSndo = where((new.MANGA_TARGET1 and 2L^11) ne 0 and $
	new.RANFLAG eq 0 and $
	new.obsved ne 0,NSndo) ; those observed
print,NSnd,NSnda,NSndo,1.*NSnda/NSnd,1.*NSndo/NSnd,1.*NSndo/NSnda
; DR14:  5503        1438         128     0.261312    0.0232600    0.0890125 
; DR15:  5503        1438         239     0.261312    0.0434309     0.166203
; DR16:  5503        1223         324     0.222242    0.0588770     0.264922
 
; Full Secondary sample
iS = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0, NS)
iSa = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and $
	max(tiled.manga_tileids,dim=1) ge 0 and $
	tiled.IFUDESIGNSIZEMAIN gt 0,NSa) ; those allocated w/ IFU	
iSo = where((new.MANGA_TARGET1 and 2L^11) ne 0 and $
	new.obsved ne 0,NSo) ; those observed
print,NS,NSa,NSo,1.*NSa/NS,1.*NSo/NS,1.*NSo/NSa
; DR14: 16868       11171         935     0.662260    0.0554304    0.0836989
; DR15: 16868       11171        1645     0.662260    0.0975219     0.147256
; DR16: 16868       10672        2233     0.632677     0.132381     0.209239

; main galaxy sample: all Primary+ & full Secondary
; regardless tiled or not
s = where((tiled.manga_target1 and 2L^10+2L^11+2L^12) ne 0) 

; comoving volumes: Vp, Vs, Ve, Vem
Vp = fltarr(nobj)
Vs = fltarr(nobj)
Ve = fltarr(nobj)
Vp[s] = (vcomoving(tiled[s].zmax,/Mpc)-vcomoving(tiled[s].zmin,/Mpc))*area
Vs[s] = (vcomoving(tiled[s].szmax,/Mpc)-vcomoving(tiled[s].szmin,/Mpc))*area > 0
Ve[s] = (vcomoving(tiled[s].ezmax,/Mpc)-vcomoving(tiled[s].ezmin,/Mpc))*area
; modified Ve for Primary+ sample
; when there are overlapping redshift ranges between Primary+ and Secondary
Vem = Ve
idx = where(tiled[s].ezmax gt tiled[s].szmin and $
	tiled[s].szmax gt 0 and tiled[s].szmin gt 0,noverlap)
Vem[s[idx]] = (vcomoving(tiled[s[idx]].szmin,/Mpc)-$
	       vcomoving(tiled[s[idx]].ezmin,/Mpc))*area > 0

; Use observed # of IFUs to compute IFU Allocation Probabilities
if keyword_set(use_obsved) then begin
	NPa = NPo
	NEa = NEo
	NSda = NSdo
	NSnda = NSndo
endif

; IFU allocation probabilities: Pp, Pe, Psd, Psnd, Ps
Pp = 1.0*NPa/NP 
Pe = 1.0*NEa/N_E 
Ps = fltarr(nobj)
Psd = fltarr(nobj)
Psnd = fltarr(nobj)
Psd[s] = tiled[s].PROBS * NSda/NSd
Psnd[s] = (1.0-tiled[s].PROBS)*NSnda/NSnd
Ps[s] = Psd[s] + Psnd[s]

; save variables into structure
for i=1,n_elements(newtags)-1 do xx=execute('new.'+newtags[i]+'='+newtags[i])

; redefine index arrays w/o considering tiled or not
iP  = where((tiled.MANGA_TARGET1 and 2L^10) ne 0)
iS  = where((tiled.MANGA_TARGET1 and 2L^11) ne 0)
iSd = where((tiled.MANGA_TARGET1 and 2L^11) ne 0 and tiled.RANFLAG eq 1)
iE = where((tiled.MANGA_TARGET1 and 2L^10+2L^12) ne 0)

; evaluate subsample weights
new[iP].PWEIGHT = VF/(Vp[iP] * Pp)
new[iS].SWEIGHT = VF/(Vs[iS] * Ps[iS])
new[iSd].SRWEIGHT = VF/(Vs[iSd] * Psd[iSd])
new[iE].EWEIGHT = VF/(Ve[iE] * Pe)

; Primary & secondary (down-sampled)
iPSR = [iP,iSd]
new[iPSR].PSRWEIGHT = VF/(Vp[iPSR]*Pp+Vs[iPSR]*Psd[iPSR])
; Primary & secondary (Full)
iPS = [iP,iS]
new[iPS].PSWEIGHT = VF/(Vp[iPS]*Pp+Vs[iPS]*Ps[iPS])
; Primary+ & secondary (down-sampled)
iESR = [iE,iSd]
new[iESR].ESRWEIGHT = VF/(Vem[iESR]*Pe+Vs[iESR]*Psd[iESR])
; Primary+ & secondary (Full)
iES = [iE,iS]
new[iES].ESWEIGHT = VF/(Vem[iES]*Pe+Vs[iES]*Ps[iES])

; save new weights structure (matched to obsved DRPALL file)
match2,obsvd.NSA_NSAID,new.NSA_NSAID,sa,sb
help,where(obsvd.NSA_NSAID ne new[sa].NSA_NSAID) ; must be -1
new2 = new[sa]

; update input drpfile
print,' write re-calculated weights in file: '+outdrp
obsvd.pweight   = new2.pweight   ; primary
obsvd.sweight   = new2.sweight   ; secondary (full)
obsvd.srweight  = new2.srweight  ; secondary (downsampled)
obsvd.eweight   = new2.eweight   ; primary+
obsvd.psrweight = new2.psrweight ; primary & secondary (downsampled)
obsvd.esrweight = new2.esrweight ; primary+ & secondary (downsampled)
obsvd.psweight  = new2.psweight  ; primary & secondary (full)
obsvd.esweight  = new2.esweight  ; primary+ & secondary (full)
mwrfits,obsvd,outdrp,/create

end
