;+
; Combine specObj catalogs from MPL-8, MPL-11, and MPL-11 extra
; 	
;	10130 plateifus
;	15592 objects
; 	167 maintarg=1 objects have dis2obj > 1"
; 
; This version sorts the objects based on their designation:
; 	plate-ifu-index
;-

; drpall file to match results to
drpall = mrdfits('0.drpall/drpall_clean.fits',1)

; output file
outfile = '0.manga_specobj.fits'

; specObj catalogs
cat08 = mrdfits('../photoObj8/6.master/spclass.fits',1)
cat11 = mrdfits('../photoObj11/Josh/13.photo_proc_sp.fits',1)
cat1b = mrdfits('../photoObj11b/2.photo_proc_sp.fits',1)

; get unique plateifu names
pifu08 = strtrim(cat08[rem_dup(cat08.plateifu)].plateifu) ; 6142
pifu11 = strtrim(cat11[rem_dup(cat11.plateifu)].plateifu) ; 3720
pifu1b = strtrim(cat1b[rem_dup(cat1b.plateifu)].plateifu) ; 721

; total number of plateifus = 6142+3720+721-286(dupl.) = 10297
; the two extra ones compared to drpall_full: 
;	11828-1902 (baddata) 
;	8138-12703 (cube existed in MPL8 but not in MPL11)
pifu = [pifu08,pifu11,pifu1b]
pifu = pifu[rem_dup(pifu)]
help,pifu

; check for duplicates
match2,pifu08,pifu11,sa,sb
help,where(sa ge 0),where(sb ge 0) ; both -1, no overlap
match2,[pifu08,pifu11],pifu1b,sa,sb
help,where(sa ge 0),where(sb ge 0) ; 286 duplications in pifu1b (71 w/ pifu08, 215 w/ pifu11)

; remove duplicates from cat1b
idx = where(sb ge 0)
pifu_over = pifu1b[idx]
flag = intarr(n_elements(cat1b))+1
for i=0,n_elements(pifu_over)-1 do begin
	s = where(strtrim(cat1b.plateifu) eq pifu_over[i])
	flag[s] = 0
endfor
cat1b = cat1b[where(flag eq 1)]

; combine all three catalogs
comb1 = [cat08,cat11,cat1b]

; keep only plateifus in drpall
drpall.plateifu = strtrim(drpall.plateifu)
comb1.plateifu = strtrim(comb1.plateifu)
match2,drpall.plateifu,comb1.plateifu,sa,sb
s = where(sb ge 0)
comb1 = comb1[s]
; sanity check - # of plateifus must match
help,drpall,rem_dup(comb1.plateifu) ; 

; unsorted cases - 2
; turn out to be good
s = where(comb1.sptype eq -99)
help,s
comb1[s].sptype = 1

; remove tags
remove_tags,comb1,['STERN12','SCHNEIDER10'],comb2

; add tag
tags = ['spname']
vals = ['""']
toadd = mrd_struct(tags,vals,n_elements(comb2))
comb = struct_combine(comb2,toadd)

; translate spec type flags to names
restore,'../photoObj11/1.sptype_flags.sav'
sptypeDef = {spname:types,sptype:value}
for i=0,n_elements(value)-1 do begin
	s = where(comb.sptype eq value[i],ct)
	if ct ne 0 then comb[s].spname = types[i]
endfor

; rectify cases where maintarg is not index-0
idx = where(comb.maintarg eq 1 and comb.index ne 0) ; 10
forprint,comb[idx].plateifu
for i=0,n_elements(idx)-1 do begin
	; find the index-0 object
	s = where(comb.plateifu eq comb[idx[i]].plateifu and comb.index eq 0,ct)
	if ct eq 0 or ct gt 1 then stop
	for j=0,ct-1 do begin
		; swap their indices
		comb[s].index = comb[idx[i]].index
		comb[idx[i]].index = 0
	endfor
endfor
; sanity checks, both should be -1
help,where(comb.maintarg eq 1 and comb.index ne 0)
help,where(comb.index eq 0 and comb.maintarg ne 1)

; count # of objects
help,comb,where(comb.maintarg eq 1 and comb.dis2obj gt 1)

; first save to get header
mwrfits,comb,outfile,/create
hdr = headfits(outfile,exten=1)

; add comments to header
comment = [ $
'---------------------------------------',$
'Structure Tags:                        ',$  
'  PLATEIFU      - MaNGA Observation ID ',$
'  INDEX         - Obj Index in IFU     ',$
'  SPTYPE        - Spectral Type Flag   ',$
'  RA            - J2000 (deg)          ',$
'  DEC           - J2000 (deg)          ',$
'  TYPE          - SDSS PhotObj Type    ',$
'  PSFMAG_R      - SDSS PSF r mag       ',$
'  PSFMAGERR_R   - mag error            ',$
'  PETROMAG_R    - SDSS Petro r mag     ',$
'  PETROMAGERR_R - mag error            ',$
'  MODELMAG_R    - SDSS Model r mag     ',$
'  MODELMAGERR_R - mag error            ',$
'  MAINTARG      - Primary Target?      ',$
'                  Closest to ObjRA/Dec ',$
'  DIS2OBJ       - Sep w/ Primary (")   ',$
'  SPNAME        - Spectral Type Name   ',$
'---------------------------------------',$
'Spectral Type Names (SPNAME) and Flags (SPTYPE):',$
'  blagn:    4  - broad line AGN w/ correct redshift                ',$ 
'  good:     1  - good fit (correct z and enough S/N)               ',$
'  star:    -1  - prominent stellar features at z=0                 ',$ 
'  zoff:    -2  - galaxy at z different from that of the target     ',$ 
'  zoff_qso:-3  - BLAGN at z different from that of the target      ',$ 
'  border:   0  - borderline S/N spectra, unable to classify        ',$
'  lowSN:   -4  - low S/N spectra, unable to classify               ',$ 
'  others:  -9  - the photoObj is part of the main galaxy (it should',$ 
'                 have been removed when cleanning PhotoObjs)       ',$   
'               - no SPFIT spectrum (aperture lies outside          ',$ 
'                 of IFU or bin SNR < 3*0.85, so got skipped)       ',$  
'               - erroneous logcube spectra                         ',$  
'               - in any of the cases above, invalid classification ',$
'  unsorted:-99 - not yet classified, placeholder                   ',$
'---------------------------------------']
sxaddhist,comment,hdr,/comment

; sort struct based on plateifu
comb.plateifu = strtrim(comb.plateifu)
plates= long(((strsplit(comb.plateifu,'-',/ex)).ToArray())[*,0])
ifus  = long(((strsplit(comb.plateifu,'-',/ex)).ToArray())[*,1])
desig = string(plates,f='(i05)')+string(ifus,f='(i05)')+string(comb.index,f='(i02)')
idx   = sort(desig)
spObj = comb[idx]

; save to FITS file
mwrfits,spObj,outfile,hdr,/create
mwrfits,sptypeDef,outfile

end
