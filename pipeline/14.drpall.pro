; MaNGA SpecObj catalog
spobj = mrdfits('../0.manga_specobj.fits',1)
spobj.plateifu = strtrim(spobj.plateifu)

;;;;;;;;;;;;;;;;
; drpall cat
;---------------
drp = mrdfits('../0.drpall/drpall_clean_weight.fits',1)
drp.plateifu = strtrim(drp.plateifu)
match2,spobj.plateifu,drp.plateifu,sa,sb
; A = B[SUBA]
drp_matched = drp[sa]
help,where(sa lt 0) ; should be -1
mwrfits,drp_matched,'drpall.fits',/create
; sanity check
help,where(spobj.plateifu ne drp_matched.plateifu) ; -1

end
