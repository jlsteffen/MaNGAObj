;+
; Run fit_manga with DS9 region files defining the extraction regions
; reg1arc/ % Time elapsed: 55756.613 seconds.
; reg1kpc/ % 
;-

; all plateifus from drpall
drpall = mrdfits('../0.drpall/drpall_clean.fits',1)
plateifu = strtrim(drpall.plateifu)

; string -> long
plates   = long(((strsplit(plateifu,'-',/ex)).ToArray())[*,0])
ifus     = long(((strsplit(plateifu,'-',/ex)).ToArray())[*,1])

for j=0,1 do begin
  if j eq 0 then indir='r1arc/' else indir='r1kpc/'
  ; run fit_manga in parallel mode
  cmd = 'fit_manga'
  extra = ',tag=''MPL-11'',mdegree=6,ds9reg='''+indir+$
  	''',/ignore_drp3qual,/overwrite,prebin=1,/quiet'
  outputs=indir+strc(plates)+'-'+strc(ifus)+'.fits' 
  manga_parallel,plates,ifus,cmd,extra,ncpu=12,outputs=outputs,/skip
endfor

end
