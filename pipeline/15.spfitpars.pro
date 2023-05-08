; combine all SPFIT files into master catalogs

; manga specobj catalog
spobj = mrdfits('../0.manga_specobj.fits',1) ; 15592 objs
pifu = strtrim(spobj[uniq(spobj.plateifu)].plateifu) ; 10130 plateifus
ct = n_elements(pifu)

; spfit directory
spdir = '../1.spfit/'

for j=0,1 do begin
  if j eq 0 then indir='r1arc' else indir='r1kpc'

  ; stack SPFIT pars using the order of plateifus in specObj catalog 
  ; so that the stacked pars match the sequence in specObj
  fmt = '(" Progress:",f6.2,"% done",$,%"\r")'
  for i=0,ct-1 do begin
  	; best-fit parameters
  	str = mrdfits(spdir+indir+'/'+pifu[i]+'.fits',1,/silent)
  	if i eq 0 then pars = str else pars = [pars,str]
  	; print progress
  	print,f=fmt,100.*i/ct
  endfor
  print,''
  help,spobj,pars ; the two arrays should have the same number of lines
  mwrfits,pars,indir+'.fits',/create
  
  ; 84 sources are too close to the IFU edges - no MaNGA spectra
  dis = sphdist(spobj.ra,spobj.dec,pars.ra,pars.dec,/deg)*3600.
  s = where(dis gt 1)
  cmt = '# PlateIFU, Index, sptype, SpecObj.RA, SpecObj.Dec, SPFIT.RA, SPFIT.Dec'
  fmt = '(a11,2(1x,i2),4(1x,f10.6))'
  forprint,strtrim(spobj[s].plateifu),spobj[s].index,spobj[s].sptype,$
  	spobj[s].ra,spobj[s].dec,pars[s].ra,pars[s].dec,$
  	text=indir+'.txt',comment=cmt,f=fmt

endfor

end
