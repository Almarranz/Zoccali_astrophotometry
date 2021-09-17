PRO calibrate_brick,band


;band='Ks'
exptime=10
folder='im_jitter_NOgains/'
indir = '/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/07.1_Reduce_aligned/058_'+band+'/dit_'+strn(exptime)+'/'+folder
pruebas='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/photometry/pruebas/'
sirius='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/SIRIUS/'
tmp='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp_bs/'
results='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results_bs/'
tmp_p=pruebas
name='NPL_058'


;chip=2

delta_r = 5.
delta_mag = 5.
SIGMA_CUT = 3.0 ; sigma cut to be applied in resistant_mean
magerr_si = 0.02 ; max acceptable sirius magnitude uncertainty for photometric calibration stars

for chip=1,4 do begin

exp=readfits(tmp+'wt_chip'+strn(chip)+'.fits',ext=1)

	readcol, tmp + 'alig_SIRUS_chip'+strn(chip)+'.txt', ra,dec,f,df,x,y,dx,dy, Format ='A,A,A,A,A,A,A,A'




	dmax_sirius=1

	ra=float(ra)
	dec=float(dec)
	x=float(x)
	y=float(y)
	f=float(f)
	df=float(df)
	dx=float(dx)
	dy=float(dy)



	;~ readcol, sirius +'the_brick_idl.txt',a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si, Format = 'A,A,A,A,A,A,A,A'
	readcol, sirius +'SIRIUS_on_NPL058_c'+strn(chip)+'.txt',a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si, Format = 'A,A,A,A,A,A,A,A'
	a_si=float(a_si)
	dec_si=float(dec_si)
	J_si=float(J_si)
	dJ_si=float(dJ_si)
	H_si=float(H_si)
	dH_si=float(dH_si)
	K_si=float(K_si)
	dK_si=float(dK_si)   

	if band eq 'J' then begin
	  m_si = J_si
	  dm_si = dJ_si
	endif

	if band eq 'H' then begin
	 m_si = H_si
	 dm_si = dH_si
	endif

	if band eq 'Ks' then begin
	 m_si = K_si
	 dm_si = dK_si
	endif


	; use only valid measurements
	;bueno = where(f gt 0)
	;ra = ra[bueno]
	;dec = dec[bueno]
	;f = f[bueno]
	;df = df[bueno]
	;x=x[bueno]
	;y=y[bueno]

	; use only valid measurements
	valid = where(m_si gt 0)
	a_si = a_si[valid]
	dec_si = dec_si[valid]
	m_si = m_si[valid]
	dm_si = dm_si[valid]

	print, n_elements(m_si), band


	;####################### want to have list of valid stars for comparation after photometry calibration############################
	imagen=readfits(tmp+'wt_chip'+strn(chip)+'.fits',EXT=1,header)
	dat=readfits(tmp+'wt_chip'+strn(chip)+'.fits',EXT=0,cabeza)
	sz = size(imagen)
	xsize_hawki = sz[1]
	ysize_hawki = sz[2]
	EXTAST, header, astr

	AD2XY, a_si, dec_si, astr, x_valid, y_valid 

	; Select  stars within HAWK-I image
	; ----------------------------------------------
	cutx_min = min(x)
	cutx_max = max(x)
	cuty_min = min(y)
	cuty_max = max(y)

	good = where(x_valid gt cutx_min and x_valid lt cutx_max and y_valid gt cuty_min and y_valid lt cuty_max)
	x_valid = x_valid[good]
	y_valid = y_valid[good]
	m_valid = m_si[good]
	dm_valid = dm_si[good]

	forprint, TEXTOUT= tmp + 'VALID_SIRUS_on_' + band + '_chip' + strn(chip) + '.txt', x_valid,y_valid,m_valid,dm_valid, format='(4(e, 8X))', /NOCOMMENT



	;####################### ####################### ####################### ####################### 
	; Select calibration stars by brightness
	if band eq 'J' then begin
	  good = where(m_si gt 11.0 and dm_si lt magerr_si) 
	endif

	if band eq 'H' then begin
	;  good = where(m_si gt 11.0  and m_si lt 14 and dm_si lt 0.04)
	  good = where(m_si gt 12.0  and dm_si lt magerr_si)
	endif

	if band eq 'Ks' then begin
	  good = where(m_si gt 12.0 and dm_si lt magerr_si)    
	endif

	; appendix _si means all valid stars in SIRIUS catalogue
	; appendic _ref means potential photmetric reference stars
	a_ref = a_si[good]
	dec_ref = dec_si[good]
	m_ref = m_si[good]
	dm_ref = dm_si[good]


	imagen=readfits(tmp+'wt_chip'+strn(chip)+'.fits',EXT=1,header)
	dat=readfits(tmp+'wt_chip'+strn(chip)+'.fits',EXT=0,cabeza)
	sz = size(imagen)
	xsize_hawki = sz[1]
	ysize_hawki = sz[2]
	EXTAST, header, astr 

	; Compute pixel positions of reference stars

	AD2XY, a_ref, dec_ref, astr, x_ref, y_ref      ; selected SIRIUS reference stars
												   ; selected SIRIUS reference stars
												   
	; Select  stars within HAWK-I image
	; ----------------------------------------------
	cutx_min = min(x)
	cutx_max = max(x)
	cuty_min = min(y)
	cuty_max = max(y)

	print,'#########################'
	print,min(x),max(x),min(y),max(y)

	good = where(x_ref gt cutx_min and x_ref lt cutx_max and y_ref gt cuty_min and y_ref lt cuty_max)
	x_ref = x_ref[good]
	y_ref = y_ref[good]
	m_ref = m_ref[good]
	dm_ref = dm_ref[good]

	forprint, TEXTOUT= tmp + 'ALL_SIRUS_on_' + band + '_chip' + strn(chip) + '.txt', x_ref,y_ref,m_ref,dm_ref, format='(4(e, 4X))', /NOCOMMENT

	; Create map of all potential referenc stars for this field
	dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
	im = image_model(x_ref,y_ref,10^((-m_ref)/2.5),xsize_hawki,ysize_hawki,'gaussian', dat)
	ima = image_model(x,y,f,xsize_hawki,ysize_hawki,'gaussian', dat)
	writefits, tmp_p +'sirius_' + band + '_0.fits',dat, cabeza
	writefits, tmp_p +'sirius_' + band + '_0.fits',im, header,/app
	writefits, tmp_p +'mean_' + band + '_0.fits', ima


	dmax_sirius=1
	print, n_elements(x_ref), n_elements(x)
	compare_lists, x, y, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
	   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
	print, 'Found ' + strn(n_elements(subc1)) + ' common stars for SIRUS and HAWK-I.',strn(n_elements(sub2))



	x_common_hawki = x[subc1]
	y_common_hawki = y[subc1]
	f_common_hawki = f[subc1]
	df_common_hawki = df[subc1]
	mag_common_hawki = - 2.5*alog10(f_common_hawki) ; instrumental magnitude of common stars
	mag = -2.5*alog10(f)                            ; instrumental magnitude of all HAWKI stars

	; Find all stars common stars that are isolated
	ISOLATED_STARS, x, y, mag, x_common_hawki, y_common_hawki, mag_common_hawki, delta_mag, delta_r, ind_iso
	print, 'Found ' + strn(n_elements(ind_iso)) + ' isolated common stars in HAWK-I.'
	x_ref_hawki = x_common_hawki[ind_iso]
	y_ref_hawki = y_common_hawki[ind_iso]
	f_ref_hawki = f_common_hawki[ind_iso]
	df_ref_hawki = df_common_hawki[ind_iso]


	; Now make sure that there is no invalid pixel within +- mask pixels
	; of the potential calibration stars

	;here I can try using the mask insteand of *exp_2.fits that i dot know what is

	sz = size(wt)
	valid = []
	mask = 50
	for i=0, n_elements(x_ref_hawki)-1 do begin
	  ;print, '###',x_ref_hawki[i],y_ref_hawki[i]
	  value = total(where(exp[round(x_ref_hawki[i])-mask:round(x_ref_hawki[i])+mask, round(y_ref_hawki[i])-mask:round(y_ref_hawki[i])+mask] eq 0))
	  if value eq -1 then begin
		;print, '···',x_ref_hawki[i],y_ref_hawki[i],i
		valid = [valid, i]  
		
	  endif
	endfor
	print, n_elements(valid),n_elements(x_ref_hawki)
	print,valid



	x_ref_hawki = x_ref_hawki[valid]
	y_ref_hawki = y_ref_hawki[valid]
	f_ref_hawki = f_ref_hawki[valid]
	df_ref_hawki = df_ref_hawki[valid]

	; Find common stars with potential SIRIUS reference sources
	compare_lists, x_ref_hawki, y_ref_hawki, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
	   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
	print, 'Found ' + strn(n_elements(subc1)) + ' stars for photometric calibration.'

	x_ref_hawki = x_ref_hawki[subc1]
	y_ref_hawki = y_ref_hawki[subc1]
	f_ref_hawki = f_ref_hawki[subc1]
	df_ref_hawki = df_ref_hawki[subc1]



	; Select sources with photometric uncertainties < 0.05 in HAWK-I

	if (band eq 'J') then good = where(df_ref_hawki/f_ref_hawki lt 0.1,count) $
	  else good = where(df_ref_hawki/f_ref_hawki lt 0.05,count)

	x_ref_hawki = x_ref_hawki[good]
	y_ref_hawki = y_ref_hawki[good]
	f_ref_hawki = f_ref_hawki[good]
	df_ref_hawki = df_ref_hawki[good]

	print,'photometric uncertainties < 0.05',n_elements(x_ref_hawki)



	; Find common stars with potential SIRIUS reference sources
	compare_lists, x_ref_hawki, y_ref_hawki, x_ref, y_ref, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax_sirius,$
	   SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2 
	print, 'Found ' + strn(n_elements(subc1)) + ' stars for photometric caibration.'

	x_ref_hawki = x_ref_hawki[subc1]
	y_ref_hawki = y_ref_hawki[subc1]
	f_ref_hawki = f_ref_hawki[subc1]
	df_ref_hawki = df_ref_hawki[subc1]
	x_ref_sirius = x_ref[subc2]
	y_ref_sirius = y_ref[subc2]
	mag_ref_sirius = m_ref[subc2]
	dmag_ref_sirius = dm_ref[subc2]


	forprint, TEXTOUT= tmp + 'ref_sirius_' + band + '_chip' + strn(chip) + '.txt', x_ref_sirius ,y_ref_sirius, mag_ref_sirius, dmag_ref_sirius, format='(4(e, 4X))', /NOCOMMENT

	; Make map of calibration stars of SIRIUS
	dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
	im = image_model(x_ref_sirius,y_ref_sirius,10^(-0.4 * mag_ref_sirius),xsize_hawki,ysize_hawki,'gaussian', dat)
	;writefits, res_dir + 'phot_ref_' + band + '.fits', im
	writefits, tmp_p + 'phot_ref_' + band + '.fits', im
	
	; Make map of calibration stars of HAWKI
	dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
	im_H = image_model(x_ref_hawki,y_ref_hawki,f_ref_hawki,xsize_hawki,ysize_hawki,'gaussian', dat)
	;writefits, res_dir + 'phot_ref_' + band + '.fits', im
	writefits, tmp + 'phot_ref_' + band +'_chip'+strn(chip)+'.fits',dat, cabeza
	writefits, tmp + 'phot_ref_' + band +'_chip'+strn(chip)+'.fits', im_H,header,/app




	; Determine ZP, eliminate outliers
	zp = mag_ref_sirius + 2.5*alog10(f_ref_hawki/exptime)
	resistant_mean, zp, SIGMA_CUT, zp_mean, zp_sigma, rej, goodvec=good_cal
	forprint, TEXTOUT= tmp + 'ZP_ref_sirius_' + band + '_chip' + strn(chip) + '.txt', zp, mag_ref_sirius, format='(2(e, 4X))', /NOCOMMENT

	zp = zp_mean
	 
	print, 'ZERO POINT: ' + strn(zp) + ' +- ' + strn(zp_sigma)
	print, 'used ' + strn(n_elements(good_cal)) + ' stars for photometric calibration.'



	x_ref_hawki = x_ref_hawki[good_cal]
	y_ref_hawki = y_ref_hawki[good_cal]
	mag_ref_sirius = mag_ref_sirius[good_cal]
	f_ref_hawki = f_ref_hawki[good_cal]

	;Plot distribution of calibration stars in the field
	set_plot,'PS', /interpolate
	device, XOFFSET=0, YOFFSET=0, $
	   ;FILENAME= res_dir + 'reference_stars_photometry' + strn(chip) +'.eps', XSIZE=20., YSIZE=10., $
	   FILENAME= tmp + 'reference_stars_photometry' + strn(chip) +'.eps', XSIZE=20., YSIZE=20., $
	   /portrait, /color, BITS_PER_PIXEL=8, encapsulated=1
	!P.MULTI=0
	!P.CHARSIZE=1.2
	!X.THICK=4
	!Y.THICK=4
	!P.THICK=4.0
	!P.CHARTHICK=4
	!P.COLOR=0
	cgplot, x_ref_hawki, y_ref_hawki, Color='black', XRANGE = [0,2548], YRANGE = [0,2548], XTITLE='X-Axis [pixels]', YTITLE='Y-Axis [pixels]', XSTYLE=1, PSYM=1, YSTYLE=1
	device, /close

	;Comparison of calibrated HAWK-I reference star magnitudes with SIRIUS magnitudes to estimate the uncertainty

	m_hawki = zp - 2.5*alog10(f_ref_hawki/exptime)
	dm = m_hawki - mag_ref_sirius

	set_plot,'PS',/interpolate
	device, XOFFSET=0, YOFFSET=0, $
		FILENAME= tmp + 'ZP_calib_' + band +'_chip'+strn(chip)+ '.eps', XSIZE=20, YSIZE=15., $
		/portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
	  !P.MULTI=0
	  !P.CHARSIZE=1.3
	  !P.THICK=4.0
	  !X.THICK=4.0
	  !Y.THICK=4.0
	  !P.CHARTHICK=4.0
	  blue = cgColor("Dodger Blue")
	  red = cgColor("Dark Red")
	  green = cgColor("Green")
	  black =  cgColor("Black")
	  plot, mag_ref_sirius, dm,  XRAN = [10,21], XTITLE='[' + band + ']', YTITLE='[d' + band + ']', XSTYLE=1, PSYM=1, YSTYLE=1, YRAN=[-0.5,0.5]

	  ; plot uncertainty in bins
	  min_m = round(min(mag_ref_sirius))
	  max_m = round(max(mag_ref_sirius))
	  nbin = max_m-min_m
	  minmag = min_m
	  mags = minmag + findgen(nbin)
	  mmag = fltarr(nbin)
	  sigmag = fltarr(nbin)
	  for j = 0, nbin - 1 do begin
		thisbin = where(abs(mag_ref_sirius - mags[j]) le 0.5)
		vals = dm[thisbin]
		if n_elements(vals) gt 1 then begin
		RESISTANT_Mean,vals,3.0,Mean,Sigma,Num_Rej
		mmag[j] = Mean
		sigmag[j] = Sigma * sqrt(n_elements(thisbin) - Num_Rej - 1)
		xyouts, mags[j]-0.5, -0.3, strn(sigmag[j],FORMAT='(f5.3)'), color = blue, charsize = 0.8
		xyouts, mags[j]-0.5, -0.35, strn(n_elements(vals)-Num_Rej,FORMAT='(i)'), color = black, charsize = 1
		xyouts, mags[j]-0.5, -0.4, strn(Num_Rej,FORMAT='(i)'), color = red, charsize = 1
		endif
	  endfor
	  oploterror, mags, mmag, sigmag, errcolor = blue, errthick = 4, psym = 3
	  xyouts, 11., 0.45, 'Number of detected stars: ' + strn(n_elements(dm),format='(I)')
	  xyouts, 16, 0.3, 'ZP =  ' + strn(zp,FORMAT='(f6.3)') + ' +- ' + strn(zp_sigma,FORMAT='(f6.3)')
	  
	device, /close



	; Now, finally, compute magnitudes for all HAWK-I stars
	; ------------------------------------------------------

	m = zp - 2.5 * alog10(f/exptime)
	dm = 2.5/alog(10) * (df/f)

	set_plot,'PS',/interpolate
	device, XOFFSET=0, YOFFSET=0, $
		FILENAME = tmp_p + 'uncertainty_phot_' + strn(chip) + '.eps', XSIZE=20., YSIZE=15., $
		  /portrait, /color,BITS_PER_PIXEL=8, encapsulated=1
		!P.MULTI=0
		!P.CHARSIZE=1.3
		!P.THICK=4.0
		!P.CHARTHICK=4.0
		cgplot, m, dm,  PSYM=3, Color = 'black', XTITLE='[' + band + strn(chip) +  ']', YTITLE='[d' + band + ']',$
		  THICK = 0.5, XRANGE = [11,22], YRANGE = [0,0.1]
	device, /close 

	; Exclude stars near borders or gaps
	mask = 3
	valid = []
	for i=0, n_elements(x)-1 do begin
	  value = total(where(exp[round(x[i])-mask:round(x[i])+mask, round(y[i])-mask:round(y[i])+mask] eq 0))
	  ;value = total(where(exp[round(x[i])-mask:round(x[i])+mask, round(y[i])-mask:round(y[i])+mask] eq 0))
	  if value eq -1 then begin
	  valid = [valid, i]  
	  
	  endif
	endfor




	print,n_elements(valid)

	x = x[valid]
	y = y[valid]
	ra = ra[valid]
	dec = dec[valid]
	m = m[valid]
	dm = dm[valid]
	f = f[valid]
	df = df[valid]
	dx=dx[valid]
	dy=dy[valid]


	  ; map with only calibrated stars in HAWK-I
	dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
	im = image_model(x,y,f,xsize_hawki,ysize_hawki,'gaussian', dat)
	writefits, tmp_p +'def' + strn(chip) + '.fits', im
	   

	;Write calibrate list of stars to file
	forprint, TEXTOUT= tmp + 'stars_calibrated_' + band + '_chip' + strn(chip) + '_sirius.txt', ra ,dec , m, dm, f, df,x,y,dx,dy, format='(10(f, 4X))', /NOCOMMENT
endfor

end

