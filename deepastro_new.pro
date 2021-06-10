PRO DEEPASTRO_new,band,exptime

;folder='im_dark/'
;folder='im_jitter_gains/'
folder='im_jitter_NOgains/'
;folder ='im_sky_ESOReflex/'

data='_NOgains'

s=500
sx=[2048,2048]

pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
outdir='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder
tmp = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp/
psf_path='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder

indir=py_pruebas
outdir=py_pruebas
tmp=py_pruebas


for chip=1, 4 do begin

readcol, indir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', x_off, y_off,x_off_s,y_off_s,Format ='A,A,A,A',COUNT=count
cube_d=count
;x_off=float(x_off)
;y_off=float(y_off)
;x_off_s=float(x_off_s)
;y_off_s=float(y_off_s)



for i=1,cube_d do begin ;###############





datos=readfits(tmp+'im'+strn(i)+'chip'+strn(chip)+'_'+band+'dit_'+strn(exptime)+'_cacho.fits',EXTEN_NO=0,header0)
im = readfits(tmp+ 'im'+strn(i)+'chip'+strn(chip)+'_'+band+'dit_'+strn(exptime)+'_cacho.fits', EXTEN_NO=1,header1)



print,'Este son dim del cacho',size(im)
	

gauss_noise_std, im, mode, std, h, v, vmean, hfit
noise = im
noise[*,*] = std
sz = size(im)
n1 = sz[1]
n2 = sz[2]
;psf = readfits('psf.fits')
psf = readfits(psf_path+'psf_im'+strn(i)+'_chip'+strn(chip)+'_.fits')
psf = psf/total(psf)   ; normalize PSF

; Settings for StarFinder you are most likely to 
; want to play with
; these parameters apply to the final StarFinder run
; not to PSF extraction
; ------------------------------------------------

sf_thresh = [5.,5.]
sf_back_box = 20
deblend = 0
deblost = 0
min_correlation = 0.7


; General settings for StarFinder
; --------------------------

correl_mag = 4
niter = 2
compbg = 1
rel_thresh = 1
guide_x = ""
guide_y = ""


; StarFinder run
;####################

  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = sf_back_box, $
        sf_thresh, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
	ESTIMATE_BG = compbg, DEBLEND = deblend, DEBLOST = deblost, $
        N_ITER = niter, SILENT=0, $
	GUIDE_X = guide_x, GUIDE_Y = guide_y, $
	SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
      	x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC

  ;writefits, 'stars.fits', stars
  ;writefits, 'bg_sf.fits', background
  ;writefits, 'resid.fits', im-stars-background
  ;subtracted =  im-stars
  ;writefits, 'subtracted.fits', im-stars
  ;writefits, 'bg.fits', median_filter(subtracted,12)
  
  writefits, outdir+ 'stars_im'+strn(i)+'_chip'+strn(chip)+'_.fits', stars
  writefits, outdir+ 'bg_sf_im'+strn(i)+'_chip'+strn(chip)+'_.fits', background
  writefits, outdir+ 'resid_im'+strn(i)+'_chip'+strn(chip)+'_.fits', im-stars-background
  subtracted =  im-stars
  writefits, outdir+'subtracted_im'+strn(i)+'_chip'+strn(chip)+'_.fits', im-stars
  writefits, outdir+'bg_deep_im'+strn(i)+'_chip'+strn(chip)+'_.fits', median_filter(subtracted,12)

  ; save list
  ; select stars in region with more than covfrac coverage
  nstars = n_elements(f)
  ;openw, outp, 'stars.txt', /get_lun
  openw, outp, outdir+'WHOLE_stars_im'+strn(i)+'_chip'+strn(chip)+'_.txt', /get_lun
  for s = 0, nstars-1 do begin
   xi = round(x[s]) & yi = round(y[s])
   printf, outp, format='(7f13.3)', x[s], y[s], f[s], sx[s], sy[s], sf[s], c[s]
  endfor
  free_lun, outp


  ;readcol, 'stars.txt', x, y, f
  readcol, outdir+'WHOLE_stars_im'+strn(i)+'_chip'+strn(chip)+'_.txt', x, y, f, Format ='A,A,A'
  x=float(x)
  y=float(y)
  f=float(f)
  dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
  im = image_model(x,y,f,n1,n2,'gaussian', dat)
  ;writefits, 'map.fits', im
  writefits, outdir+'WHOLE_map_im'+strn(i)+'_chip'+strn(chip)+'_.fits',datos,header0
  writefits, outdir+'WHOLE_map_im'+strn(i)+'_chip'+strn(chip)+'_.fits', im,header1,/app 
  
  print, '#################### DONE with image',i,',band', band,'chip',chip,' ####################'
  print, '#################### DONE with image',i,',band', band,'chip',chip,' ####################'
  print, '#################### DONE with image',i,',band', band,'chip',chip,' ####################'

endfor
endfor

  print, 'Finished.'

END
