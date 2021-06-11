PRO EXTRACTPSF_new,band,exptime


;band='H'
;exptime=10
folder='im_jitter_NOgains/'
;folder='im_sky_ESOReflex/'

pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
outdir='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder

tmpdir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp/


for chip=1,4 do begin
;for i=imagen,imagen do begin ;###############
readcol, indir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', x_off, y_off,x_off_s,y_off_s,Format ='A,A,A,A',COUNT=count
cube_d=fix(count)
print, 'Number of images',count

for i=1, cube_d do begin ;###############

cube = readfits(tmpdir+'im'+strn(i)+'chip'+strn(chip)+'_'+band+'dit_'+strn(exptime)+'_cacho.fits', EXTEN_NO=1)
l=size(cube)
lado=l[1]
im=cube
;writefits, pruebas+'cacho_im'+strn(imagen)+band+'.fits',im
;stop
print,'Este son dim del cacho',size(im)

;tmpdir = 'tmp/'

if band eq 'H' then begin
ZP = 26.33
endif else begin
ZP=25.66
endelse
print,'###########################################'
print, 'Estas usando ZP=',ZP,' para la banda ', band
print,'###########################################'


;im = readfits('../Reduced/im_minus_sky/H/dit_2/im1_chip1.fits',EXT=1)
gauss_noise_std, im, mode, std, h, v, vmean, hfit
noise = im
noise[*,*] = std
sz = size(im)
n1 = sz[1]
n2 = sz[2]

support = fltarr(n1,n2)
valid = where(FINITE(im),complement=invalid)
support[invalid] = 0
support[valid] = 1
;writefits, tmpdir + 'support.fits', support

; Parameters that need to be edited frequently
min_correlation = 0.8
n_ref_min = 5
n_ref_max = 100
maskrad = 24
satlevel = 1.1e5
unweighted = 1
DIT = exptime
mag_min = 14
mag_max = 8


; circular mask to avoid reference sources near saturated stars or
; near edge
dummy = fltarr(2*maskrad+1,2*maskrad+1)
dummy[*,*] = 1
circmask = circ_mask(dummy,maskrad,maskrad,maskrad)

psf_fwhm = 3
delta_r = psf_fwhm * 3.
delta_mag = 5.

psf_size = 2 * (maskrad + 10) + 1
n_pix = float(psf_size)^2


; Parameters for PSF estimation and StarFinder
; ------------------------------------------------

correl_mag = 4.0
deblend = 0
deblost = 0
niter = 2
rel_thresh = 1
guide_x = ""
guide_y = ""
weigh_med = 0
norm_max = 1


; 1) First estimate of PSF
; ----------------------------

psf = psf_gaussian(NPIXEL=psf_size,FWHM=psf_fwhm,/NORMALIZE,/DOUBLE)

threshold = 10.*noise
background = estimate_background(im,2*maskrad)
search_objects, im, LOW_SURFACE = background, threshold, $
	                PRE_SMOOTH = 1, MINIF = 2, $
	                n, x, y, f
good = where(f gt 0,n)
x = x[good]
y = y[good]
f = f[good]

m = ZP - 2.5 * alog10(f/DIT) 

; -----------------------------
; Select reference stars
; For first PSF extraction I do not set any  magnitude limits
; because the fluxes returned by search_objects are mcuh lower and
; also more inaccurate than what StarFinder returns
; -----------------------------

 ; stars must be sorted by decreasing magnitude
 ; so that the brightest ones will be used to extract the PSF
 ord = sort(m)
 x_psf = x[ord]
 y_psf = y[ord]
 m_psf = m[ord]
 n_ref = n_elements(m_psf)
 print, 'Found '+ strn(n_ref) + ' stars.'


; Exclude saturated sources or sources close to saturated ones
; ----------------------------------------------------------------

 fov_mask = support
 sat_pixels = where(im gt satlevel, complement=not_saturated,n_saturated)
 if (n_saturated gt 0) then fov_mask[sat_pixels] = 0
 writefits, tmpdir + 'fov_mask_sat'+strn(i)+'_chip'+strn(chip)+'_.fits', fov_mask
 fov_mask = round(CONVOLVE(fov_mask,circmask)) ; ROUND is important, othrwise there will be lose pixels
; writefits, tmpdir + 'fov_mask_convolved.fits', fov_mask
 maskpix = where(fov_mask lt total(circmask),complement=goodpix)
 fov_mask[maskpix] = 0
 fov_mask[goodpix] = 1
 writefits, tmpdir + 'fov_mask'+strn(i)+'_chip'+strn(chip)+'_.fits', fov_mask
 
 accept = replicate(1,n_ref)
 for s = 0, n_ref-1 do begin
   xx = round(x_psf[s])
   yy = round(y_psf[s])
   xx = 0 > xx & xx = (n1 - 1) < xx
   yy = 0 > yy & yy = (n2 - 1) < yy
   if (fov_mask[xx,yy] lt 1) then accept[s] = 0
 endfor
good = where(accept gt 0, n_ref)
x_psf = x_psf[good]
y_psf = y_psf[good]
m_psf = m_psf[good] 
print, 'Found '+ strn(n_ref) + '  unsaturated reference stars far from saturated stars.'

if (n_ref lt n_ref_min) then begin
  print, 'Could not find enough reference stars with sufficient support.'
  STOP
endif

; Select isolated stars
isolated_stars, x, y, m, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
                
x_psf = x_psf[ind_iso]
y_psf = y_psf[ind_iso]
m_psf = m_psf[ind_iso]
n_ref = n_elements(m_psf)
print, 'Found '+ strn(n_ref) +  ' supported, unsaturated and isolated reference stars.'

; Use the brightest stars
; At least n_ref_min stars
  if n_ref gt n_ref_max then begin
     n_ref = n_ref_max
  endif else begin
   if n_ref gt n_ref_min then n_ref = n_ref
   if n_ref lt n_ref_min then begin
        print, 'Cannot find enough PSF reference stars.'
        STOP
   endif
  endelse
       
 x_psf = x_psf[0:n_ref-1]
 y_psf = y_psf[0:n_ref-1]
 m_psf = m_psf[0:n_ref-1]  
 f_psf = 10^(-0.4*m_psf)
 print, 'using '+ strn(n_ref) + ' PSF reference stars.'
 
 print
 print, 'Magnitudes of reference stars: '
 print, m_psf

; -----------------------------

debug = 0
nrad = 4
iter = 0                        ; more iterations do not necessarily make it better
mindist = 4
MAKEPSF, x_psf, y_psf, x, y, f, im, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_noise, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=1

 psf = centroider(psf)
 MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked
 psf = psf_masked  
 psf = psf/total(psf)  ; normalization of PSF
 writefits, tmpdir + 'tmppsf0'+strn(i)+'_chip'+strn(chip)+'_.fits', psf



  ; Run StarFinder and iterate search for PSF reference stars and PSF extraction
    ; --------------------------------------------------------------------------

  Threshold = [10.]
  estim_bg = 1
  back_box = maskrad
  starfinder, im, psf, X_BAD=x_bad, Y_BAD = y_bad, $
        BACKGROUND = background, BACK_BOX = back_box, $
        threshold, REL_THRESHOLD = rel_thresh, /PRE_SMOOTH, $
        NOISE_STD = noise, min_correlation, $
        CORREL_MAG = correl_mag, INTERP_TYPE = 'I', $
        ESTIMATE_BG = estim_bg, DEBLEND = deblend, N_ITER = niter, SILENT=1, $
        GUIDE_X = guide_x, GUIDE_Y = guide_y, $
        SV_SIGMA_R = 0.0085, SV_SIGMA_A= 0.0050, $
        x, y, f, sx, sy, sf, c, STARS = stars, $
        LOGFILE = logfilename, /CUBIC;, /NO_SLANT
    writefits, tmpdir + 'stars'+strn(i)+'_chip'+strn(chip)+'_.fits', stars
    writefits, tmpdir + 'bg'+strn(i)+'_chip'+strn(chip)+'_.fits', background

 ; 2) Repeat PSF extraction
 ; ---------------------

  m = ZP - 2.5 * alog10(f/DIT)

; -----------------------------
; Select reference stars
; -----------------------------

 good = where(m lt mag_min and m gt mag_max,n_ref)
 x_psf = x[good]
 y_psf = y[good]
 m_psf = m[good]
 ; stars must be sorted by decreasing magnitude
 ; so that the brightest ones will be used to extract the PSF
 ord = sort(m_psf)
 x_psf = x_psf[ord]
 y_psf = y_psf[ord]
 m_psf = m_psf[ord]

 ; We can only accept PSF reference stars that lie in the suport region
 ; --------------------------------------------------------------------
print, 'Found '+ strn(n_ref) + ' stars in proposed magnitude range.'

; Now exclude saturated sources or sources close to saturated ones
; ----------------------------------------------------------------
 accept = replicate(1,n_ref)
 for s = 0, n_ref-1 do begin
   xx = round(x_psf[s])
   yy = round(y_psf[s])
   xx = 0 > xx & xx = (n1 - 1) < xx
   yy = 0 > yy & yy = (n2 - 1) < yy
   if (fov_mask[xx,yy] lt 1) then accept[s] = 0
endfor
good = where(accept gt 0, n_ref)
x_psf = x_psf[good]
y_psf = y_psf[good]
m_psf = m_psf[good] 
print, 'Found '+ strn(n_ref) + ' unsaturated reference stars far from saturated ones.'

if (n_ref lt n_ref_min) then begin
  print, 'Could not find enough reference stars with sufficient support.'
  STOP
endif

isolated_stars, x, y, m, x_psf, y_psf, m_psf, delta_mag, delta_r, ind_iso
x_psf = x_psf[ind_iso]
y_psf = y_psf[ind_iso]
m_psf = m_psf[ind_iso]
n_ref = n_elements(m_psf)
print, 'Found  '+ strn(n_ref) + ' supported, unsaturated and isolated reference stars.'

; Use the brightest stars
; At least n_ref_min stars
  if n_ref gt n_ref_max then begin
     n_ref = n_ref_max
  endif else begin
   if n_ref gt n_ref_min then n_ref = n_ref
   if n_ref lt n_ref_min then begin
        print, 'Cannot find enough PSF reference stars.'
        STOP
   endif
  endelse
       
 x_psf = x_psf[0:n_ref-1]
 y_psf = y_psf[0:n_ref-1]
 m_psf = m_psf[0:n_ref-1]  
 f_psf = 10^(-0.4*m_psf)
 print, 'Found '+ strn(n_ref) + ' PSF reference stars.'
 dat = ptr_new({X_size: 40, Y_size: 40, Sigma_x: 2., Sigma_y: 2., Angle: 0.0})
 refstars = image_model(x_psf,y_psf,f_psf,n1,n2,'gaussian', dat)
 writefits, tmpdir + 'refstars'+strn(i)+'_chip'+strn(chip)+'_.fits', refstars
 forprint, TEXTOUT= outdir+'psfstars'+strn(i)+'_chip'+strn(chip)+'_.txt', x_psf, y_psf, m_psf, /NOCOMMENT

 print
 print, 'Magnitudes of reference stars: '
 print, m_psf

; -----------------------------

  debug = 1
  nrad = 4
  mindist = 4
  iter = 1  ; more iterations do not necessarily make it better

  MAKEPSF, x_psf, y_psf, x, y, f, im, nrad, FOVMASK = fov_mask, PSF=psf,  BACKGROUND=background, DEBUG = debug, ITER = iter, MINDIST = mindist, NOISE_PSF = psf_noise, MASKRAD = maskrad, UNWEIGHTED=unweighted, TMPDIR = tmpdir, LOCAL_SKY=1

  writefits, tmpdir + 'tmppsf1'+strn(i)+'_chip'+strn(chip)+'_.fits', psf

  psf = centroider(psf)
  writefits, tmpdir + 'tmppsf'+strn(i)+'_chip'+strn(chip)+'_nomask.fits', psf
  MASK_PSF, psf, maskrad, PSF_MASKED=psf_masked
  psf = psf_masked  
  writefits, tmpdir + 'tmppsf'+strn(i)+'_chip'+strn(chip)+'_masked.fits', psf  
  psf = psf/total(psf)  ; normalization of PSF
  ;writefits, 'psf.fits', psf
  writefits, outdir +'psf_im'+strn(i)+'_chip'+strn(chip)+'_.fits', psf
print, '#################### DONE with image',i,'chip',chip,'DIT',DIT,' ####################'
print, '#################### DONE with image',i,'chip',chip,'DIT',DIT,' ####################'
print, '#################### DONE with image',i,'chip',chip,'DIT',DIT,' ####################'
endfor

endfor
END
