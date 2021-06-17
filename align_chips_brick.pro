pro align_chips_brick,chip


;chip=2
band='H'
exptime=10
folder='im_jitter_NOgains/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results/'
tmp_p=pruebas
name='NPL_054'

magerr_si = 0.02 
;rot_angle = 0.4 * !PI/180. ; manually estimated angle of rotation between HAWK-I observations and VVV field
rot_angle = 0
; mark stars Y/N
markstars = 1

;----------------------- Pixel coordinates from SIRUS on the brick -----------------------
; Read list of reference stars(the list from SIRIUS the brick)
readcol, sirius +'the_brick_idl.txt',a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si, Format = 'A,A,A,A,A,A,A,A'
a_si=float(a_si)
dec_si=float(dec_si)
J_si=float(J_si)
dJ_si=float(dJ_si)
H_si=float(H_si)
dH_si=float(dH_si)
K_si=float(K_si)
dK_si=float(dK_si)   

readcol, tmp + 'NPL_054fluxes_chip'+strn(chip)+'.txt', a,d,f,df,x,y, Format ='A,A,A,A,A,A'
a=float(a)
d=float(d)
f=float(f)
df=float(df)
x=float(x)
y=float(y)
 

imagen=readfits(indir+'cube_chip'+strn(chip)+'_canvas.fits',EXT=1,header)
dat=readfits(indir+'cube_chip'+strn(chip)+'_canvas.fits',EXT=0,cabeza)

sz = size(imagen)
xsize_ref = sz[1]
ysize_ref= sz[2]
EXTAST, header, astr 

if band eq 'H' then begin
 m_si = H_si
 dm_si = dH_si
endif

if band eq 'Ks' then begin
 m_si = K_si
 dm_si = dK_si
endif

; use only valid measurements
valid = where(m_si gt 0)
a_si = a_si[valid]
dec_si = dec_si[valid]
m_si = m_si[valid]
dm_si = dm_si[valid]

; Select calibration stars by brightness
if band eq 'J' then begin
  good = where(m_si gt 11.0 and dm_si lt magerr_si) 
endif

if band eq 'H' then begin
;  good = where(m_si gt 11.0  and m_si lt 14 and dm_si lt 0.04)
  good = where(m_si gt 11.0  and dm_si lt magerr_si)
endif

if band eq 'Ks' then begin
  good = where(m_si gt 11.0 and dm_si lt magerr_si)    
endif

; appendix _si means all valid stars in SIRIUS catalogue
; appendic _ref means potential photmetric reference stars
a_ref = a_si[good]
dec_ref = dec_si[good]
m_ref = m_si[good]
dm_ref = dm_si[good]  

AD2XY, a_ref, dec_ref, astr, x_ref, y_ref
AD2XY, a_si, dec_si, astr, x_all, y_all



; mark VVV stars first
 if (markstars gt 0) then begin
  RETURNMARKED, xsize_ref, ysize_ref, x_ref, y_ref, 10^((-m_ref)/2.5), XM = xm_ref, YM = ym_ref, FM = fm_ref, BOXSIZE = 31, dmax = 2.0, g_sigma = 2.0, DISP_STRETCH = 'linear', DISP_LARGE=1
  SAVE, xm_ref, ym_ref, FILENAME= tmp + 'Refstars_SIR_' + strn(chip)
 endif
 
 RESTORE, tmp + 'Refstars_SIR_' + strn(chip)
 ;stop
 
 ; now mark HAWK-I stars
 if (markstars gt 0) then begin
  RETURNMARKED, xsize_ref, ysize_ref, x, y, f, XM = xm, YM = ym, FM = fm, BOXSIZE = 21, dmax = 2., g_sigma = 3.0, DISP_STRETCH = 'linear', DISP_LARGE=1
  SAVE, xm, ym, FILENAME= tmp + 'Refstars_HAWKI_' + strn(chip)
 endif
 RESTORE, tmp + 'Refstars_HAWKI_' + strn(chip)
 
 ;stop
 ; preliminary offset and rotation
 ; -------------------------------

 xm = xm * cos(rot_angle) - ym * sin(rot_angle)
 ym = xm * sin(rot_angle) + ym * cos(rot_angle)
 if (n_elements(xm) gt 1) then xoff = median(xm_ref - xm) else xoff = xm_ref - xm
 if (n_elements(ym) gt 1) then yoff = median(ym_ref - ym) else yoff = ym_ref - ym
 xi = x * cos(rot_angle) - y * sin(rot_angle)
 yi = x * sin(rot_angle) + y * cos(rot_angle)
 xi = xi + xoff
 yi = yi + yoff
 print,'This is x_off and y_off',xoff,yoff

 x0 = xi
 y0 = yi
 

; dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
; map = image_model(xi,yi,f,xsize_quad,ysize_quad,'gaussian', dat)
; writefits, tmp_path + 'align_sources.fits', map

 dmax = 1.0
 compare_lists, x_ref, y_ref, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
 nc = n_elements(subc1)
 print, 'Found ' + strn(nc) + ' common stars.'
 
 
 ; iterative degree 1 alignment
 ; ------------------------------

 for it = 1, 20 do begin
  degree = 1
  polywarp, x_ref[subc1], y_ref[subc1], x[subc2], y[subc2], degree, Kx, Ky
  print, Kx
  print, Ky
  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
  compare_lists, x_ref, y_ref, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
  nc = n_elements(subc1)
  print, 'Iteration ' + strn(it)
  print, 'Found ' + strn(nc) + ' common stars.'
endfor

dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
im = image_model(x_ref,y_ref,10^((-m_ref)/2.5),xsize_ref,xsize_ref,'gaussian', dat)
ima = image_model(xi,yi,f,xsize_ref,xsize_ref,'gaussian', dat)


xy2ad, xi,yi,astr,ra_alig,dec_alig



;forprint, TEXTOUT= tmp+'alig_SIRUS_chip'+strn(chip)+'.txt',ra_alig,dec_alig,f,df,xi,yi, /NOCOMMENT 
forprint,format='(7f13.5)', TEXTOUT= tmp+'alig_SIRUS_chip'+strn(chip)+'.txt',ra_alig,dec_alig,f,df,xi,yi, /NOCOMMENT 

writefits, tmp_p +'sirius_chip' + strn(chip)+ '_0.fits',dat, cabeza
writefits, tmp_p +'sirius_chip' + strn(chip) + '_0.fits',im, header,/app

writefits, tmp_p +'mean1_chip' + strn(chip) + '_0.fits', dat, cabeza
writefits, tmp_p +'mean1_chip' + strn(chip) + '_0.fits', ima,header,/app





;STOP
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 end
