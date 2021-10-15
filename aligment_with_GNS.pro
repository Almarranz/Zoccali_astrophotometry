pro aligment_with_GNS,lst,degree3

;~ NOTE:
;~ Lists 1 to 3 are on the brick , Zoc chip 3
;~ Lists 0 are on the brick , Zoc chip 2
;~ Lists 10 is on chip 2 out of brick
;~ List 16 and 12 are on chip 3 out of brick

;~ On this branch (astro_aling) we dont need to click on stars, we get the matriz transformation with aa in python

if (lst eq 10) || (lst eq 0) then chip=2 else chip=3
band='H'
exptime=10
folder='im_jitter_NOgains/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
GNS_ori='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field_out/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp_bs/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results/'
tmp_p=pruebas
name='NPL_054'
markstars=0
rot_angle=0

;####################
;degree3=1 ; if degree3 = 1 uses a deg3 polinomial in the alignment, if =0 only a second degree pol. is used
;####################
;~ lst=1
;~ if lst eq 16 then readcol, GNS_ori+'field16_out_of_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
if lst gt 4 then readcol, GNS_ori+'field'+strn(lst)+'_out_of_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
;~ if lst gt 4 then GNS_ori+'field'+strn(lst)+'_out_of_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
if (lst eq 10) || (lst eq 0) then readcol, GNS+'field12_on_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
if lst eq 1 then readcol, GNS+'field12_on_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
if lst eq 2 then readcol, GNS+'field12_on_brick_accu.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
if lst eq 3 then readcol, GNS+'field12_on_brick_reduced.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
print, '#######################'
print, 'Reading lst = ',lst
print, '#######################'




if lst gt 4 then begin
	readcol, tmp+'OUT'+strn(lst)+'_stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
;~ readcol, tmp+'stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
endif else begin
    readcol, tmp+'aa_BRICK_stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
    print,'.....',n_elements(x),'.....'
    ;~ readcol, tmp+'BRICK_stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
    ;~ print,'.....',n_elements(x),'.....'
endelse
a=float(a)
d=float(d)
f=float(f)
df=float(df)
m=float(m)
dm=float(dm)

;~ xi=float(xi)
;~ yi=float(yi)
x=float(x)
y=float(y)

dx=float(dx)
dy=float(dy)




raH=float(raH)
decH=float(decH)

mH=float(mH)
mK=float(mK)

;coorecte the rebinning on GNS

x_gns=x_gns*0.5
y_gns=y_gns*0.5

print,'Elemnts in GNS',n_elements(x_gns)


valid_H=where(mH lt 90 and mK lt 90)
raH=raH[valid_H]
decH=decH[valid_H]
mH=mH[valid_H]
mK=mK[valid_H]
x_gns=x_gns[valid_H]
y_gns=y_gns[valid_H]

print,'Elemnts in GNS',n_elements(x_gns)

H_Ks=where((mH-mK) gt 1.3)

raH=raH[H_Ks]
decH=decH[H_Ks]
mH=mH[H_Ks]
mK=mK[H_Ks]
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]

print,'Elemnts in GNS',n_elements(x_gns)


imagen=readfits(GNS+'field12_'+band+'.fits',header)
sz = size(imagen)
xsize_ref = sz[1]*0.5
ysize_ref= sz[2]*0.5
EXTAST, header, astr 




	 

	; dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
	; map = image_model(xi,yi,f,xsize_quad,ysize_quad,'gaussian', dat)
	; writefits, tmp_path + 'align_sources.fits', map

	 dmax = 1
	 compare_lists, x_gns, y_gns, x, y, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	 nc = n_elements(subc1)
	 print, 'Found ' + strn(nc) + ' common stars.'
     ;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c-xoff ,dx[subc2] , y2c-yoff, dy[subc2], x1c,dx_gns[subc1]/0.106,y1c,dy_gns[subc1]/0.106 ,format='(10(f, 4X))', /NOCOMMENT 
     ;~ stop
     ; iterative degree 1 alignment
	 ; ------------------------------
     count=0
     comm=[]
     it=0
     lim_it=1
	 while count lt lim_it do begin
	  it=it+1
	  degree = 1
	  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], degree, Kx, Ky
	  print, Kx
	  print, Ky
	  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
	  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
	  compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	  nc = n_elements(subc1)
	  print, 'Iteration ' + strn(it)
	  print, 'Found ' + strn(nc) + ' common stars.'
	  comm=[comm,nc]
	  if (n_elements(comm) gt 2) then begin
	   if comm[-2] eq comm[-1] then begin
	   count=count+1
	  endif else begin
	   count=0
	  endelse
	  endif
	 endwhile

     ;~ ; iterative degree 2 alignment
 ;~ ; ------------------------------
     ;~ print, '#######################'
	 ;~ print, 'Now Degree 2 alignment.'
	 ;~ print, '#######################'
	 ;~ count=0
     ;~ comm=[]
     ;~ it=0
	 ;~ while count lt lim_it do begin
	  ;~ it=it+1
	  ;~ degree = 2
	  ;~ polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], degree, Kx, Ky
	  ;~ print, Kx
	  ;~ print, Ky
	  ;~ xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
	  ;~ yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x
	  ;~ compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	  ;~ nc = n_elements(subc1)
	  ;~ print, 'Iteration ' + strn(it)
	  ;~ print, 'Found ' + strn(nc) + ' common stars.'
	  ;~ comm=[comm,nc]
	  ;~ if (n_elements(comm) gt 2) then begin
	   ;~ if comm[-2] eq comm[-1] then begin
	   ;~ count=count+1
	  ;~ endif else begin
	   ;~ count=0
	  ;~ endelse
	  ;~ endif
	;~ endwhile
	
	if degree3 eq 1 then begin
		 ; iterative degree 3 alignment
	 ; ------------------------------
		 print, '#######################'
		 print, 'Now Degree 3 alignment.'
		 print, '#######################'
		 count=0
		 comm=[]
		 it=0
		 while count lt lim_it && it lt 101 do begin
		  it=it+1
		  degree = 3
		  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], degree, Kx, Ky
		  print, Kx
		  print, Ky
		  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x +$
			Kx[0,3]*x^3 + Kx[1,3]*x^3*y + Kx[2,3]*x^3*y^2 + Kx[3,0]*y^3 + Kx[3,1]*x*y^3 + Kx[3,2]*x^2*y^3 + Kx[3,3]*x^3*y^3
		  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x +$
		   Ky[0,3]*x^3 + Ky[1,3]*x^3*y + Ky[2,3]*x^3*y^2 + Ky[3,0]*y^3 + Ky[3,1]*x*y^3 + Ky[3,2]*x^2*y^3 + Ky[3,3]*x^3*y^3

		  compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
		  nc = n_elements(subc1)
		  print, 'Iteration ' + strn(it)
		  print, 'Found ' + strn(nc) + ' common stars.'
		  comm=[comm,nc]
		  if (n_elements(comm) gt 2) then begin
		   if comm[-2] eq comm[-1] then begin
		   count=count+1
		  endif else begin
		   count=0
		  endelse
		  endif
		endwhile
	endif
	
	
    
    distancia=1
    compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=distancia, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	nc = n_elements(subc1)
    print, 'Found after '+strn(degree)+' degree alignment ' + strn(nc) + ' common stars.'
    
    print, N_ELEMENTS(f)
    f=f[subc2]
    print,N_ELEMENTS(f)
    
    dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
	im = image_model(x_gns,y_gns,10^((-mH)/2.5),xsize_ref,xsize_ref,'gaussian', dat)
	ima = image_model(xi,yi,f,xsize_ref,xsize_ref,'gaussian', dat)
	
	
	writefits, tmp_p +'gns12_chip' + strn(chip) + '.fits',im

	
	writefits, tmp_p +'brick_chip' + strn(chip) + '.fits', ima
	
    forprint, TEXTOUT= tmp_p+'2lists_IDL_GNS.txt',x1c ,x2c , y1c, y2c,format='(10(f, 4X))', /NOCOMMENT 
    
    x_dis=(x2c-x1c)*0.106/4.1*1000
    y_dis=(y2c-y1c)*0.106/4.1*1000
    ;~ Adding velocities uncertanties for x and y directions
    dx_gns=dx_gns[subc1]*0.106*1000
    dy_gns=dy_gns[subc1]*0.106*1000
    
    dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
    dy=dy[subc2]*1000
    
    dvx=sqrt(dx^2+dx_gns^2)/4.1
    dvy=sqrt(dy^2+dy_gns^2)/4.1
    
    mH=mH[subc1]
    mK=mK[subc1]
    m=m[subc2]
    
    forprint, TEXTOUT= tmp+'aa_IDL_xdis_ydis_chip'+strn(chip)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,format='(10(f, 4X))', /NOCOMMENT 
    if lst gt 4 then begin
		forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+'aa_IDL_arcsec_vx_vy_chip'+strn(chip)+'_out_Brick'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,format='(10(f, 4X))', /NOCOMMENT 
		;~ forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+'IDL_arcsec_vx_vy_chip3.txt',x_dis,y_dis,mH,mK,format='(10(f, 4X))', /NOCOMMENT 
    endif else begin
		forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+'aa_IDL_arcsec_vx_vy_chip'+strn(chip)+'.txt',x_dis,y_dis,dvx,dvy,mH,m,format='(10(f, 4X))', /NOCOMMENT
    endelse
    forprint, TEXTOUT= tmp +'aa_IDL_lst_chip'+strn(chip)+'.txt',lst, format='I', /NOCOMMENT 
    
    

	;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c ,dx , y2c, dy, x1c,dx_gns/0.106,y1c,dy_gns/0.106 ,format='(10(f, 4X))', /NOCOMMENT 

end

