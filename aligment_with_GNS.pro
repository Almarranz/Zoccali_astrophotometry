pro aligment_with_GNS,field,lst
;lst  can be 1 or 4 (refers to the chip on GNS fields)
;~ field=20 ; fields can be 3 or 20 (refers to GNS fields)
;~ Esto es una prueba para git pull
;~ NOTE:
;~ Lists 1 to 3 are on the brick , Zoc chip 3
;~ Lists 0 are on the brick , Zoc chip 2
;~ Lists 10 is on chip 2 out of brick
;~ List 16 and 12 are on chip 3 out of brick


band='H'

 
exptime=10
folder='im_jitter_NOgains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field'+strn(field)+'/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp_bs/'
gaussian='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'Gaussian_fit/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results/'
tmp_p=pruebas
name='NPL058_'

markstars=0
rot_angle=0


readcol, GNS + 'cat_Ban_'+strn(field)+'_'+strn(lst)+'.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A';,SKIPLINE = 1

if field eq 3 then begin
	readcol,tmp+'stars_calibrated_H_chip2_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
endif

if field eq 20 then begin
	readcol,tmp+'stars_calibrated_H_chip3_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
endif

if field eq 16 then begin
	readcol,tmp+'stars_calibrated_H_chip2_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
endif
if field eq 10 then begin
	readcol,tmp+'stars_calibrated_H_chip1_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
endif
a=float(a)
d=float(d)
f=float(f)
df=float(df)
m=float(m)
dm=float(dm)
x=float(x)
y=float(y)
dx=float(dx)
dy=float(dy)





if markstars eq 0 then begin
    
    if field eq 3 then begin
    
		if lst eq 1 then begin
		       

	
			xm_ref= 1703.93  ; xm_ref is GNS
			ym_ref= 1161.27  

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1000.1831055
			ym=330.6050415
    
		endif
		
		if lst eq 4 then begin
		       
 

		
		    xm_ref= 2639.51   ; xm_ref is GNS
			ym_ref= 1032.54

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1478.5396729
			ym=1109.546875
		
		endif 
    endif
    
    if field eq 20 then begin
    
		if lst eq 1 then begin
		       

	       
 

			xm_ref= 843.598   ; xm_ref is GNS
			ym_ref= 805.357  

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=595.5708008
			ym=751.9746094
    
		endif
		
		if lst eq 4 then begin
		
		
		    xm_ref= 1785.32    ; xm_ref is GNS
			ym_ref= 986.934 

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1080.1188965
			ym=1684.3498535
		
		endif 
	endif
		
	if field eq 16 then begin
    
		if lst eq 3 then begin

			xm_ref=  1607.22   ; xm_ref is GNS mg: 18.860020  16.975086
			ym_ref= 1269.01 

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1496.7492676
			ym=1823.2287598
    
		endif
		
		if lst eq 2 then begin
		
		        
 

		    xm_ref= 1432.55    ; xm_ref is GNS 16.829983  15.191847
			ym_ref= 1541.23

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1397.8381348
			ym=1110.1224365
		
		endif 
    endif
    
    if field eq 10 then begin
    
		if lst eq 3 then begin
		        
 


			xm_ref=  731.097  ; xm_ref is GNS mg: 15.879695  13.938910
			ym_ref= 1875.66

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1022.0848389
			ym=1543.3221436
    
		endif
		
		if lst eq 2 then begin
		        
 

		    xm_ref=  879.621     ; xm_ref is GNS 117.055305           15.188314
			ym_ref= 1578.72 

			xm_ref=xm_ref*0.5
			ym_ref=ym_ref*0.5
			

			xm=1079.5262451
			ym=544.4326172
		
		endif 
    endif
    
    
	print, '#######################'
	print, 'xoff, yoo ',xm_ref-xm,ym_ref -ym
	print, '#######################'
endif



raH=float(raH)
decH=float(decH)

mH=float(mH)
mK=float(mK)

;coorecte the rebinning on GNS

x_gns=x_gns*0.5
y_gns=y_gns*0.5





valid_H=where(mH lt 90 and mK lt 90)
raH=raH[valid_H]
decH=decH[valid_H]
mH=mH[valid_H]
mK=mK[valid_H]
x_gns=x_gns[valid_H]
y_gns=y_gns[valid_H]
print,n_elements(raH)

H_Ks=where(mH-mK>1.3)

raH=raH[H_Ks]
decH=decH[H_Ks]
mH=mH[H_Ks]
mK=mK[H_Ks]
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]



if (markstars gt 0) then begin
	imagen=readfits(GNS+'field'+strn(field)+'_chip'+strn(lst)+'.fits',header)
	sz = size(imagen)
	xsize_ref = sz[1]*0.5
	ysize_ref= sz[2]*0.5
	EXTAST, header, astr 



	AD2XY,raH, decH, astr, x_gns, y_gns
	AD2XY,a,d,astr,x_brick,y_brick

	x=x_brick
	y=y_brick
endif

; mark GNS stars first
	 if (markstars gt 0) then begin
	  RETURNMARKED, xsize_ref, ysize_ref, x_gns, y_gns, 10^((-mH)/2.5)*10, XM = xm_ref, YM = ym_ref, FM = fm_ref, BOXSIZE = 31, dmax = 1.0, g_sigma = 2.0, DISP_STRETCH = 'linear', DISP_LARGE=2
	  SAVE, xm_ref, ym_ref, FILENAME= tmp + 'Refstars_GNS_' + band
	  
	  RESTORE, tmp + 'Refstars_GNS_' + band
	 endif
	 
	 

; now mark HAWK-I stars
	 if (markstars gt 0) then begin
	  RETURNMARKED, xsize_ref, ysize_ref, x_brick, y_brick, f, XM = xm, YM = ym, FM = fm, BOXSIZE = 21, dmax = 1., g_sigma = 3.0, DISP_STRETCH = 'linear', DISP_LARGE=2
	  SAVE, xm, ym, FILENAME= tmp + 'Refstars_BRICK_' + band
	  
	  RESTORE, tmp + 'Refstars_BRICK_' + band
	 endif
	 



; preliminary offset and rotation
	 ; -------------------------------

	 xm = xm * cos(rot_angle) - ym * sin(rot_angle)
	 ym = xm * sin(rot_angle) + ym * cos(rot_angle)
	 if (n_elements(xm) gt 1) then xoff = median(xm_ref - xm) else xoff = xm_ref - xm
	 if (n_elements(ym) gt 1) then yoff = median(ym_ref - ym) else yoff = ym_ref - ym
	 xi = x * cos(rot_angle) - y * sin(rot_angle)
	 yi = x * sin(rot_angle) + y * cos(rot_angle)
	 ;~ xoff=-640.8315405000001 
	 ;~ yoff=217.90724489999997 
	 xi = xi + xoff
	 yi = yi + yoff
	 print,'This is x_off and y_off',xoff,yoff
		


     
	 x0 = xi
	 y0 = yi
	 

	; dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
	; map = image_model(xi,yi,f,xsize_quad,ysize_quad,'gaussian', dat)
	; writefits, tmp_path + 'align_sources.fits', map

	 dmax = 1
	 compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	 nc = n_elements(subc1)
	 print, 'Found ' + strn(nc) + ' common stars.'
	 
	 
	 
	 
     ;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c-xoff ,dx[subc2] , y2c-yoff, dy[subc2], x1c,dx_gns[subc1]/0.106,y1c,dy_gns[subc1]/0.106 ,format='(10(f, 4X))', /NOCOMMENT 
     ;~ stop
     ; iterative degree 1 alignment
	 ; ------------------------------
     count=0
     comm=[]
     it=0
     lim_it=2
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

     ; iterative degree 2 alignment
 ; ------------------------------
     print, '#######################'
	 print, 'Now Degree 2 alignment.'
	 print, '#######################'
	 count=0
     comm=[]
     it=0
	 while count lt lim_it do begin
	  it=it+1
	  degree = 2
	  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], degree, Kx, Ky
	  print, Kx
	  print, Ky
	  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
	  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x
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
	

	;~ readcol, GNS+'field12_on_brick.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
	;~ raH=float(raH)
	;~ decH=float(decH)

	;~ mH=float(mH)
	;~ mK=float(mK)

	;~ ;coorecte the rebinning on GNS

	;~ x_gns=x_gns*0.5
	;~ y_gns=y_gns*0.5
    
    distancia=1
    compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=distancia, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	nc = n_elements(subc1)
    print, 'Found after '+strn(degree)+' degree alignment ' + strn(nc) + ' common stars.'
    
    print, N_ELEMENTS(f)
    f=f[subc2]
    print,N_ELEMENTS(f)
    
    if (markstars gt 0) then begin
		dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
		im = image_model(x_gns,y_gns,10^((-mH)/2.5),xsize_ref,xsize_ref,'gaussian', dat)
		ima = image_model(xi,yi,f,xsize_ref,xsize_ref,'gaussian', dat)
	
	
	writefits, tmp_p +'gns'+strn(field)+'_chip' + strn(lst) + '.fits',im

	
	writefits, tmp_p +'zoc_on_GNS_field'+ strn(field) + '.fits', ima
	endif
	
    ;~ forprint, TEXTOUT= tmp_p+'2lists_IDL_GNS.txt',x1c ,x2c , y1c, y2c,format='(10(f, 4X))', /NOCOMMENT 
    if (field eq 3) or (field eq 10) or (field eq 16) then begin ; different dates for different fields
		x_dis=(x2c-x1c)*0.106/4.3*1000
		y_dis=(y2c-y1c)*0.106/4.3*1000
		;~ Adding velocities uncertanties for x and y directions
		dx_gns=dx_gns[subc1]*0.106*1000
		dy_gns=dy_gns[subc1]*0.106*1000
		
		dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		dy=dy[subc2]*1000
		
		dvx=sqrt(dx^2+dx_gns^2)/4.3
		dvy=sqrt(dy^2+dy_gns^2)/4.3
	endif
	
	if field eq 20 then begin
		x_dis=(x2c-x1c)*0.106/4.2*1000
		y_dis=(y2c-y1c)*0.106/4.2*1000
		;~ Adding velocities uncertanties for x and y directions
		dx_gns=dx_gns[subc1]*0.106*1000
		dy_gns=dy_gns[subc1]*0.106*1000
		
		dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		dy=dy[subc2]*1000
		
		dvx=sqrt(dx^2+dx_gns^2)/4.2
		dvy=sqrt(dy^2+dy_gns^2)/4.2
		
	endif
    
    mH=mH[subc1]
    mK=mK[subc1]
    
    forprint, TEXTOUT= tmp+name+'IDL_xdis_ydis_field'+strn(field)+'_chip'+strn(lst)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,format='(10(f, 4X))', /NOCOMMENT 
   
	;~ forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,format='(10(f, 4X))', /NOCOMMENT 
	forprint, TEXTOUT= gaussian+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,format='(10(f, 4X))', /NOCOMMENT 
		
    
    
    

	;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c ,dx , y2c, dy, x1c,dx_gns/0.106,y1c,dy_gns/0.106 ,format='(10(f, 4X))', /NOCOMMENT 

end

