pro aligment_with_GNS,field,lst,degree,color
;lst  can be 1 to 4 (refers to the chip on GNS fields)
;~ field=20 ; fields can be 3 or 20 or 16 (refers to GNS fields)
;~ Esto es una prueba para git pull
;~ NOTE:
;~ degree: degree of the poly. fit. (1 and 2=2. 1,2, and 3=3)

;~ IMPORTAT!!!!!!!!!!!!!!
;~ Check if Z1 is in front of the lists names


band='H'

exptime=10
folder='im_jitter_NOgains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field'+strn(field)+'/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp_bs/'
gaussian='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'Gaussian_fit/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results/'
tmp_p=pruebas
;~ name='aa_NPL058_'
;~ name='Z1_aa_NPL058'
zone='Z1' ;Zone refers to the little squares the same size of the one of the brick on zone B; Z1 is left, z2 is right
name=strn(zone)+'_aa_NPL058'
print, name

markstars=0

;~ readcol, GNS + 'cat_Ban_'+strn(field)+'_'+strn(lst)	+'.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A';,SKIPLINE = 1
readcol, GNS + strn(zone)+'_cat_Ban_'+strn(field)+'_'+strn(lst)+'.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1

;~ readcol,tmp+'stars_calibrated_H_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
;~ readcol,tmp+'aa_stars_calibrated_H_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1
readcol,tmp+strn(zone)+'_aa_stars_calibrated_H_on_field'+strn(field)+'_'+strn(lst)+'.txt',a ,d , m, dm, f, df,x,y,dx,dy,Format ='A,A,A,A,A,A,A,A,A,A',SKIPLINE = 1



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
dx_gns=float(dx_gns)
dy_gns=float(dy_gns)  	









raH=float(raH)
decH=float(decH)

mH=float(mH)
mK=float(mK)

;coorecte the rebinning on GNS

x_gns=x_gns*0.5
y_gns=y_gns*0.5


;~ dx_gns=dx_gns*0.5
;~ dy_gns=dy_gns*0.5





valid_H=where(mH lt 90 and mK lt 90)
raH=raH[valid_H]
decH=decH[valid_H]
mH=mH[valid_H]
mK=mK[valid_H]
x_gns=x_gns[valid_H]
y_gns=y_gns[valid_H]

;~ dx_gns=dx_gns[valid_H]
;~ dy_gns=dy_gns[valid_H]


print,n_elements(raH)

H_Ks=where((mH-mK) gt color)

raH=raH[H_Ks]
decH=decH[H_Ks]
mH=mH[H_Ks]
mK=mK[H_Ks]
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]

;~ dx_gns=dx_gns[H_Ks]
;~ dy_gns=dy_gns[H_Ks]




	 

	; dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
	; map = image_model(xi,yi,f,xsize_quad,ysize_quad,'gaussian', dat)
	; writefits, tmp_path + 'align_sources.fits', map

	 dmax = 1
	 compare_lists, x_gns, y_gns, x, y, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	 nc = n_elements(subc1)
	 print, 'Found ' + strn(nc) + ' common stars.'
	 
	 if (degree eq 0) then begin
	 
		x_dis=(x2c-x1c)*0.106/4.3*1000
		y_dis=(y2c-y1c)*0.106/4.3*1000
		;~ Adding velocities uncertanties for x and y directions
		dx_gns=dx_gns[subc1]*0.106*1000
		dy_gns=dy_gns[subc1]*0.106*1000
		
		dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		dy=dy[subc2]*1000
		
		dvx=sqrt(dx^2+dx_gns^2)/4.3
		dvy=sqrt(dy^2+dy_gns^2)/4.3
		
		distancia=1
		compare_lists, x_gns, y_gns, x, y, x1c, y1c, x2c, y2c, MAX_DISTANCE=distancia, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
		nc = n_elements(subc1)
		print, 'Found after '+strn(degree)+' degree alignment ' + strn(nc) + ' common stars.'
		
		mH=mH[subc1]
		mK=mK[subc1]
		m=m[subc2]
		a=a[subc2]
		d=d[subc2]
		raH=raH[subc1]
		decH=decH[subc1]
		
		dx_gns=dx_gns[subc1]
		dy_gns=dy_gns[subc1]
		
		dx=dx[subc2]
		dy=dy_[subc2]
		
		
		
    
		;~ forprint, TEXTOUT= tmp+name+'IDL_xdis_ydis_field'+strn(field)+'_chip'+strn(lst)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
		forprint, TEXTOUT= tmp+name+'_IDL_xdis_ydis_field'+strn(field)+'_chip'+strn(lst)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
	   
		;~ forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,format='(10(f, 4X))', /NOCOMMENT 
		;~ forprint, TEXTOUT= gaussian+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,m,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
		forprint, TEXTOUT= gaussian+name+'_IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,m,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
	
	stop
    
	 
	 endif 
	 
	 
     ;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c-xoff ,dx[subc2] , y2c-yoff, dy[subc2], x1c,dx_gns[subc1]/0.106,y1c,dy_gns[subc1]/0.106 ,format='(10(f, 4X))', /NOCOMMENT 
     ;~ stop	
     ; iterative degree 1 alignment
	 ; ------------------------------
	if (degree eq 1) || (degree eq 2) || (degree eq 3)  then begin
     print, '#######################'
	 print, 'Now Degree 1 alignment.'
	 print, '#######################' 
     count=0
     comm=[]
     it=0
     lim_it=1
	 while count lt lim_it do begin
	  it=it+1
	 
	  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], 1, Kx, Ky
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
   endif
     ; iterative degree 2 alignment
 ; ------------------------------
  if (degree eq 2) || (degree eq 3) then begin
     print, '#######################'
	 print, 'Now Degree 2 alignment.'
	 print, '#######################'
	 count=0
     comm=[]
     it=0
	 while count lt lim_it do begin
	  it=it+1
	  
	  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], 2, Kx, Ky
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
 endif
 
 if degree eq 3 then begin
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
	  
	  polywarp, x_gns[subc1], y_gns[subc1], x[subc2], y[subc2], 3, Kx, Ky
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
    ;~ if (field eq 3) or (field eq 10) or (field eq 16) or (field eq 12) or (field eq 7) then begin ; different dates for different fields
		x_dis=(x2c-x1c)*0.106/4.3*1000
		y_dis=(y2c-y1c)*0.106/4.3*1000
		;~ Adding velocities uncertanties for x and y directions
		dx_gns=dx_gns[subc1]*0.106*1000
		dy_gns=dy_gns[subc1]*0.106*1000
		
		dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		dy=dy[subc2]*1000
		
		dvx=sqrt(dx^2+dx_gns^2)/4.3
		dvy=sqrt(dy^2+dy_gns^2)/4.3
		
		
        
		
	;~ endif
	
	;~ if field eq 20 then begin
		;~ x_dis=(x2c-x1c)*0.106/4.2*1000
		;~ y_dis=(y2c-y1c)*0.106/4.2*1000
		;~ ;Adding velocities uncertanties for x and y directions
		;~ dx_gns=dx_gns[subc1]*0.106*1000
		;~ dy_gns=dy_gns[subc1]*0.106*1000
		
		;~ dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		;~ dy=dy[subc2]*1000
		
		;~ dvx=sqrt(dx^2+dx_gns^2)/4.2
		;~ dvy=sqrt(dy^2+dy_gns^2)/4.2
		
	;~ endif
    
    mH=mH[subc1]
    mK=mK[subc1]
    m=m[subc2]
    a=a[subc2]
    d=d[subc2]
    raH=raH[subc1]
    decH=decH[subc1]
    
    forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/dxy_GNS_vs_ZOC/'+'out_comm_GNS_ZOC.txt',mH,dx_gns,dy_gns,m,dx,dy,format='(6(f, 8X))', /NOCOMMENT 
    ;~ forprint, TEXTOUT= pruebas +'dvx_mag_OUT1.txt',mH,dvx,format='(2(f, 8X))', /NOCOMMENT 
    
    ;~ stop
    
		;~ forprint, TEXTOUT= tmp+name+'IDL_xdis_ydis_field'+strn(field)+'_chip'+strn(lst)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
		forprint, TEXTOUT= tmp+name+'_IDL_xdis_ydis_field'+strn(field)+'_chip'+strn(lst)+'_degree'+strn(degree)+'.txt',x2c-x1c,y2c-y1c,dvx,dvy,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
	   
		;~ forprint, TEXTOUT= '/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,format='(10(f, 4X))', /NOCOMMENT 
		;~ forprint, TEXTOUT= gaussian+name+'IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'.txt',x_dis,y_dis,dvx,dvy,mH,m,a,d,raH,decH,format='(10(f, 4X))', /NOCOMMENT 
		forprint, TEXTOUT= gaussian+name+'_IDL_mas_vx_vy_field'+strn(field)+'_chip'+strn(lst)+'_degree'+strn(degree)+'.txt',x_dis,y_dis,dvx,dvy,mH,mK,m,a,d,raH,decH,format='(11(f, 4X))', /NOCOMMENT 
	
	
    
    

	;~ forprint, TEXTOUT= tmp_p+'checking_lits.txt',x2c ,dx , y2c, dy, x1c,dx_gns/0.106,y1c,dy_gns/0.106 ,format='(10(f, 4X))', /NOCOMMENT 

end

