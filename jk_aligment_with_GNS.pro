pro jk_aligment_with_GNS,field,lst,degree
;lst  can be 1 to 4 (refers to the chip on GNS fields)
;~ field=20 ; fields can be 3 or 20 or 16 (refers to GNS fields)
;~ Esto es una prueba para git pull
;~ NOTE:
;~ degree: degree of the poly. fit. (1 and 2=2. 1,2, and 3=3)

;~ IMPORTAT!!!!!!!!!!!!!!
;~ Check if Z1 is in front of the lists names


band='H'



field=16
lst=3

exptime=10
folder='im_jitter_NOgains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field'+strn(field)+'/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp_bs/'
gaussian='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'Gaussian_fit/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'/results/'
jackknive='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/jack_knive/'
;~ name='aa_NPL058_'
;~ name='Z1_aa_NPL058'
zone='Z1' ;Zone refers to the little squares the same size of the one of the brick on zone B; Z1 is left, z2 is right
name=strn(zone)+'_aa_NPL058'
print, name

markstars=0

;~ readcol, GNS + 'cat_Ban_'+strn(field)+'_'+strn(lst)+'.txt',x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, Format='A,A,A,A,A,A,A,A,A,A,A,A,A,A';,SKIPLINE = 1
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

H_Ks=where((mH-mK) gt 1.3)

raH=raH[H_Ks]
decH=decH[H_Ks]
mH=mH[H_Ks]
mK=mK[H_Ks]
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]




	 

	; dat = ptr_new({X_size: 20, Y_size: 20, Sigma_x: 1.5, Sigma_y: 1.5, Angle: 0.0})
	; map = image_model(xi,yi,f,xsize_quad,ysize_quad,'gaussian', dat)
	; writefits, tmp_path + 'align_sources.fits', map

	 
	
	 
	 nsample = 1477
for i=0, nsample - 1 do begin

	
 
	;~ x_gns_jk = x_gns[indb] ; after removing	
	;~ y_gns_jk=  y_gns[indb]
	
	 dmax = 1
	 compare_lists, x_gns, y_gns, x, y, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	 nc = n_elements(subc1)
	 print, 'Found ' + strn(nc) + ' common stars.ALL'
	 print, n_elements(x2c),n_elements(x[subc2]),n_elements(x)
	 
	 
	imax = n_elements(x1c)
	;removing randomly
	;~ indb = RANDOM_SAMPLE(seed, FINDGEN(imax), imax-1)
	
	;~ x_gns_jk = x1c[indb] ; after removing	
	;~ y_gns_jk=  y1c[indb]
	
	x_gns_jk=x1c
	y_gns_jk=y1c
	
	
	remove, i, x_gns_jk
	remove, i, y_gns_jk
	
	print, n_elements(x_gns_jk),n_elements(x1c)
		
	
	
	dmax = 1
	 compare_lists, x_gns_jk, y_gns_jk, x, y, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	 nc = n_elements(subc1)
	 print, 'Found ' + strn(nc) + ' common stars.'
	 print, n_elements(x2c),n_elements(x[subc2]),n_elements(x)
	 
	
	
     ; iterative degree 1 alignment
	 ; ------------------------------
	if (degree eq 1) || (degree eq 2) || (degree eq 3)  then begin
     print, '##############################################'
	 print, 'Iteration',i+1,'Degree 1 alignment.'
	 print, '##############################################' 
     count=0
     comm=[]
     it=0
     lim_it=1
	 while count lt lim_it do begin
	  it=it+1
	 
	  polywarp, x_gns_jk[subc1], y_gns_jk[subc1], x[subc2], y[subc2], 1, Kx, Ky
	  print, Kx
	  print, Ky		
	  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y
	  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y
	  ;~ openw, out, jackknive + 'x_new_y_new_jackknife_deg1_cut_grid_im'+ strtrim(string(image+1, format='(I02)'), 2)+ '_'+ strn(i) +'.txt', /get_lun ;2013
	  if degree eq 1 then begin
		openw, out, jackknive + name+'_jackknife_degree_'+strn(degree)+'_'+ strn(i+1) +'.txt', /get_lun ;2013
		for k = 0, n_elements(xi) -1 do begin
			printf, out, xi[k], yi[k], format='(2(F, 4X))';FORMAT='(3f13.5)'
		endfor
		free_lun, out
	  endif 
	  compare_lists, x_gns_jk, y_gns_jk, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	  nc = n_elements(subc1)
	  print, 'Iteration ' + strn(i)
	  print, 'Found ' + strn(nc) + ' common stars. after iteration ',strn(i)
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
	  
	  polywarp, x_gns_jk[subc1], y_gns_jk[subc1], x[subc2], y[subc2], 2, Kx, Ky
	  print, Kx
	  print, Ky
	  xi = Kx[0,0] + Kx[0,1]*x + Kx[1,0]*y + Kx[1,1]*x*y + Kx[0,2]*x^2 + Kx[1,2]*x^2*y + Kx[2,2]*x^2*y^2 + Kx[2,0]*y^2 + Kx[2,1]*y^2*x 
	  yi = Ky[0,0] + Ky[0,1]*x + Ky[1,0]*y + Ky[1,1]*x*y + Ky[0,2]*x^2 + Ky[1,2]*x^2*y + Ky[2,2]*x^2*y^2 + Ky[2,0]*y^2 + Ky[2,1]*y^2*x
	  if degree eq 2 then begin
		openw, out, jackknive + name+'_jackknife_degree_'+strn(degree)+'_'+ strn(i+1) +'.txt', /get_lun ;2013
		for k = 0, n_elements(xi) -1 do begin
			printf, out, xi[k], yi[k], format='(2(F, 4X))';FORMAT='(3f13.5)'
		endfor
		free_lun, out
	  endif 
	  compare_lists, x_gns_jk, y_gns_jk, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=dmax, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
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
  
 
endfor 
	
    
    ;~ distancia=1
    ;~ compare_lists, x_gns, y_gns, xi, yi, x1c, y1c, x2c, y2c, MAX_DISTANCE=distancia, SUBSCRIPTS_1=subc1, SUBSCRIPTS_2 = subc2, SUB1 = sub1, SUB2 = sub2
	;~ nc = n_elements(subc1)
    ;~ print, 'Found after '+strn(degree)+' degree alignment ' + strn(nc) + ' common stars.'
    
    ;~ print, N_ELEMENTS(f)
    ;~ f=f[subc2]
    ;~ print,N_ELEMENTS(f)
    
    ;~ if (markstars gt 0) then begin
		;~ dat = ptr_new({X_size: 10, Y_size: 10, Sigma_x: 2, Sigma_y: 2, Angle: 0.0})
		;~ im = image_model(x_gns,y_gns,10^((-mH)/2.5),xsize_ref,xsize_ref,'gaussian', dat)
		;~ ima = image_model(xi,yi,f,xsize_ref,xsize_ref,'gaussian', dat)
	
	
	;~ writefits, tmp_p +'gns'+strn(field)+'_chip' + strn(lst) + '.fits',im

	
	;~ writefits, tmp_p +'zoc_on_GNS_field'+ strn(field) + '.fits', ima
	;~ endif
	
		;~ x_dis=(x2c-x1c)*0.106/4.3*1000
		;~ y_dis=(y2c-y1c)*0.106/4.3*1000
		;Adding velocities uncertanties for x and y directions
		;~ dx_gns=dx_gns[subc1]*0.106*1000
		;~ dy_gns=dy_gns[subc1]*0.106*1000
		
		;~ dx=dx[subc2]*1000 ; uncertanties in Zoc's data are already in arcsec
		;~ dy=dy[subc2]*1000
		
		;~ dvx=sqrt(dx^2+dx_gns^2)/4.3
		;~ dvy=sqrt(dy^2+dy_gns^2)/4.3
    
    ;~ mH=mH[subc1]
    ;~ mK=mK[subc1]
    ;~ m=m[subc2]
    ;~ a=a[subc2]
    ;~ d=d[subc2]
    ;~ raH=raH[subc1]
    ;~ decH=decH[subc1]
    
    
end

