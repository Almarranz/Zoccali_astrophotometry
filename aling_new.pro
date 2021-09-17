PRO aling_new, band,exptime

;probando cambio para script_git dir.

folder='im_jitter_NOgains/'
py_pruebas='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/photometry/pruebas/'
indir = '/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/07_Cleancubes/058_'+band+'/dit_'+strn(exptime)+'/'+folder
mask_path='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/04_Makemask/058_'+band+'/dit_'+strn(exptime)+'/im/'
outdir='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/07.1_Reduce_aligned/058_'+band+'/dit_'+strn(exptime)+'/'+folder
file_txt ='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/06_Reduce/058_'+band+'/dit_'+strn(exptime)+'/'+folder+'/'
;outdir=py_pruebas

mask_path='/Users/alvaromartinez/Desktop/PhD/HAWKI/The_Brick/04_Makemask/058_'+band+'/dit_'+strn(exptime)+'/im/'

wx=1000
wy=1000

readcol, file_txt + 'reduced_dit'+strn(exptime)+'.txt',files, Format ='A',COUNT=count
cube_d=fix(count/4)
print, 'Cube dimension',cube_d


s=500
sizex_new = 2048 + s
sizey_new = 2048 + s
total = fltarr(sizex_new,sizey_new,cube_d)
mask_cube=fltarr(sizex_new,sizey_new,cube_d);---------------MASK

for chip=1,4 do begin

  for i=1, cube_d do begin


	cube_H=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXT=0,h0)
	im1=readfits(indir+'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXT=1,header)
	mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits')
	
	sz=size(im1)
	print,sz[2]
	ra=sxpar(h0,'RA')
	dec=sxpar(h0,'DEC')

	x_off = strsplit(h0[601],'HIERARCH ESO SEQ CUMOFFSETX = ', ESCAPE = '/', /extract)
	y_off = strsplit(h0[602],'HIERARCH ESO SEQ CUMOFFSETY = ', ESCAPE = '/', /extract)

	x_off=float(x_off[0])
	y_off=float(y_off[0])
	print,x_off, y_off


	CRPIX1=sxpar(header,'CRPIX1');pixel
	CRPIX2=sxpar(header,'CRPIX2')

	CRVAL1=sxpar(header,'CRVAL1');degrees
	CRVAL2=sxpar(header,'CRVAL2')

	sxaddpar, header, 'CRPIX1',CRPIX1+s/2-x_off
	sxaddpar, header, 'CRPIX2',CRPIX2+s/2-x_off

	;sxaddpar, header, 'CRVAL1',CRVAL1+250*0.106/3600 ;este parametro lo cambia automaticamente al cambiar CRPIX, pero no viceversa
	;sxaddpar, header, 'CRVAL2',CRVAL2+250*0.106/3600

	;sxaddpar, h0, 'RA',ra+250*0.106/3600 ; should I change this one??? 
	;sxaddpar, h0, 'DEC',dec+250*0.106/3600


 if i eq 1 then begin

	total[s/2-x_off:s/2+sz[2]-1-x_off,s/2-y_off:s/2+sz[2]-1-y_off,i-1]=im1
    mask_cube[s/2-x_off:s/2+sz[2]-1-x_off,s/2-y_off:s/2+sz[2]-1-y_off,i-1]=mask
    
	small_ref = fltarr(wx,wy)
	a = sizex_new/2-wx/2
	b = sizey_new/2-wy/2
	small_ref = im1[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
	small_ref[0:2,*]=1
	small_ref[wx-2:wx,*]=1
	small_ref[*,0:2]=1
	small_ref[*,wy-2:wy]=1

	openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun
	printf, outp, format='(7f13.3)', x_off, y_off, 0, 0
	free_lun, outp
    
    cube_ref=cube_H
    h0_ref=h0
    header_ref=header
    
	


 endif else begin

	   small_new = fltarr(wx,wy)
	   small_new = im1[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
	   small_new[0:2,*]=1
	   small_new[wx-2:wx,*]=1
	   small_new[*,0:2]=1
	   small_new[*,wy-2:wy]=1
	   
	   print,'comparando con.....','im'+strn(i)+'_chip'+strn(chip)+'_dit'+strn(exptime)
	   correl_optimize, small_ref, small_new, x_off_s, y_off_s, MAGNIFICATION=4, /NUMPIX
       print, 'Offsets from correlation: ' + strn(x_off_s) + ', ' + strn(y_off_s)
       
       openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun, /APPEND
       printf, outp, format='(7f13.3)', x_off, y_off, x_off_s, y_off_s
       free_lun, outp
       
       total[s/2-x_off+x_off_s:sz[1]+s/2-1-x_off+x_off_s,s/2-y_off+y_off_s:sz[1]+s/2-1-y_off+y_off_s,i-1]=im1
       mask_cube[s/2-x_off+x_off_s:sz[1]+s/2-1-x_off+x_off_s,s/2-y_off+y_off_s:sz[1]+s/2-1-y_off+y_off_s,i-1]=mask
       
      
 endelse
 
       writefits, outdir+'cube_chip'+strn(chip)+'_canvas.fits',cube_ref,h0_ref
	   writefits, outdir+'cube_chip'+strn(chip)+'_canvas.fits',total,header_ref,/app
       
       writefits, outdir+'mask_cube_chip'+strn(chip)+'_canvas.fits',cube_ref,h0_ref
	   writefits, outdir+'mask_cube_chip'+strn(chip)+'_canvas.fits',mask_cube,header_ref,/app

 endfor
print, '########### Done with chip',chip,'###########'

endfor 

end
