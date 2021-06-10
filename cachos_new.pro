PRO cachos_new,band,exptime

;folder='im_dark/'
;folder='im_jitter_gains/'
folder='im_jitter_NOgains/'
;folder ='im_sky_ESOReflex/'



s=500
sx=[2048,2048]

pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
outdir='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder
tmp = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+strn(exptime)+'/'+folder+'tmp/



for chip=1, 4 do begin

readcol, indir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', x_off, y_off,x_off_s,y_off_s,Format ='A,A,A,A',COUNT=count

cube_d=count
print, count

x_off=float(x_off)
y_off=float(y_off)
x_off_s=float(x_off_s)
y_off_s=float(y_off_s)



for i=1,cube_d do begin ;###############





datos=readfits(indir+'cube_chip'+strn(chip)+'_canvas.fits', EXTEN_NO=0,header0)
cube = readfits(indir+'cube_chip'+strn(chip)+'_canvas.fits', EXTEN_NO=1,header1)
l=size(cube)
lado=l[1]
x_pix=sxpar(header1,'CRPIX1')
y_pix=sxpar(header1,'CRPIX2')

sxaddpar, header1, 'CRPIX1 ',x_pix+x_off[i-1]-x_off_s[i-1]-s/2
sxaddpar, header1, 'CRPIX2 ',y_pix+y_off[i-1]-y_off_s[i-1]-s/2

im=cube[s/2-x_off[i-1]+x_off_s[i-1]:sx[1]+s/2-1-x_off[i-1]+x_off_s[i-1],s/2-y_off[i-1]+y_off_s[i-1]:sx[1]+s/2-1-y_off[i-1]+y_off_s[i-1],i-1]
print,size(im)
writefits, tmp+ 'im'+strn(i)+'chip'+strn(chip)+'_'+band+'dit_'+strn(exptime)+'_cacho.fits',datos,header0
writefits, tmp+ 'im'+strn(i)+'chip'+strn(chip)+'_'+band+'dit_'+strn(exptime)+'_cacho.fits', im,header1,/app 

endfor
print, '############## Done witn chip',chip,'##############'
print, '############## Done witn chip',chip,'##############'
print, '############## Done witn chip',chip,'##############'
endfor



END
