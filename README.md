# Zoccali_astrophotometry
## USED on LETTER
********************************
on branch: NPL058
********************************

at /Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/WHOLE_im/Zoccali_astrophotometry
README
Photometry for WHOLE images

### IDL
1. aling_new.pro. Aligns the images with correl_optimize, put them in a cube and generate a list for offsets.
2. cachos_new.pro. cut off the images from the cube using the offset list generated in aling_new
3. extractpsf_new.pro. Gets the psf over cachos.
4. deepastro_new. pro. Runs starfinder on the reduced aligned cachos (extracted from the cube)

### PYTHON
1. coordinates_on_the_cube.py. Add the offsets of align_newrpro back to the coordinates´ lists
2. cube_lists_alignment_improved?.py.We use this one in stars_selection branch. Also cube_lists_alignment.py can be used
3. _divide_the_liits- Divide this lists of stars in regions that have this same wt and generates a txt with the images that makes up each list
4. _ZPs_and_meanoffset_im1calibrator_WHOLE. Gets the zp of each image using ima1 as a calibrator
5. _mag_and_position_plots. Plots mag vs drag and mag vs dx and dy for all list for each chip. Each uncertainty is the propagation error formula for the mean. Only data with al least four paintings. So NO list I.
 Also returns some txt files and .json files with mean mag and mag (text files) and x and y position for the common stars (stars that are presentes within 1.5 pixel apart in all lists)
6. _make_txt_of_chips. Make a text file for each chip with 'ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b' (dx and dy are in arcsec).

7. DO NOT USE

8. _aling_chips_brick.pro. Aligns the chips with SIRIUS.This is instead of SIRIUS_alignment
9. calibrate_brick.pro Calibrate the potometry comparing with SIRIUS 
10. _photometry_plotting.py. Plots dmag vs mag, ZPs, ZPs uncertainties in bins, and x diff with SIRIUS 
10. aa_GNS_brick.py. Alings GNs with NPL058, out of brick. You have to choose field and chip forn GNS and it generates a trasnformed list of Zoc data to being used in the next script
11. aligment_with_GNS.pro. It does align with the transdformed GNS using IDL code.  after this one you can jump to Gassuian_fit reposositoy.
12. jk_aligment_with_GNS.pro. Repit the alignment process but getting rid of one of the common stars each time. i.e. repit the alignment one time for each common star, and store the list of stars in a list.
13. jk_error_alignment.py. Compute the mu and sigma for the position of the stars after the differentes alignments (the ones done in the previous step). The sigma is the error of the aligment.
