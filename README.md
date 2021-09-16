# Zoccali_astrophotometry
********************************
on branch: iMac
********************************

at /Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/WHOLE_im/Zoccali_astrophotometry
README
Photometry for WHOLE images

#####IDL#####
1. aling_new.pro. Aligns the images with correl_optimize, put them in a cube and generate a list for offsets.
2. cachos_new.pro. cut off the images from the cube using the offset list generated in aling_new
3. extractpsf_new.pro. Gets the psf over cachos.
4. deepastro_new. pro. Runs starfinder on the reduced aligned cachos (extracted from the cube)

#####PYTHON#####
1. coordinates_on_the_cube.py. Add the offsets of align_newrpro back to the coordinates´ lists
2. cube_lists_alignment_improved?.py.We use this one in stars_selection branch. Also cube_lists_alignment.py can be used
3. _divide_the_liits- Divide this lists of stars in regions that have this same wt and generates a txt with the images that makes up each list
4. _ZPs_and_meanoffset_im1calibrator_WHOLE. Gets the zp of each image using ima1 as a calibrator
5. _mag_and_position_plots. Plots mag vs drag and mag vs dx and dy for all list for each chip. Each uncertainty is the propagation error formula for the mean. Only data with al least four paintings. So NO list I.
 Also returns some txt files and .json files with mean mag and mag (text files) and x and y position for the common stars (stars that are presentes within 1.5 pixel apart in all lists)
6. _make_txt_of_chips. Make a text file for each chip with 'ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b' (dx and dy are in arcsec).

7. DO NOT USE

################## WARNING ######################## 

DO NOT USE 8_SIRIUS_alignment.py NO READY yet

###############################################

8. _aling_chips_brick.pro. Aligns the chips with SIRIUS.This is instead of SIRIUS_alignment
9. calibrate_brick.pro Calibrate the potometry comparing with SIRIUS 
10. _photometry_plotting.py. Plots dmag vs mag, ZPs, ZPs uncertainties in bins, and x diff with SIRIUS 
11. aligment_with_GNS.pro. It does align with GNS using IDL code. It seems to be better some how that python alygnmnt. after this you can go to Gaussian_fir repository
######################################################################################################
WARNING for here on is python alignment. Eventhoug they work, it is better to stop here. Use carefully
######################################################################################################
11. _aligment_with_GNS.py'you have to open the lists on Aladin locate a two common stars and copy and paste their X and Y coord$
12. _proper_motions_plots.py With the scritp like this, making the aligment at previos script using only stars h-ks>1.3 from GN$
13. _align_and_prmot.py Combine the alignment with GNS and the proper motion analisys in a single script. You can set where you$


################# NOTE ##########################

2. cube_lists_alignment. can be improved by looping the degree 1 polynomial over list_E and then apply the kx ky to the whole list.
When looping with distancia =1 the improvement is huge from the first loop to the second one. This way I can reach the same number of aligment 
stars using a distancia = 1 that I got using a distancia = 2. Anyway, the final improvent for the uncertainty in the position doesnt seem to be
significt.
A test script is (almost?) ready at:
Desktop/PhD/HAWK/The_Brick/photometry/WHOLE_im/Zoccali_astrophotometry/cube_lists_alignment_improved?.py.

7_plotting_commons is not useful anymora. A plot of common stars is done on 11.

out_of_Brick_alignment. This uses the same procedure for alignment a GNS list (of field16) that is out of the brick.
