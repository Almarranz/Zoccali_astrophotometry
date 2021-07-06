# Zoccali_astrophotometry
********************************
on branch: MAIN
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
2. cube_lists_alignment.
3. _divide_the_liits- Divide this lists of stars in regions that have this same wt and generates a txt with the images that makes up each list
4. _ZPs_and_meanoffset_im1calibrator_WHOLE. Gets the zp of each image using ima1 as a calibrator
5. _mag_and_position_plots. Plots mag vs drag and mag vs dx and dy for all list for each chip. Each uncertainty is the propagation error formula for the mean. Only data with al least four paintings. So NO list I.
 Also returns some txt files and .json files with mean mag and mag (text files) and x and y position for the common stars (stars that are presentes within 1.5 pixel apart in all lists)
6. _make_txt_of_chips. Make a text file for each chip with 'ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b' (dx and dy are in degrees)
7. _plotting_commmons.  Make plots of the stars on the overplottings areas of both surveys

################## WARNING ######################## 

DO NOT USE 8_SIRIUS_alignment.py NO READY yet

###############################################

8. _aling_chips_brick.pro. Aligns the chips with SIRIUS.This is instead of SIRIUS_alignment
9. calibrate_brick.pro Calibrate the potometry comparing with SIRIUS 
10. _photometry_plotting.py. Plots dmag vs mag, ZPs, ZPs uncertainties in bins, and x diff with SIRIUS 
11. _aligment_with_GNS.py Align with GNS 2deg polynomial and makes some plots. I have tried the aligments with the method of initial offset, the clinking on the same star one. The outcome is the same'

################# NOTE ##########################

2. cube_lists_alignment. can be improved by looping the degree 1 polynomial over list_E and then apply the kx ky to the whole list.
When looping with distancia =1 the improvement is huge from the first loop to the second one. This way I can reach the same number of aligment 
stars using a distancia = 1 that I got using a distancia = 2. Anyway, the final improvent for the uncertainty in the position doesnt seem to be
significt.
A test script is (almost?) ready at:
Desktop/PhD/HAWK/The_Brick/photometry/WHOLE_im/Zoccali_astrophotometry/cube_lists_alignment_improved?.py
