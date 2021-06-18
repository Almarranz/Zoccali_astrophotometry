# Zoccali_astrophotometry
at /Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/WHOLE_im/Zoccali_astrophotometry
README
Photometry for WHOLE images
#####IDL#####
1. aling_new.pro. Aligns the images with correl_optimize, put them in a cube and generate a list for offsets.
2.cachos_new.pro. cut off the images from the cube using the offset list generated in aling_new
3.extractpsf_new.pro. Gets the psf over cachos.
4.deepastro_new. pro. Runs starfinder on the reduced aligned cachos (extracted from the cube)

#####PYTHON#####
1. coordinates_on_the_cube.py. Add the offsets of align_newrpro back to the coordinates´ lists
2. cube_lists_alignment.
3_divide_the_liits- Divide this lists of stars in regions that have this same wt and generates a txt with the images that makes up each list
4_ZPs_and_meanoffset_im1calibrator_WHOLE. Gets the zp of each image using ima1 as a calibrator
5_mag_and_position_plots. Plots mag vs drag and mag vs dx and dy for all list for each chip. Each uncertainty is the propagation error formula for the mean. Only data with al least four paintings. So NO list I.
 Also returns some txt files and .json files with mean mag and mag (text files) and x and y position for the common stars (stars that are presentes within 1.5 pixel apart in all lists)
6_make_txt_of_chips. Make a text file for each chip with 'ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b' (dx and dy are in degrees)
7_plotting_commmons.  Make plots of the stars on the overplottings areas of both surveys
##################NOTE######################## 
DO NOT USE 8_SIRIUS_alignment.py NO READY yet
###############################################
8_aling_chips_brick.pri. Aligns the chips with SIRIUS. 
