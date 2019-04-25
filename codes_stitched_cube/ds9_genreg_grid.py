from __future__ import division

import numpy as np

import sys
import os

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

def genreg_orig():

    # Open file to write regions file to
    file_cat = open(taffy_extdir + 'HI_grid.reg', 'wa')

    # Write header string
    hdrstr1 = "# Region file format: DS9 version 4.1" 
    hdrstr2 = "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" "
    hdrstr3 = "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1"

    file_cat.write(hdrstr1 + '\n' + hdrstr2 + hdrstr3 + '\n')

    # Needs some basic coordinate info and grid info to know how to draw the grid
    # MY IFU data has dimensions 58x58.
    # 1. I'm going to make this grid in pixel space so it is easier to impose on the IFU grid.
    # 2. Then I'll use ds9 to convert the resulting regions file .
    # to wcs and impose that on the HI cube. 
    # 3. Convert that regions file back to image space 
    # (of HI cube) to work with Polygon and extract spectra.

    # DS9 starts counting at (1,1) which it considers 
    # to be the center of the lower left pixel.
    # Therefore, the center of the starting cell 
    # at lower left will be at an edge which is (6.5,6.5).
    starting_x = 6.5  # starting x lower left
    starting_y = 6.5  # starting y lower left
    width = 12
    height = 12

    # Loop over all grid cells
    for i in range(5):  # x
        for j in range(5):  # y

            # Only for starting cell
            if i==0 and j==0:
                center_x = starting_x
                center_y = starting_y

            # Generate string to write and write
            write_str = str(center_x) + ',' + str(center_y) + ',' + str(width) + ',' + str(height)
            file_cat.write('image;box(' + write_str + ',0) # color=red width=2;' + '\n')

            # Step over to next cell over
            if j==3:  # i.e., if you're about to go to the last column then cell size changes
                center_x += 11
                width = 10
            else:
                center_x += width
                width = 12

        # ----------------------------
        # Update coordinates for next row up
        center_x = starting_x

        if i==3:  # i.e., if you're about to go to the last row then the cell size changes
            center_y += 11
            height = 10
        else:
            center_y += height
            height = 12
    
    # Close file so that it can actually be written
    file_cat.close()

    return None

def genreg():

    # Open file to write regions file to
    file_cat = open(taffy_extdir + 'HI_grid_new.reg', 'wa')

    # Write header string
    hdrstr1 = "# Region file format: DS9 version 4.1" 
    hdrstr2 = "global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" "
    hdrstr3 = "select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1"

    file_cat.write(hdrstr1 + '\n' + hdrstr2 + hdrstr3 + '\n')

    # Needs some basic coordinate info and grid info to know how to draw the grid
    # MY IFU data has dimensions 58x58.
    # 1. I'm going to make this grid in pixel space so it is easier to impose on the IFU grid.
    # 2. Then I'll use ds9 to convert the resulting regions file .
    # to wcs and impose that on the HI cube. 
    # 3. Convert that regions file back to image space 
    # (of HI cube) to work with Polygon and extract spectra.

    # DS9 starts counting at (1,1) which it considers 
    # to be the center of the lower left pixel.
    # Therefore, the center of the starting cell 
    # at lower left will be the center of the (4,4) pixel.
    starting_x = 4  # starting x lower left
    starting_y = 4  # starting y lower left
    width = 7
    height = 7

    # Loop over all grid cells
    for i in range(8):  # x
        for j in range(8):  # y

            # Only for starting cell
            if i==0 and j==0:
                center_x = starting_x
                center_y = starting_y

            # Generate string to write and write
            write_str = str(center_x) + ',' + str(center_y) + ',' + str(width) + ',' + str(height)
            file_cat.write('image;box(' + write_str + ',0) # color=red width=2;' + '\n')

            # Step over to next cell over
            center_x += width

        # ----------------------------
        # Update coordinates for next row up
        center_x = starting_x
        center_y += height
    
    # Close file so that it can actually be written
    file_cat.close()

    return None

def main():

    # Generate the "grid" -- the grid is essentially a set of regions
    # that are box shaped and look like a grid, saved to a single 
    # ds9 regions file.
    # Check the file by loading into DS9 and making sure it looks like
    # a grid. 
    genreg()

    return None

if __name__ == '__main__':
    main()
    sys.exit(0)