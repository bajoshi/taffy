from __future__ import division

import numpy as np

import sys
import os

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

def genreg():

    # Open file to write regions file to
    file_cat = open(taffy_extdir + 'HI_grid.reg', 'wa')

    # Needs some basic coordinate info and grid info 
    # to know how to draw the grid
    # MY IFU data has dimensions 58x58
    box_size = 5  # e.g., a box_size = 5 will make a 5x5 grid
    center_ra = 
    center_dec = 

    # Loop over all grid cells
    for i in range(box_size*box_size):
        file_cat.write('fk5;box(' + str(current_ra) + ',' + str(current_dec) + ',0.5") # color=red width=1;' + '\n')
    
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