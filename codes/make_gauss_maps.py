from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products_big_cube_velsort/'

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in lzifu output file
    h = fits.open(taffy_products + 'big_cube_2_comp_velsort.fits')

    # assign line arrays
    # Each one of these cubes has the shape 
    # Blue: (2227,58,58)
    # Red: (2350,58,58)
    b_cont = h['B_CONTINUUM'].data
    b_resid = h['B_RESIDFIT'].data
    b_line = h['B_LINE'].data
    b_line_comp1 = h['B_LINE_COMP1'].data
    b_line_comp2 = h['B_LINE_COMP2'].data

    r_cont = h['R_CONTINUUM'].data
    r_resid = h['R_RESIDFIT'].data
    r_line = h['R_LINE'].data
    r_line_comp1 = h['R_LINE_COMP1'].data
    r_line_comp2 = h['R_LINE_COMP2'].data

    # I'm choosing to fit a gaussian to just hte [OIII]5007 line.
    # This is because I think this line has high enough SNR 
    # and also is not contaminated by any other lines. I would've 
    # chosen H-alpha but that has the problem of contamination 
    # from the [NII] satellite lines.

    # For the [OIII]5007 line --
    # It appears that both the northern and southern galaxies have 
    # double peaked profiles. The southern galaxy more so than the 
    # northern one. The southern galaxy has double peaked lines all
    # above its nucleus and not so much below the nucleus. The S.Nuc. 
    # itself seems to have really broad lines.
    # The northern galaxy mostly has double peaked lines closer to 
    # the HII region. It might also have double peaked lines near 
    # the nucleus and in the NW region but this isn't obvious.
    # 

    # loop over all spaxels and fit a gaussian to the individual line fits
    for i in range(58):
        for j in range(58):


    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)