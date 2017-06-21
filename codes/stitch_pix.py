from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

home = os.getenv('HOME')  # does not have a trailing slash

taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products/'
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in cube that is edited each time a pixel is fixed
    big_cube_hlist = fits.open(taffy_extdir + 'big_cube_2_comp_velsort.fits')
    
    # read in most recent pixel run
    pix_hdu = fits.open(taffy_products + 'Taffy_2_comp.fits')

    curr_pix_x,curr_pix_y = 33,35
    arr_x,arr_y = curr_pix_y-1,curr_pix_x-1
    # 2 lines above arent' quite pythonic syntax but 
    # it is easier to read coordinates this way.

    print 'Stitching pixel (ds9 coords)', curr_pix_x, curr_pix_y
    
    # get extanmes
    extnames = []
    
    for i in range(1,39):
        extnames.append(pix_hdu[i].header['EXTNAME'])

    # define shapes because the extensions in the fits file have diff shapes
    blue_shape = (2227, 58, 58)
    red_shape = (2350, 58, 58)
    line_shape = (3, 58, 58)

    # loop over each extension and replace nan data with new fit data
    for i in range(38):
        shape = big_cube_hlist[extnames[i]].data.shape
        if (shape == blue_shape) or (shape == red_shape) or (shape == line_shape):
            big_cube_hlist[extnames[i]].data[:,arr_x,arr_y] = pix_hdu[extnames[i]].data[:,arr_x,arr_y]
            #big_cube_hlist[extnames[i]].data[:,arr_x,arr_y] = np.zeros(shape[0])
            # use this second line above to replce the line fits 
            # with zeros if LZIFU refuses to give you an answer for a pixel.
        elif shape == (58, 58):
            big_cube_hlist[extnames[i]].data[arr_x,arr_y] = pix_hdu[extnames[i]].data[arr_x,arr_y]

    big_cube_hlist.writeto(taffy_extdir + 'big_cube_2_comp_velsort.fits', clobber=True, output_verify='fix')

    big_cube_hlist.close()
    pix_hdu.close()

    print "Done."
    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)