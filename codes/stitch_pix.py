from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

home = os.getenv('HOME')  # does not have a trailing slash

taffy_products = home + '/Desktop/ipac/taffy_lzifu/products/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in cube that is edited each time a pixel is fixed
    patched_filename = 'Taffy_2_comp_patched.fits'
    patched_cube_hlist = fits.open(taffy_extdir + patched_filename)
    
    # read in most recent pixel run
    pix_hdu = fits.open(taffy_products + patched_filename.replace('_patched.fits', '.fits'))

    curr_pix_x,curr_pix_y = 31,36
    arr_x,arr_y = curr_pix_y-1,curr_pix_x-1
    # 2 lines above arent' quite pythonic syntax but 
    # it is easier to read coordinates this way.

    print 'Stitching pixel (ds9 coords)', curr_pix_x, curr_pix_y
    
    # get extanmes
    extnames = []
    
    if '1' in patched_filename:
        for i in range(1,37):
            extnames.append(pix_hdu[i].header['EXTNAME'])
        total_ext = 36
        line_shape = (2, 58, 58)
    elif '2' in patched_filename:
        for i in range(1,39):
            extnames.append(pix_hdu[i].header['EXTNAME'])
        total_ext = 38
        line_shape = (3, 58, 58)

    # define shapes because the extensions in the fits file have diff shapes
    blue_shape = (2227, 58, 58)
    red_shape = (2350, 58, 58)

    # loop over each extension and replace nan data with new fit data
    for i in range(total_ext):
        shape = patched_cube_hlist[extnames[i]].data.shape
        if (shape == blue_shape) or (shape == red_shape) or (shape == line_shape):
            patched_cube_hlist[extnames[i]].data[:,arr_x,arr_y] = pix_hdu[extnames[i]].data[:,arr_x,arr_y]
            #patched_cube_hlist[extnames[i]].data[:,arr_x,arr_y] = np.zeros(shape[0])
            # use this second line above to replce the line fits 
            # with zeros if LZIFU refuses to give you an answer for a pixel.
            # comment out the first line if you use the second
        elif shape == (58, 58):
            patched_cube_hlist[extnames[i]].data[arr_x,arr_y] = pix_hdu[extnames[i]].data[arr_x,arr_y]

    patched_cube_hlist.writeto(taffy_extdir + patched_filename, clobber=True, output_verify='fix')

    patched_cube_hlist.close()
    pix_hdu.close()

    print "Done."
    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)