from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products =  home + '/Desktop/ipac/taffy_lzifu/products_work/'
taffy_data =  home + '/Desktop/ipac/taffy_lzifu/data/'
taffy_extdir =  home + '/Desktop/ipac/taffy_lzifu/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir =  home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

def single_comp_stitch(extnames, stitched_cube, single_comp, arr_x, arr_y, blue_shape, red_shape, line_shape):

    for k in range(38):
        if 'COMP1' in extnames[k]:
            stitched_cube[extnames[k]].data[:,arr_x,arr_y] = single_comp[extnames[k]].data[:,arr_x,arr_y]
        elif 'COMP2' in extnames[k]:
            if 'B_LINE_COMP2' in extnames[k]:
                stitched_cube[extnames[k]].data[:,arr_x,arr_y] = np.ones(blue_shape[0]) * -9999.0
            elif 'R_LINE_COMP2' in extnames[k]:
                stitched_cube[extnames[k]].data[:,arr_x,arr_y] = np.ones(red_shape[0]) * -9999.0
        else:
            shape = stitched_cube[extnames[k]].data.shape
            if (shape == blue_shape) or (shape == red_shape):
                stitched_cube[extnames[k]].data[:,arr_x,arr_y] = single_comp[extnames[k]].data[:,arr_x,arr_y]
            elif (shape == (3, 58, 58)) or (shape == (2, 58, 58)):
                # This elif part was put in because the single comp fits have line cubes which are shaped
                # as (2, 58, 58) and the two comp fits have lines cubes shaped as (3, 58, 58).
                # What I'm trying to do here is to first take all the 3 values that are present in the 
                # stitched cube at a given pixel and replace them by -9999.0. Then I keep the last one 
                # (i.e. zeroth index) as -9999.0 and replace the first two with data from the single comp
                # fit which can now be done with the correct shape. 
                stitched_cube[extnames[k]].data[:,arr_x,arr_y] = np.ones(shape[0]) * -9999.0
                stitched_cube[extnames[k]].data[:2,arr_x,arr_y] = single_comp[extnames[k]].data[:,arr_x,arr_y]
            elif shape == (58, 58):
                stitched_cube[extnames[k]].data[arr_x,arr_y] = single_comp[extnames[k]].data[arr_x,arr_y]

    return stitched_cube

def two_comp_stitch(extnames, stitched_cube, two_comp, arr_x, arr_y, blue_shape, red_shape, line_shape):

    # loop over each extension
    for k in range(38):
        shape = stitched_cube[extnames[k]].data.shape
        if (shape == blue_shape) or (shape == red_shape) or (shape == line_shape):
            stitched_cube[extnames[k]].data[:,arr_x,arr_y] = two_comp[extnames[k]].data[:,arr_x,arr_y]
        elif shape == (58, 58):
            stitched_cube[extnames[k]].data[arr_x,arr_y] = two_comp[extnames[k]].data[arr_x,arr_y]

    return stitched_cube

def get_nan_arr():

    # define the list of spaxels which are nan in the 2-comp fit but
    # the data only shows a single comp so the code stitches the 1-comp fit
    # these are ds9 referenced coordinates
    nan_single_comp_list = [[33, 15],[35, 13],[34, 12],[35, 12],[36, 12],[34, 11],\
    [35, 11],[36, 11],[37, 11],[38, 11],[35, 10],[36, 10],[35, 9],[36, 9],[37, 9],\
    [35, 8],[40, 8],[41, 9],[39, 5],[39, 6],[40, 5],[40, 6],[41, 5],[40, 4],\
    [39, 3],[38, 4],[45, 12],[44, 11],[44, 10],[45, 10],[45, 14],[29, 4],[29, 9],\
    [30, 11],[26, 9],[25, 18],[26, 16],[26, 17],[34, 56], [35, 56],[21, 56],\
    [17, 53],[17, 54],[19, 54],[19, 55],[36, 51],[36, 52],[36, 53],[36, 54],\
    [36, 55],[36, 56],[37, 48],[37, 49],[37, 50],[37, 51],[37, 52],[37, 53],\
    [37, 54],[37, 55],[37, 56],[37, 57],[9, 39],[10, 39],[11, 39],[13, 38],\
    [9, 38],[10, 38],[11, 38],[12, 38],[10, 37],[11, 37],[12, 37],[10, 36],\
    [11, 36],[9, 36],[9, 37],[24, 17],[42, 10],[19, 54],[20, 55],[35, 54],\
    [12, 36],[13, 36],[21, 18],[22, 18],[21, 17],[22, 17],[23, 17],[21, 14],\
    [22, 14],\
    # adding some some spaxels which are really noisy and my 2-comp fit
    # while satisfying all conditions is actually bogus. So these will 
    # have the 1-comp stitched in.
    [28, 20],[29, 20],[30, 20],[29, 21],[30, 21],[29, 22],[30, 22]]

    nan_arr_list_x = []
    nan_arr_list_y = []

    for u in range(len(nan_single_comp_list)):

        nan_arr_list_x.append(nan_single_comp_list[u][1] - 1)
        nan_arr_list_y.append(nan_single_comp_list[u][0] - 1)

    nan_single_comp_arr = zip(nan_arr_list_x, nan_arr_list_y)

    return nan_single_comp_arr

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    #plt.imshow(single_idx, origin='lower', cmap='Greys')
    #plt.show()
    #plt.imshow(diffmean_idx, origin='lower', cmap='Greys')
    #plt.show()
    #plt.imshow(diffstd_idx, origin='lower', cmap='Greys')
    #plt.show()
    #plt.imshow(diffboth_idx, origin='lower', cmap='Greys')
    #plt.show()
    #sys.exit(0)

    # read in both 1 and 2 comp fitting results
    single_comp = fits.open(taffy_extdir + 'Taffy_1_comp_patched.fits')
    two_comp = fits.open(taffy_extdir + 'Taffy_2_comp_patched.fits')

    stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

    # loop over all pixels and stitch them according 
    # to whether they need a one or two comp fit
    # get extnames
    extnames = []
    
    for u in range(1,39):
        extnames.append(two_comp[u].header['EXTNAME'])

    # define shapes because the extensions in the fits file have diff shapes
    blue_shape = (2227, 58, 58)
    red_shape = (2350, 58, 58)
    line_shape = (3, 58, 58)

    nan_single_comp_arr = get_nan_arr()

    # set up counters
    nan_count = 0
    onecomp_count = 0
    invalid_fit_count = 0
    twocomp_count = 0

    # loop over all pixels 
    # My goal here is to keep the shapes of the stitched cube 
    # extensions to be the same as those of the two comp fit results.
    for i in range(58):
        for j in range(58):

            if (i,j) in nan_single_comp_arr:
                stitched_cube = single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                nan_count += 1
                continue

            if single_idx[i,j]:
                stitched_cube = single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                onecomp_count += 1
                continue

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    stitched_cube = single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                    invalid_fit_count += 1

                elif comp1_inv_idx[i,j] and ~comp2_inv_idx[i,j]:
                    stitched_cube = two_comp_stitch(extnames, stitched_cube, two_comp, i, j, blue_shape, red_shape, line_shape)
                    twocomp_count += 1

                elif ~comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    stitched_cube = two_comp_stitch(extnames, stitched_cube, two_comp, i, j, blue_shape, red_shape, line_shape)
                    twocomp_count += 1

                else:
                    stitched_cube = two_comp_stitch(extnames, stitched_cube, two_comp, i, j, blue_shape, red_shape, line_shape)
                    twocomp_count += 1

    print "total nan spaxels in two comp result", nan_count
    print "total single comp classified spaxels", onecomp_count
    print "total invalid fit spaxels", invalid_fit_count
    print "total two comp classified spaxels", twocomp_count

    stitched_cube.writeto(taffy_extdir + 'stitched_cube.fits', clobber=True, output_verify='fix')
    stitched_cube.close()

    single_comp.close()
    two_comp.close()

    print "Done Stitching."
    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)