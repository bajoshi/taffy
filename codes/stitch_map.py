from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products_big_cube_velsort/'
taffy_data = '/Volumes/Bhavins_backup/ipac/TAFFY/data/'
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

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
            elif (shape == (3, 58, 38)) or (shape == (2, 58, 58)):
                # This elif part was put in because the single comp fits have line cubes which are shaped
                # as (2, 58, 58) and the two comp fits have lines cubes shaped as (3, 58, 58).
                # What I'm trying to do here is to first take all the 3 values that are present in the 
                # stitched cube at a given pixel and replace them by -9999.0. Then I keep the first one 
                # (i.e. zeroth index) as -9999.0 and replace the next two with data from the single comp
                # fit which can now be done with the correct shape. 
                stitched_cube[extnames[k]].data[:,arr_x,arr_y] = np.ones(shape[0]) * -9999.0
                stitched_cube[extnames[k]].data[1:,arr_x,arr_y] = single_comp[extnames[k]].data[:,arr_x,arr_y]
            elif shape == (58, 58):
                stitched_cube[extnames[k]].data[arr_x,arr_y] = single_comp[extnames[k]].data[arr_x,arr_y]

    return None

if __name__ == '__main__':
    
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

    # loop over all pixels 
    # My goal here is to keep the shapes of the stitched cube 
    # extensions to be the same as those of the two comp fit results.
    for i in range(58):
        for j in range(58):

            if [i,j] in nan_single_comp_list:
                single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                continue

            if single_idx[i,j]:
                single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                continue

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    single_comp_stitch(extnames, stitched_cube, single_comp, i, j, blue_shape, red_shape, line_shape)
                    continue

                else:
                    # loop over each extension
                    for k in range(38):
                        shape = stitched_cube[extnames[k]].data.shape
                        if (shape == blue_shape) or (shape == red_shape) or (shape == line_shape):
                            stitched_cube[extnames[k]].data[:,i,j] = two_comp[extnames[k]].data[:,i,j]
                        elif shape == (58, 58):
                            stitched_cube[extnames[k]].data[i,j] = two_comp[extnames[k]].data[i,j]

    stitched_cube.writeto(savedir + 'stitched_cube.fits', clobber=True, output_verify='fix')
    stitched_cube.close()

    single_comp.close()
    two_comp.close()

    sys.exit(0)