from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products_big_cube_velsort/'
taffy_data = '/Volumes/Bhavins_backup/ipac/TAFFY/data/'
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

if __name__ == '__main__':
    
    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    # read in both 1 and 2 comp fitting results
    single_comp = fits.open(taffy_extdir + 'products/Taffy_1_comp.fits')
    two_comp = fits.open(taffy_products + 'big_cube_2_comp_velsort.fits')

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
    #line_shape = (3, 58, 58)

    # loop over all pixels 
    for i in range(58):
        for j in range(58):

            if single_idx[i,j]:
                for k in range(38):
                    if 'COMP1' in extnames[k]:
                        stitched_cube[extnames[k]].data = single_comp[extnames[k]].data
                    elif 'COMP2' in extnames[k]:
                        stitched_cube[extnames[k]].data = np.ones((58,58)) * np.nan
                        # this should really be blue or red shape depending on what it is but I'm just going to leave it as (58,58) 
                    else:
                        stitched_cube[extnames[k]].data = single_comp[extnames[k]].data

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    # next few lines are the same as above
                    if 'COMP1' in extnames[k]:
                        stitched_cube[extnames[k]].data = single_comp[extnames[k]].data
                    elif 'COMP2' in extnames[k]:
                        stitched_cube[extnames[k]].data = np.ones((58,58)) * np.nan 
                    else:
                        stitched_cube[extnames[k]].data = single_comp[extnames[k]].data

                else:
                    # loop over each extension and replace nan data with new fit data
                    for k in range(38):
                        stitched_cube[extnames[k]].data = two_comp[extnames[k]].data
                        #shape = stitched_cube[extnames[k]].data.shape
                        #if (shape == blue_shape) or (shape == red_shape) or (shape == line_shape):
                        #    stitched_cube[extnames[k]].data = two_comp[extnames[k]].data
                        #elif shape == (58, 58):
                        #    stitched_cube[extnames[k]].data = pix_hdu[extnames[k]].data


    sys.exit(0)