from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products_big_cube_velsort/'
taffy_data = '/Volumes/Bhavins_backup/ipac/TAFFY/data/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

if __name__ == '__main__':
    
    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    print h.info()
    sys.exit(0)

    single_idx = h[1].data
    diffmean_idx = h[2].data
    diffstd_idx = h[3].data
    diffboth_idx = h[4].data

    # read in both 1 and 2 comp fitting results
    single_comp = fits.open(taffy_products + 'Taffy_1_comp.fits')
    two_comp = fits.open(taffy_products + 'big_cube_2_comp_velsort.fits')

    # loop over all pixels and stitch them according 
    # to whether they need a one or two comp fit
    #for i in range(58):
    #    for j in range(58):



    sys.exit(0)