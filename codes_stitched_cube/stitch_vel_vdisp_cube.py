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

if __name__ == '__main__':
    
    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    # read in velocity and vdisp maps for each component
    # and also read in lzifu result for single comp fit
    one_comp = fits.open(taffy_extdir + 'products/Taffy_1_comp.fits')
    # put red line fit in array
    r_line = one_comp['R_LINE_COMP1'].data

    # read in linefits
    mapname = 'vdisp'
    if mapname == 'vel':
        vel_comp1_hdu = fits.open(savedir + 'vel_halpha_comp1.fits')
        vel_comp2_hdu = fits.open(savedir + 'vel_halpha_comp2.fits')

        map_comp1 = vel_comp1_hdu[0].data
        map_comp2 = vel_comp2_hdu[0].data
    elif mapname == 'vdisp':
        vdisp_comp1_hdu = fits.open(savedir + 'std_halpha_comp1.fits')
        vdisp_comp2_hdu = fits.open(savedir + 'std_halpha_comp2.fits')

        map_comp1 = vdisp_comp1_hdu[0].data
        map_comp2 = vdisp_comp2_hdu[0].data

    # loop over all spaxels and check case and save in array
    map_cube = np.ones((58,58)) * -9999.0
    comp = '2'

    for i in range(58):
        for j in range(58):

            if single_idx[i,j]:
                map_cube[i,j] = np.mean([map_comp1[i,j], map_comp2[i,j]])

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    map_cube[i,j] = np.mean([map_comp1[i,j], map_comp2[i,j]])

                else:
                    if comp == '1':
                        map_cube[i,j] = map_comp1[i,j]
                    elif comp == '2':
                        map_cube[i,j] = map_comp2[i,j]

    hdu = fits.PrimaryHDU(data=map_cube)
    hdu.writeto(savedir + mapname + '_cube_comp' + comp + '.fits', clobber=True)

    sys.exit(0)