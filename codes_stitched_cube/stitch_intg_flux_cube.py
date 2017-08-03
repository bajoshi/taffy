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

    # read in linefits
    amp_comp1_hdu = fits.open(savedir + 'amp_halpha_comp1.fits')
    amp_comp2_hdu = fits.open(savedir + 'amp_halpha_comp2.fits')
    amp_onecomp_hdu = fits.open(savedir + 'amp_halpha_onecomp.fits')

    amp_comp1 = amp_comp1_hdu[0].data
    amp_comp2 = amp_comp2_hdu[0].data
    amp_onecomp = amp_onecomp_hdu[0].data

    vdisp_comp1_hdu = fits.open(savedir + 'std_halpha_comp1.fits')
    vdisp_comp2_hdu = fits.open(savedir + 'std_halpha_comp2.fits')
    vdisp_onecomp_hdu = fits.open(savedir + 'std_halpha_onecomp.fits')

    vdisp_comp1 = vdisp_comp1_hdu[0].data
    vdisp_comp2 = vdisp_comp2_hdu[0].data
    vdisp_onecomp = vdisp_onecomp_hdu[0].data

    # loop over all spaxels and check case and save in array
    intg_flux_map = np.ones((58,58)) * -9999.0
    comp = '2'

    for i in range(58):
        for j in range(58):

            if single_idx[i,j]:
                amp = amp_onecomp[i,j]
                vdisp = vdisp_onecomp[i,j]
                intg_flux_map[i,j] = amp * np.sqrt(2 * np.pi * vdisp**2)

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    amp = amp_onecomp[i,j]
                    vdisp = vdisp_onecomp[i,j]
                    intg_flux_map[i,j] = amp * np.sqrt(2 * np.pi * vdisp**2)

                else:
                    if comp == '1':
                        intg_flux_map[i,j] = amp_comp1[i,j] * np.sqrt(2 * np.pi * vdisp_comp1[i,j]**2)
                    elif comp == '2':
                        intg_flux_map[i,j] = amp_comp2[i,j] * np.sqrt(2 * np.pi * vdisp_comp2[i,j]**2)

    hdu = fits.PrimaryHDU(data=intg_flux_map)
    hdu.writeto(savedir + 'intg_flux_cube_comp' + comp + '.fits', clobber=True)

    sys.exit(0)