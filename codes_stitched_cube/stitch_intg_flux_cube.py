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

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

speed_of_light = 299792.458
redshift = 0.0145
halpha_air_wav = 6562.8

if __name__ == '__main__':
    
    # read in lzifu result for single comp fit
    one_comp = fits.open(taffy_extdir + 'products/Taffy_1_comp.fits')

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
    comp = '1'

    for i in range(58):
        for j in range(58):

            if single_idx[i,j]:
                amp = amp_onecomp[i,j]
                vdisp = vdisp_onecomp[i,j]
                vdisp = ((vdisp * 0.3) / halpha_air_wav) * speed_of_light
                intg_flux_map[i,j] = amp * np.sqrt(2 * np.pi * vdisp**2)

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    amp = amp_onecomp[i,j]
                    vdisp = vdisp_onecomp[i,j]
                    vdisp = ((vdisp * 0.3) / halpha_air_wav) * speed_of_light
                    intg_flux_map[i,j] = amp * np.sqrt(2 * np.pi * vdisp**2)

                else:
                    if comp == '1':
                        vdisp_comp1[i,j] = ((vdisp_comp1[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                        intg_flux_map[i,j] = amp_comp1[i,j] * np.sqrt(2 * np.pi * vdisp_comp1[i,j]**2)
                    elif comp == '2':
                        vdisp_comp2[i,j] = ((vdisp_comp2[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                        intg_flux_map[i,j] = amp_comp2[i,j] * np.sqrt(2 * np.pi * vdisp_comp2[i,j]**2)

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels')

    # apply mask
    # first nan all spaxels outside of region of interest
    all_mask = np.logical_not(all_mask)
    intg_flux_map = np.multiply(intg_flux_map, all_mask)
    mask_idx = np.where(intg_flux_map == 0)
    intg_flux_map[mask_idx] = np.nan

    # now nan all spaxels that are not within an acceptable physical range
    min_idx = np.where(intg_flux_map < 0)
    max_idx = np.where(intg_flux_map > 1e5)

    intg_flux_map[min_idx] = np.nan
    intg_flux_map[max_idx] = np.nan

    # write it out
    hdr = one_comp['CHI2'].header
    hdr['EXTNAME'] = 'INTG_FLUX_COMP' + comp
    hdu = fits.PrimaryHDU(data=intg_flux_map, header=hdr)
    hdu.writeto(savedir + 'intg_flux_cube_comp' + comp + '.fits', clobber=True)

    sys.exit(0)