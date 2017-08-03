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
    mapname = 'vel'
    comp = '2'
    if mapname == 'vel':
        vel_comp1_hdu = fits.open(savedir + 'vel_halpha_comp1.fits')
        vel_comp2_hdu = fits.open(savedir + 'vel_halpha_comp2.fits')
        vel_onecomp_hdu = fits.open(savedir + 'vel_halpha_onecomp.fits')

        map_comp1 = vel_comp1_hdu[0].data
        map_comp2 = vel_comp2_hdu[0].data
        map_onecomp = vel_onecomp_hdu[0].data

    elif mapname == 'vdisp':
        vdisp_comp1_hdu = fits.open(savedir + 'std_halpha_comp1.fits')
        vdisp_comp2_hdu = fits.open(savedir + 'std_halpha_comp2.fits')
        vdisp_onecomp_hdu = fits.open(savedir + 'std_halpha_onecomp.fits')

        map_comp1 = vdisp_comp1_hdu[0].data
        map_comp2 = vdisp_comp2_hdu[0].data
        map_onecomp = vdisp_onecomp_hdu[0].data

    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # also define line wav
    halpha_air_wav = 6562.8

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels')

    # loop over all spaxels and check case and save in array
    map_cube = np.ones((58,58)) * -9999.0

    for i in range(58):
        for j in range(58):

            if single_idx[i,j]:
                map_cube[i,j] = map_onecomp[i,j]

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] or comp2_inv_idx[i,j]:
                    map_cube[i,j] = map_onecomp[i,j]

                else:
                    if comp == '1':
                        map_cube[i,j] = map_comp1[i,j]
                    elif comp == '2':
                        map_cube[i,j] = map_comp2[i,j]

            # convert velocities to physical units
            if mapname == 'vel':
                try:
                    wavidx = int(map_cube[i,j])
                    # heliocentric
                    map_cube[i,j] = ((red_wav_arr[wavidx] - halpha_air_wav) / halpha_air_wav) * speed_of_light
                    # relative to systemic
                    map_cube[i,j] -= speed_of_light * redshift

                except IndexError as e:
                    map_cube[i,j] = np.nan

            # convert vel dispersions to physical units
            if mapname == 'vdisp':
                try:
                    map_cube[i,j] = ((map_cube[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                except IndexError as e:
                    map_cube[i,j] = np.nan

    # apply mask
    # first nan all spaxels outside of region of interest
    all_mask = np.logical_not(all_mask)
    map_cube = np.multiply(map_cube, all_mask)
    mask_idx = np.where(map_cube == 0)
    map_cube[mask_idx] = np.nan

    if mapname == 'vel':
        # now nan all spaxels that are not within an acceptable physical range
        min_idx = np.where(map_cube < -400)
        max_idx = np.where(map_cube > 400)

        map_cube[min_idx] = np.nan
        map_cube[max_idx] = np.nan

        #plt.imshow(map_cube, origin='lower', vmin=-250, vmax=250)
        #plt.colorbar()
        #plt.show()
        #sys.exit(0)

    if mapname == 'vdisp':
        # now nan all spaxels that are not within an acceptable physical range
        min_idx = np.where(map_cube < 0)
        max_idx = np.where(map_cube > 350)

        map_cube[min_idx] = np.nan
        map_cube[max_idx] = np.nan

    # write out map
    # get header from lzifu output
    hdr = one_comp['CHI2'].header
    hdr['EXTNAME'] = mapname + '_comp' + comp
    hdu = fits.PrimaryHDU(data=map_cube, header=hdr)
    hdu.writeto(savedir + mapname + '_cube_comp' + comp + '.fits', clobber=True)

    sys.exit(0)