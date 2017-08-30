from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = home + '/Desktop/ipac/taffy_lzifu/products_work/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

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
    one_comp = fits.open(taffy_products + 'Taffy_1_comp_patched.fits')
    # put red line fit in array
    r_line = one_comp['R_LINE_COMP1'].data

    # read in linefits
    mapname = 'vel'
    comp = 1
    if mapname == 'vel':
        map_comp1 = np.load(savedir + 'vel_halpha_comp1.npy')
        map_comp2 = np.load(savedir + 'vel_halpha_comp2.npy')
        map_onecomp = np.load(savedir + 'vel_halpha_onecomp.npy')

    elif mapname == 'vdisp':
        map_comp1 = np.load(savedir + 'std_halpha_comp1.npy')
        map_comp2 = np.load(savedir + 'std_halpha_comp2.npy')
        map_onecomp = np.load(savedir + 'std_halpha_onecomp.npy')

    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # also define line wav
    halpha_air_wav = 6562.8

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

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
                    if comp == 1:
                        map_cube[i,j] = map_comp1[i,j]
                    elif comp == 2:
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
    hdr['EXTNAME'] = mapname + '_comp' + str(comp)
    hdu = fits.PrimaryHDU(data=map_cube, header=hdr)
    hdu.writeto(savedir + mapname + '_cube_comp' + str(comp) + '.fits', clobber=True)

    # close fits files
    one_comp.close()
    h.close()

    sys.exit(0)