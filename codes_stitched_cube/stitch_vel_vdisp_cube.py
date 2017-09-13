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
import stitch_map as sm

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
    comp = 2
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

    # get nana spaxel array from stich_map code
    nan_single_comp_arr = sm.get_nan_arr()

    # loop over all spaxels and check case and save in array
    map_cube = np.ones((58,58)) * -9999.0

    for i in range(58):
        for j in range(58):

            if all_mask[i,j]:
                continue

            if (i,j) in nan_single_comp_arr:
                map_cube[i,j] = map_onecomp[i,j]
                continue

            if single_idx[i,j]:
                map_cube[i,j] = map_onecomp[i,j]

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    map_cube[i,j] = map_onecomp[i,j]

                elif comp1_inv_idx[i,j] and not comp2_inv_idx[i,j]:
                    map_cube[i,j] = map_comp2[i,j]

                elif not comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    map_cube[i,j] = map_comp1[i,j]

                else:
                    if comp == 1:
                        map_cube[i,j] = map_comp1[i,j]
                    elif comp == 2:
                        map_cube[i,j] = map_comp2[i,j]

    # convert velocities and velocity dispersions to physical units
    # this has to be done separately because of the continue statement
    # in the previous for loop that causes it to skip converting units 
    # for pixels that were in the nan_single_comp_arr. This means that
    # the values in those pixels are not within the physical range and
    # later on get replaced by np.nan
    for u in range(58):
        for v in range(58):
            if mapname == 'vel':
                try:
                    current_wavidx = int(map_cube[u,v])
                    # heliocentric
                    map_cube[u,v] = ((red_wav_arr[current_wavidx] - halpha_air_wav) / halpha_air_wav) * speed_of_light
                    # relative to systemic
                    map_cube[u,v] -= speed_of_light * redshift

                except IndexError as e:
                    map_cube[u,v] = np.nan

            if mapname == 'vdisp':
                try:
                    map_cube[u,v] = ((map_cube[u,v] * 0.3) / halpha_air_wav) * speed_of_light
                except IndexError as e:
                    map_cube[u,v] = np.nan

    # apply mask
    # first nan all spaxels that still have -9999.0
    mask_idx = np.where(map_cube == -9999.0) 
    map_cube[mask_idx] = np.nan

    #if mapname == 'vel':
    #    # now nan all spaxels that are not within an acceptable physical range
    #    min_idx = np.where(map_cube < -400)
    #    max_idx = np.where(map_cube > 400)

    #    map_cube[min_idx] = np.nan
    #    map_cube[max_idx] = np.nan

    #    #plt.imshow(map_cube, origin='lower', vmin=-250, vmax=250)
    #    #plt.colorbar()
    #    #plt.show()
    #    #sys.exit(0)

    if mapname == 'vdisp':
        # now nan all spaxels that are not within an acceptable physical range
        min_idx = np.where(map_cube < 0)
        #max_idx = np.where(map_cube > 2000)

        map_cube[min_idx] = np.nan
        #map_cube[max_idx] = np.nan

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