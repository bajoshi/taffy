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
halpha_air_wav = 6562.8

def one_comp_flux(i, j, amp_onecomp, vdisp_onecomp):

    amp = amp_onecomp[i,j]
    vdisp = vdisp_onecomp[i,j]
    vdisp = ((vdisp * 0.3) / halpha_air_wav) * speed_of_light
    intg_flux = amp * np.sqrt(2 * np.pi * vdisp**2)

    return intg_flux

def get_interp_val_spaxel(maparr, arr_x, arr_y):

    if (arr_x > 0) and (arr_y > 0):
        ll_x = arr_x - 1
        ll_y = arr_y - 1

    elif (arr_x == 0) or (arr_y == 0):

        if arr_x == 0:
            ll_x = 0
            ll_y = arr_y - 1

        if arr_y == 0:
            ll_x = arr_x - 1
            ll_y = 0

        if (arr_x == 0) and (arr_y == 0):
            ll_x = 0
            ll_y = 0

    mean_arr = []

    for v in range(ll_x, ll_x+3):
        for w in range(ll_y, ll_y+3):
            try:
                if np.isfinite(maparr[v,w]):
                    mean_arr.append(maparr[v,w])
            except IndexError as e:
                continue

    mean_arr = np.asarray(mean_arr)
    nan_idx = np.where(mean_arr == -9999.0)
    mean_arr[nan_idx] = np.nan

    min_idx = np.where(mean_arr < 0)
    mean_arr[min_idx] = np.nan

    max_idx = np.where(mean_arr > 1e6)
    mean_arr[max_idx] = np.nan

    mean_val = np.nanmedian(mean_arr)

    return mean_val

def get_interp_spaxel_list(comp):

    if comp == 1:
        interp_spaxel_list = [[27,43],[44,12]]
    elif comp == 2:
        interp_spaxel_list = [[51,26]]

    interp_arr_x = []
    interp_arr_y = []

    for u in range(len(interp_spaxel_list)):

        interp_arr_x.append(interp_spaxel_list[u][1] - 1)
        interp_arr_y.append(interp_spaxel_list[u][0] - 1)

    interp_spaxel_arr = zip(interp_arr_x, interp_arr_y)

    return interp_spaxel_arr

if __name__ == '__main__':
    
    # read in lzifu result for single comp fit
    one_comp = fits.open(taffy_products + 'Taffy_1_comp_patched.fits')

    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    # read in linefits
    amp_comp1 = np.load(savedir + 'amp_halpha_comp1.npy')
    amp_comp2 = np.load(savedir + 'amp_halpha_comp2.npy')
    amp_onecomp = np.load(savedir + 'amp_halpha_onecomp.npy')

    vdisp_comp1 = np.load(savedir + 'std_halpha_comp1.npy')
    vdisp_comp2 = np.load(savedir + 'std_halpha_comp2.npy')
    vdisp_onecomp = np.load(savedir + 'std_halpha_onecomp.npy')

    # loop over all spaxels and check case and save in array
    intg_flux_map = np.ones((58,58)) * -9999.0
    comp = 2

    nan_single_comp_arr = sm.get_nan_arr()
    interp_spaxel_arr = get_interp_spaxel_list(comp)

    # if you want to check a pixel
    pix_x = 35
    pix_y = 30
    arr_x = pix_y - 1
    arr_y = pix_x - 1

    for i in range(58):
        for j in range(58):

            if (i,j) in nan_single_comp_arr:
                intg_flux_map[i,j] = one_comp_flux(i, j, amp_onecomp, vdisp_onecomp)
                continue

            if single_idx[i,j]:
                intg_flux_map[i,j] = one_comp_flux(i, j, amp_onecomp, vdisp_onecomp)

            elif diffmean_idx[i,j] or diffstd_idx[i,j] or diffboth_idx[i,j]:
                if comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    intg_flux_map[i,j] = one_comp_flux(i, j, amp_onecomp, vdisp_onecomp)

                elif comp1_inv_idx[i,j] and not comp2_inv_idx[i,j]:
                    vdisp_comp2[i,j] = ((vdisp_comp2[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                    intg_flux_map[i,j] = amp_comp2[i,j] * np.sqrt(2 * np.pi * vdisp_comp2[i,j]**2)

                elif not comp1_inv_idx[i,j] and comp2_inv_idx[i,j]:
                    vdisp_comp1[i,j] = ((vdisp_comp1[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                    intg_flux_map[i,j] = amp_comp1[i,j] * np.sqrt(2 * np.pi * vdisp_comp1[i,j]**2)

                else:
                    if comp == 1:
                        vdisp_comp1[i,j] = ((vdisp_comp1[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                        intg_flux_map[i,j] = amp_comp1[i,j] * np.sqrt(2 * np.pi * vdisp_comp1[i,j]**2)
                    elif comp == 2:
                        vdisp_comp2[i,j] = ((vdisp_comp2[i,j] * 0.3) / halpha_air_wav) * speed_of_light
                        intg_flux_map[i,j] = amp_comp2[i,j] * np.sqrt(2 * np.pi * vdisp_comp2[i,j]**2)

    # the interpolation of spaxels has to be done AFTER the map has been filled once
    # because if there aren't any values around the spaxel then nothing gets filled in
    # e.g. think about interpolating 0,0 without actually having any values around it.
    for w in range(len(interp_spaxel_arr)):
        i = interp_spaxel_arr[w][0]
        j = interp_spaxel_arr[w][1]
        intg_flux_map[i,j] = get_interp_val_spaxel(intg_flux_map, i, j)

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # apply mask
    # first nan all spaxels outside of region of interest
    all_mask = np.logical_not(all_mask)
    intg_flux_map = np.multiply(intg_flux_map, all_mask)
    mask_idx = np.where(intg_flux_map == 0)
    intg_flux_map[mask_idx] = np.nan

    # now nan all spaxels that are not within an acceptable physical range
    min_idx = np.where(intg_flux_map < 0)
    max_idx = np.where(intg_flux_map > 1e6)

    intg_flux_map[min_idx] = np.nan
    intg_flux_map[max_idx] = np.nan

    # write it out
    hdr = one_comp['CHI2'].header
    hdr['EXTNAME'] = 'INTG_FLUX_COMP' + str(comp)
    hdu = fits.PrimaryHDU(data=intg_flux_map, header=hdr)
    hdu.writeto(savedir + 'intg_flux_cube_comp' + str(comp) + '.fits', clobber=True)

    # close all opened fits files
    one_comp.close()
    h.close()

    sys.exit(0)