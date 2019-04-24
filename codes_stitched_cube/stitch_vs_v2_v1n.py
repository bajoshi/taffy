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
gaussfits_dir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import stitch_map as sm

speed_of_light = 299792.458  # km/s
redshift = 0.0145

def genmap():

    # ---------------- Read in Inputs ---------------- #
    # read in indices file
    h = fits.open(gaussfits_dir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    # read in velocity and vdisp maps for each component
    # and also read in lzifu result for single comp fit
    one_comp = fits.open(taffy_extdir + 'Taffy_1_comp_patched.fits')
    one_comp_vel = one_comp['V'].data[1]
    # put red line fit in array
    #r_line = one_comp['R_LINE_COMP1'].data

    # Read in velocities from fitting results
    map_comp1 = np.load(gaussfits_dir + 'vel_halpha_comp1.npy')
    map_comp2 = np.load(gaussfits_dir + 'vel_halpha_comp2.npy')
    map_onecomp = np.load(gaussfits_dir + 'vel_halpha_onecomp.npy')

    # ---------------- Prep for stitching ---------------- #
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # also define line wav
    halpha_air_wav = 6562.8  # Angstroms # in air

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # get nan spaxel array from stich_map code
    nan_single_comp_arr = sm.get_nan_arr()

    # Read in both velocity component maps
    # Comp2 map already has Vs + V2 + V2b stitched
    # You only need to replace the V2b spaxels wtih V1n
    # You can simply replace those spaxels from the comp1 map
    map_vel_comp2_hdu = fits.open(gaussfits_dir + 'vel_cube_comp2.fits')
    map_vel_comp1_hdu = fits.open(gaussfits_dir + 'vel_cube_comp1.fits')

    map_vel_comp1 = map_vel_comp1_hdu[0].data
    map_vel_comp2 = map_vel_comp2_hdu[0].data

    # Create empty new map
    map_vel_comp2_new = np.zeros(map_vel_comp2.shape)

    for i in range(58):
        for j in range(58):

            map_vel_comp2_new[i,j] = map_vel_comp2[i,j]

            if diffstd_idx[i,j]:
                map_vel_comp2_new[i,j] = map_vel_comp1[i,j]

    # write out map
    # get header from lzifu output
    hdr = one_comp['CHI2'].header
    hdr['EXTNAME'] = 'vel_comp_vs_v2_v1n'
    hdu = fits.PrimaryHDU(data=map_vel_comp2_new, header=hdr)
    hdu.writeto(gaussfits_dir + 'vel_cube_vs_v2_v1n.fits', overwrite=True)

    return None

def main():
    """
    The purpose of this code is to compare the Vs+V1+V1n map
    (this is currently the top middle panel in Fig 7) with the
    Vs+V2+V1n map. This is one of the plots the referee wanted 
    to see. I think they think that this will display the rotation
    better. I'm not entirely sure.

    This code is completely based on the codes --
    1. stitche_vel_vdisp_cube.py and 2. contours_from_linefits.py.
    """

    # Generate the Vs + V2 + V1n map
    #genmap()

    # Read in the new map and make contours
    map_vel_comp2_hdu = fits.open(gaussfits_dir + 'vel_cube_comp2.fits')
    map_vel_comp2_new_hdu = fits.open(gaussfits_dir + 'vel_cube_vs_v2_v1n.fits')

    map_vel_comp2 = map_vel_comp2_hdu[0].data
    map_vel_comp2_new = map_vel_comp2_new_hdu[0].data

    # Diagnostic figures
    # Do not delete this code block. Useful for checking.
    """
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.imshow(map_vel_comp2, origin='lower', vmin=-350, vmax=350)
    ax2.imshow(map_vel_comp2_new, origin='lower', vmin=-350, vmax=350)

    plt.show()
    """

    # ------- Contours --------- # 


    return None

if __name__ == '__main__':
    main()
    sys.exit(0)