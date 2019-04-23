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

def main():
    """
    The purpose of this code is to compare the Vs+V1+V1n map
    (this is currently the top middle panel in Fig 7) with the
    Vs+V2+V1n map. This is one of the plots the referee wanted 
    to see. I think they think that this will display the rotation
    of Taffy-N better. Not entirely sure.

    This code is completely based on the codes --
    1. stitche_vel_vdisp_cube.py and 2. contours_from_linefits.py.
    """

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

    # Create final map array
    # loop over all spaxels and check case and save in array
    map_cube = np.ones((58,58)) * -9999.0
    # This is set to -9999.0 by default and NOT zeros because 
    # zeros would be a valid measurement which is not what we want.

    # ---------------- Loop over all spaxels and stitch as required ---------------- #
    # In the two nested loops below, I'm using only the results of my OWN fits.
    for i in range(58):
        for j in range(58):

            # If spaxel is outside what I defined earlier as the 
            # "could be NOT NaN" area then skip it. I defined this 
            # area to be outside both galaxies and bridge.
            if all_mask[i,j]:
                continue

            # If spaxel is one of those that are defined as nan in 
            # the above array then it is a spaxel I want to force 
            # to a single compoennt fit.
            if (i,j) in nan_single_comp_arr:
                map_cube[i,j] = map_onecomp[i,j]
                continue

            # Now use the indices defined previously to stitch as required
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


    return None

if __name__ == '__main__':
    main()
    sys.exit(0)