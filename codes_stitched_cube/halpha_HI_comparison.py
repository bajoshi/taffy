from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties

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

    # Read in HI data from Condon et al. 1993
    hi_cube_hdu = fits.open(home + '/Dropbox/Taffy/HI_Condon_et_al/taffy_HIcube.fits')
    hi_cube = hi_cube_hdu[0].data

    # Now read in optical IFU data
    # Only need red channel for Halpha
    obs_r_hdu = fits.open(taffy_data + 'Taffy_R.fits')
    obs_r = obs_r_hdu[0].data

    # Using my own line fits instead of lzifu
    r_line_total = np.load(taffy_extdir + 'halpha_profile_mylinefits.npy')

    # create wavelength array
    # I read these data from the corresponding headers
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    print r_line_total.shape
    print obs_r.shape
    print red_wav_arr.shape
    sys.exit(0)

    return None

if __name__ == '__main__':
    main()
    sys.exit(0)