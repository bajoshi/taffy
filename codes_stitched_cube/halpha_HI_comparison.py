from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
import Polygon as pg

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

def extract_spec(reg_pn, cube):
    """
    This function will extract and return the spectrum 
    corresponding to the given region from the given cube.
    Arguments:
    --- reg: ds9 region polygon in image coordinates.
    i.e., convert the region from the ds9 file to a polygon object.

    --- cube: data cube from which spectrum is to be extracted.
    """

    # 

    return spec_wav, spec_flux

def main():

    # ------------------------- Read in data ------------------------- #
    # Read in HI data from Condon et al. 1993
    hi_cube_hdu = fits.open(home + '/Dropbox/Taffy/HI_Condon_et_al/taffy_HIcube.fits')
    hi_cube = hi_cube_hdu[0].data

    # Now read in optical IFU data
    # Only need red channel for Halpha
    obs_r_hdu = fits.open(taffy_data + 'Taffy_R.fits')
    obs_r = obs_r_hdu[0].data

    # Using my own line fits instead of lzifu
    r_line_total = np.load(taffy_extdir + 'halpha_profile_mylinefits.npy')

    # ------------------------- Other preliminaries ------------------------- #
    # create wavelength array
    # I read these data from the corresponding headers
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # Make sure that rest wavelength is Halpha wavelength in air!! the IFU data was taken on the ground
    # find line index in wavelength array
    redshift = 0.0145  # average z 
    sys_vel = redshift * speed_of_light  # avg systemic velocity is z*c = 4350 km/s
    halpha_air_wav = 6562.80
    halpha_wav = halpha_air_wav*(1+redshift)
    halpha_idx = np.argmin(abs(red_wav_arr - halpha_wav))

    # ------------------------- Plotting ------------------------- #
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss) 

    # Now read in regions file and convert to polygons
    regfile = open()
    extract_spec(reg_pn, cube)

    plt.show()

    return None

if __name__ == '__main__':
    main()
    sys.exit(0)