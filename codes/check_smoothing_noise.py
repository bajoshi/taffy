from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'
taffy_data = taffy_extdir + 'data/'

if __name__ == '__main__':

    # -------------------------------------------------- Read in data -------------------------------------------------- #
    # read in observed data
    obs_b = fits.open(taffy_data + 'Taffy_unsmoothed_B.fits')
    obs_r = fits.open(taffy_data + 'Taffy_unsmoothed_R.fits')

    # read in smoothed data
    smoothed_b = fits.open(taffy_data + 'Taffy_B.fits')
    smoothed_r = fits.open(taffy_data + 'Taffy_R.fits')

    # -------------------------------------------------- Create arrays for plotting -------------------------------------------------- #
    # create lambda arrays
    b_cont_hdr = obs_b[0].header
    tot_spec_elem = b_cont_hdr['NAXIS3']
    spec_start = b_cont_hdr['CRVAL3']
    spec_delta = b_cont_hdr['CDELT3']
    lam_end = spec_start + tot_spec_elem*spec_delta - spec_delta
    lam_b = np.linspace(spec_start, lam_end, tot_spec_elem)

    r_cont_hdr = obs_r[0].header
    tot_spec_elem = r_cont_hdr['NAXIS3']
    spec_start = r_cont_hdr['CRVAL3']
    spec_delta = r_cont_hdr['CDELT3']
    lam_end = spec_start + tot_spec_elem*spec_delta - spec_delta
    lam_r = np.linspace(spec_start, lam_end, tot_spec_elem)

    # -------------------------------------------------- Plotting -------------------------------------------------- #
    # Plot some spaxels from both cubes side by side
    fig = plt.figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    # Choose spaxel here
    pix_x = 23
    pix_y = 46
    arr_x = pix_y - 1
    arr_y = pix_x - 1

    print 'Current pixel (ds9 coords)', pix_x, pix_y

    ax1.plot(lam_b, obs_b[0].data[:,arr_x,arr_y])
    ax2.plot(lam_r, obs_r[0].data[:,arr_x,arr_y])

    ax3.plot(lam_b, smoothed_b[0].data[:,arr_x,arr_y])
    ax4.plot(lam_r, smoothed_r[0].data[:,arr_x,arr_y])

    # Text to indicate spaxel and channel

    plt.show()

    # Close all fits files and exit
    obs_b.close()
    obs_r.close()
    smoothed_b.close()
    smoothed_r.close()

    sys.exit(0)