from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

if __name__ == '__main__':

    # Read in stitched cube
    stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

    # read in [SII] maps from stitched cube
    sii6716 = stitched_cube['SII6716'].data[0]
    sii6731 = stitched_cube['SII6731'].data[0]

    sii6716_err = stitched_cube['SII6716_ERR'].data[0]
    sii6731_err = stitched_cube['SII6731_ERR'].data[0]

    # Choose only spaxels with significant detections
    sig_cut = 3
    sii6716_valid_ind = np.where((sii6716/sii6716_err) > sig_cut)
    sii6731_valid_ind = np.where((sii6731/sii6731_err) > sig_cut)

    zipped_sii6716_idx = zip(sii6716_valid_ind[0], sii6716_valid_ind[1])
    zipped_sii6731_idx = zip(sii6731_valid_ind[0], sii6731_valid_ind[1])

    # create empty array and then populate only the significant spaxels
    sii6716_arr = np.ones(sii6716.shape) * -9999.0
    sii6731_arr = np.ones(sii6731.shape) * -9999.0

    # using the same process as in bpt_plots.py to find the intersection 
    # of 2D arrays. There isn't a numpy built-in func to do this.
    for i in range(len(zipped_sii6716_idx)):

        current_idx = zipped_sii6716_idx[i]

    	if current_idx in zipped_sii6731_idx:
    		sii6716_arr[current_idx] = sii6716[current_idx]
    		sii6731_arr[current_idx] = sii6731[current_idx]

    # NaN out -9999.0
    inval_idx1 = np.where(sii6716_arr == -9999.0)
    inval_idx2 = np.where(sii6731_arr == -9999.0)

    sii6716_arr[inval_idx1] = np.nan
    sii6731_arr[inval_idx2] = np.nan

    #plot ratio 6716/6731
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cax = ax.imshow(sii6716_arr/sii6731_arr, origin='lower', vmin=0, vmax=1.8)

    ax.minorticks_on()
    fig.colorbar(cax)

    plt.show()

    # Close HDU
    stitched_cube.close()

    sys.exit(0)