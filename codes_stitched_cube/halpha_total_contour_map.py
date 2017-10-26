from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

sys.path.append(taffy_dir + 'codes/')
import vel_channel_map as vcm

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in halpha map from the stitched cube
    halpha = fits.open(taffy_extdir + 'stitched_cube_HALPHA.fits')
    halpha_total = halpha[0].data[0]

    # Plotting
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    colorbrewer_cm = vcm.get_colorbrewer_cm()

    # plot contours
    levels = np.array([0, 100, 500, 1000])
    # try smoothing the map to get smoother contours
    # define kernel
    #kernel = Gaussian2DKernel(stddev=0.9)
    #halpha_total = convolve(halpha_total, kernel, boundary='extend')

    c = ax.contour(X, Y, halpha_total, transform=ax.get_transform(wcs_lzifu),\
     levels=levels, cmap=colorbrewer_cm, linewidths=2.0, interpolation='None')
    ax.clabel(c, inline=True, inline_spacing=0, fontsize=5, fmt='%1.1f', lw=4, ls='-')

    plt.show()

    sys.exit(0)