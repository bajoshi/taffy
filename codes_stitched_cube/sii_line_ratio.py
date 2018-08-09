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