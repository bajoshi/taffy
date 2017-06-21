from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

if __name__ == '__main__':
    
    # read in lzifu output file
    h = fits.open(taffy_extdir + 'big_cube_2_comp.fits')

    halpha_total = h['HALPHA'].data[0]
    hbeta_total = h['HBETA'].data[0]

    extinc_map = hbeta_total / halpha_total

    plt.imshow(extinc_map, vmin=0, vmax=1, origin='lower')
    plt.colorbar()
    plt.show()

    sys.exit(0)