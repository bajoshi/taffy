from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from astropy.wcs import WCS
from astropy.visualization import ManualInterval, ZScaleInterval, LogStretch, ImageNormalize

import os
import sys
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText
from matplotlib.colors import LinearSegmentedColormap

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffy_dir + 'codes/')
import bpt_plots as bpt
import vel_channel_map as vcm

if __name__ == '__main__':
	
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # constants:
    lightspeed = 299792.458  # km/s

    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    h, wcs_lzifu = vcm.get_lzifu_products()
    del h 
    h = fits.open(savedir + 'stitched_cube.fits')

    # put fit in arr
    r_line_total = 

	sys.exit(0)