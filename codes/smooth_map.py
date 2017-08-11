from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel

import sys

if __name__ == '__main__':
    
    # define kernel
    kernel = Gaussian2DKernel(stddev=1.47)

    #read in taffy data
    #taffy_datadir = '/Users/baj/Desktop/ipac/taffy_lzifu/data/'
    taffy_datadir = '/Volumes/Bhavins_backup/ipac/TAFFY/data/'
    taffy_b = fits.open(taffy_datadir + 'Taffy_unsmoothed_B.fits')
    taffy_r = fits.open(taffy_datadir + 'Taffy_unsmoothed_R.fits')

    taffy_cubes = [taffy_b, taffy_r]

    # loop over each wavelength slice in both cubes and smooth it
    for taffy_map in taffy_cubes:

        for i in range(len(taffy_map[0].data)):

            taffy_map[0].data[i] = convolve(taffy_map[0].data[i], kernel, boundary='extend')
    
        # currently assuming that the noise decreases because of the smoothing by a blanket factor of 2
        taffy_map[1].data /= 2

    taffy_b.writeto(taffy_datadir + 'Taffy_B.fits', output_verify='ignore', clobber=True)
    taffy_r.writeto(taffy_datadir + 'Taffy_R.fits', output_verify='ignore', clobber=True)

    sys.exit(0)