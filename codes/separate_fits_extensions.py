from __future__ import division

import numpy as np
from astropy.io import fits

import sys
import os

home = os.getenv('HOME')  # does not have a trailing slash
stacking_analysis_dir = home + "/Desktop/FIGS/stacking-analysis-pears/"

taffy_products = '/Volumes/Bhavins_backup/phil/TAFFY/products/'
ipac_taffy_dir = home + '/Desktop/ipac/taffy/'

sys.path.append(stacking_analysis_dir + 'codes/')
import fast_chi2_jackknife as fcj

def print_extnames(hdu, total_ext):

    for i in range(total_ext):
        extname = hdu[i+1].header['EXTNAME']
        print extname

    return None

if __name__ == '__main__':
    
    # read in taffy lzifu product
    # make sure this is the correct one!!
    filepath = taffy_products + 'TAFFY_2_comp.fits'
    hdulist = fits.open(filepath)
    filename = os.path.basename(filepath)
    filename_noext = filename.split('.')[0]

    total_ext = fcj.get_total_extensions(hdulist)

    for i in range(total_ext-1):  # the -1 is there because i don't need the last extension

        extname = hdulist[i+1].header['EXTNAME']
        print i, extname
        ext_filename = filename_noext + '_' + extname + '.fits'

        new_hdulist = fits.HDUList()
        new_hdulist.append(fits.ImageHDU(data=hdulist[i+1].data, header=hdulist[i+1].header))
        new_hdulist.writeto(taffy_products + ext_filename, clobber=True)

    hdulist.close()
    sys.exit(0)