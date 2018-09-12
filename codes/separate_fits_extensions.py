from __future__ import division

import numpy as np
from astropy.io import fits

import sys
import os

home = os.getenv('HOME')  # does not have a trailing slash
stacking_analysis_dir = home + "/Desktop/FIGS/stacking-analysis-pears/"

taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_dir = home + '/Desktop/ipac/taffy/'

sys.path.append(stacking_analysis_dir + 'codes/')
import fast_chi2_jackknife as fcj

def get_extnames(hdu, total_ext):

    extname_list = []

    for i in range(total_ext):
        extname = hdu[i+1].header['EXTNAME']
        extname_list.append(extname)

    return extname_list

if __name__ == '__main__':
    
    # read in taffy lzifu product
    # make sure this is the correct one!!
    filepath = taffy_extdir + 'stitched_cube.fits'
    hdulist = fits.open(filepath)
    filename = os.path.basename(filepath)
    filename_noext = filename.split('.')[0]

    basedir = os.path.dirname(filepath) + '/'
    print "Separating extensions from file:", filename_noext, "in directory:", basedir

    total_ext = fcj.get_total_extensions(hdulist)

    # Put in a list here to extract only the specific extensions 
    # which are in the list or change the get_specific_ext 
    # variable to 'all' which will extract all extensions.
    get_specific_ext = 'all'  #['V']
    if get_specific_ext == 'all':
        get_specific_ext = get_extnames(hdulist, total_ext)

    for i in range(total_ext-1):  # the -1 is there because i don't need the last extension

        extname = hdulist[i+1].header['EXTNAME']
        
        if extname in get_specific_ext:
            print i, extname
            ext_filename = filename_noext + '_' + extname + '.fits'

            new_hdulist = fits.HDUList()
            new_hdulist.append(fits.ImageHDU(data=hdulist[i+1].data, header=hdulist[i+1].header))
            new_hdulist.writeto(basedir + ext_filename, overwrite=True)
        else:
            continue

    hdulist.close()
    sys.exit(0)