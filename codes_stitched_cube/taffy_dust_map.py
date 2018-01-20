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

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # read in halpha and hbeta maps from stitched cube
    halpha = fits.open(taffy_extdir + 'stitched_cube_HALPHA.fits')
    hbeta = fits.open(taffy_extdir + 'stitched_cube_HBETA.fits')
    halpha_total = halpha[0].data[0]
    hbeta_total = hbeta[0].data[0]

    # this assumes that all the halpha and hbeta comes from
    # case B which is not the case. We are explicitly saying 
    # htat some of the ionized gas emission comes from shocks.
    # you will have to break down the halpha emission into
    # the fractions coming from shocks and from SF separately.

    # loop over all pixels and get a ebv value for each
    # ebv[i,j] = 1.97 * np.log10(halpha[i,j]/hbeta[i,j] / 2.86)
    ebv_map = np.zeros((58,58))

    for i in range(58):
        for j in range(58):

            ebv_map[i,j] = 1.97 * np.log10(halpha_total[i,j]/hbeta_total[i,j] / 2.86)

    # apply mask
    ebv_map = ma.array(ebv_map, mask=all_mask)

    # figure
    print np.nanmin(ebv_map)
    print np.nanmax(ebv_map)

    # read in i band SDSS image
    #sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    #lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    #fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.imshow(ebv_map, vmin=0, vmax=2, origin='lower')
    #plt.colorbar(ax)

    fig.savefig(ipac_taffy_figdir + 'ebv_map.eps', dpi=150, bbox_inches='tight')

    plt.show()

    sys.exit(0)