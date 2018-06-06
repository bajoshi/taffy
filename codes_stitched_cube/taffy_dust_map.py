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
    print np.nanmin(4.05*ebv_map)
    print np.nanmax(4.05*ebv_map)

    av_map = 4.05 * ebv_map

    # Get average values in the different regions
    # get the region masks first
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # Combine all masks with checker board mask
    # See checkerboard mask comment in BPT velo comp code
    checkerboard_mask = np.zeros((58,58), dtype=np.int)
    checkerboard_mask[::2, 1::2] = 1
    checkerboard_mask[1::2, ::2] = 1

    # Combine checkerboard and two comp masks and then apply
    # First convert all masks to boolean
    checkerboard_mask = checkerboard_mask.astype(bool)
    bridge_mask = bridge_mask.astype(bool)
    north_mask = north_mask.astype(bool)
    south_mask = south_mask.astype(bool)

    bridge_mask = np.ma.mask_or(checkerboard_mask, bridge_mask)
    north_mask = np.ma.mask_or(checkerboard_mask, north_mask)
    south_mask = np.ma.mask_or(checkerboard_mask, south_mask)

    # Now apply masks
    av_map_bridge = ma.array(av_map, mask=bridge_mask)
    av_map_north = ma.array(av_map, mask=north_mask)
    av_map_south = ma.array(av_map, mask=south_mask)

    print "Mean extinction in north galaxy:", np.nanmean(av_map_north)
    print "Mean extinction in south galaxy:", np.nanmean(av_map_south)
    print "Mean extinction in bridge:", np.nanmean(av_map_bridge)

    plt.imshow(av_map_bridge, vmin=0, vmax=5, origin='lower', interpolation='None')
    plt.show()

    # read in i band SDSS image
    #sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    #lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    #fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cax = ax.imshow(av_map, vmin=0, vmax=5, origin='lower', interpolation='None')
    fig.colorbar(cax)
    ax.minorticks_on()

    fig.savefig(ipac_taffy_figdir + 'ebv_map.png', dpi=150, bbox_inches='tight')

    sys.exit(0)