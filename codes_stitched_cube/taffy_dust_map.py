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
    """
    TO-DO:
    1. Dust correction on spaxel-by-spaxel basis
    """

    # read in halpha and hbeta maps from stitched cube
    halpha = fits.open(taffy_extdir + 'stitched_cube_HALPHA.fits')
    hbeta = fits.open(taffy_extdir + 'stitched_cube_HBETA.fits')
    halpha_err = fits.open(taffy_extdir + 'stitched_cube_HALPHA_ERR.fits')
    hbeta_err = fits.open(taffy_extdir + 'stitched_cube_HBETA_ERR.fits')

    # assign to 2D numpy arrays
    halpha_total = halpha[0].data[0]
    hbeta_total = hbeta[0].data[0]

    halpha_total_err = halpha_err[0].data[0]
    hbeta_total_err = hbeta_err[0].data[0]

    # make snr maps
    halpha_snr = halpha_total / halpha_total_err
    hbeta_snr = hbeta_total / hbeta_total_err

    # find all spaxels iwth significant measurements for both lines
    val_idx1 = np.where(halpha_snr >= 3.0)
    val_idx2 = np.where(hbeta_snr >= 3.0)

    halpha_snrmask = np.ones(halpha_total.shape, dtype=bool)
    hbeta_snrmask = np.ones(hbeta_total.shape, dtype=bool)
    halpha_snrmask[val_idx1[0], val_idx1[1]] = False
    hbeta_snrmask[val_idx2[0], val_idx2[1]] = False
    val_idx = ma.mask_or(halpha_snrmask, hbeta_snrmask)

    # Now find the ratio and get dust extinction
    # this assumes that all the halpha and hbeta comes from
    # case B which is not the case. We are explicitly saying 
    # htat some of the ionized gas emission comes from shocks.
    # you will have to break down the halpha emission into
    # the fractions coming from shocks and from SF separately.

    # loop over all pixels and get a ebv value for each

    # ebv[i,j] = 1.97 * np.log10(halpha[i,j]/hbeta[i,j] / 2.86)
    # this formula with the constant 1.97 out front came from Phil
    # I'm not sure of that constant. So I'm usign the formula I derived for now.

    ebv_map = np.zeros((58,58))
    halpha_frac_from_sf = 0.5  # assumed 0.5 for now
    # You need to do this using the [NII] BPT

    for i in range(58):
        for j in range(58):

            ebv_map[i,j] = 2.37 * np.log10(halpha_total[i,j]/hbeta_total[i,j] / 2.85)

    # apply snr mask
    ebv_map = ma.array(ebv_map, mask=val_idx)

    # figure
    print np.nanmin(4.05*ebv_map)
    print np.nanmax(4.05*ebv_map)

    av_map = 3.1 * ebv_map

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

    #plt.imshow(av_map_bridge, vmin=0, vmax=4, origin='lower', interpolation='None')
    #plt.show()

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

    fig.savefig(ipac_taffy_figdir + 'ebv_map.png', dpi=300, bbox_inches='tight')

    # Close all open HDUs
    halpha.close()
    hbeta.close()
    halpha_err.close()
    hbeta_err.close()

    sys.exit(0)