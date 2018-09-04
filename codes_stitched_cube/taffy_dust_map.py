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

    ebv_map = np.zeros((58,58))
    halpha_frac_from_sf = 0.62
    # This fraction will probably be different for 
    # the three regions: north, south, and bridge

    for i in range(58):
        for j in range(58):
            ebv_map[i,j] = 1.97 * np.log10(halpha_total[i,j]/hbeta_total[i,j] / 2.86)

    # apply snr mask
    ebv_map = ma.array(ebv_map, mask=val_idx)

    # figure
    print np.nanmin(4.05*ebv_map)
    print np.nanmax(4.05*ebv_map)

    av_map = 4.05 * ebv_map
    a_halpha_map = 3.33 * ebv_map

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

    print "\n", "Mean visual extinction in north galaxy:", np.nanmean(av_map_north)
    print "Mean visual extinction in south galaxy:", np.nanmean(av_map_south)
    print "Mean visual extinction in bridge:", np.nanmean(av_map_bridge)

    # ------------- Extinction at Halpha -------------- #
    a_halpha_map_bridge = ma.array(a_halpha_map, mask=bridge_mask)
    a_halpha_map_north = ma.array(a_halpha_map, mask=north_mask)
    a_halpha_map_south = ma.array(a_halpha_map, mask=south_mask)

    print "\n", "Mean Halhpa extinction in north galaxy:", np.nanmean(a_halpha_map_north)
    print "Mean Halhpa extinction in south galaxy:", np.nanmean(a_halpha_map_south)
    print "Mean Halhpa extinction in bridge:", np.nanmean(a_halpha_map_bridge)

    # Intrinsic Halpha fluxes using results from other code (halpha_frac_from_sf.py)
    aha_n = np.nanmean(a_halpha_map_north)
    aha_s = np.nanmean(a_halpha_map_south)
    aha_br = np.nanmean(a_halpha_map_bridge)

    # Luminostiy from other code
    lha_low_obs_n = 2.829e40
    lha_low_obs_s = 5.895e40
    lha_low_obs_br = 1.92e40

    lha_high_obs_n = 3.254e40
    lha_high_obs_s = 2.419e40
    lha_high_obs_br = 2.113e40

    lum_ha_low_int = lha_low_obs_n * 10**(0.4*aha_n) + lha_low_obs_s * 10**(0.4*aha_s) + lha_low_obs_br * 10**(0.4*aha_br)
    lum_ha_high_int = lha_high_obs_n * 10**(0.4*aha_n) + lha_high_obs_s * 10**(0.4*aha_s) + lha_high_obs_br * 10**(0.4*aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity in the low vel comp:", "{:.3e}".format(lum_ha_low_int)
    print "Intrinsic H-alpha luminosity in the high vel comp:", "{:.3e}".format(lum_ha_high_int)
    lum_ha_from_sf_int = 0.62*lum_ha_low_int + 0.48*lum_ha_high_int
    print "Intrinsic H-alpha luminosity from star-formation:", "{:.3e}".format(lum_ha_from_sf_int)

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

    cax = ax.imshow(av_map, vmin=0, vmax=3, origin='lower', interpolation='None')
    fig.colorbar(cax)
    ax.minorticks_on()

    fig.savefig(ipac_taffy_figdir + 'av_map.png', dpi=300, bbox_inches='tight')

    # Save as FITS file
    # Fix mask
    # this line sets all masked entries to NaN because 
    # the fits file cannot handle masked arrays yet
    av_map = ma.filled(av_map, np.nan)

    # Get original header so that we can use the same WCS
    # I'm simply copying the header from the stitched file and pasting it here
    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()
    hdr = lzifu_hdulist['B_LINE'].header

    # Create file
    hdu = fits.PrimaryHDU(av_map, header=hdr)
    # write
    hdu.writeto(taffy_extdir + 'av_map.fits', overwrite=True)

    # Close all open HDUs
    halpha.close()
    hbeta.close()
    halpha_err.close()
    hbeta_err.close()

    sys.exit(0)