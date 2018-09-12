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

    # ----------------------------------- Read in arrays ----------------------------------- #
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

    # ----------------------------------- Cut on significance ----------------------------------- #
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

    # ----------------------------------- Make E(B-V) map ----------------------------------- #
    # Now find the ratio and get dust extinction
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

    # apply snr mask
    ebv_map = ma.array(ebv_map, mask=val_idx)

    # figure
    print np.nanmin(4.05*ebv_map)
    print np.nanmax(4.05*ebv_map)

    # ----------------------------------- Make Av maps from EBV maps ----------------------------------- #
    # Visual extinction
    av_map = 4.05 * ebv_map

    # Extinction for all emission lines
    line_name_list = ['HBETA', 'OIII5007', 'OI6300', 'HALPHA', 'NII6583', 'SII6716', 'SII6731']
    # These values came from the code get_klambda.py
    klam_Hbeta = 4.598
    klam_OIII_5007 = 4.464
    klam_OI_6300 = 3.49
    klam_Halpha = 3.326
    klam_NII_6583 = 3.313
    klam_SII_6716 = 3.23
    klam_SII_6731 = 3.221

    a_Hbeta_map = klam_Hbeta * ebv_map
    a_OIII_5007_map = klam_OIII_5007 * ebv_map
    a_OI_6300_map = klam_OI_6300 * ebv_map
    a_Halpha_map = klam_Halpha * ebv_map
    a_NII_6583_map = klam_NII_6583 * ebv_map
    a_SII_6716_map = klam_SII_6716 * ebv_map
    a_SII_6731_map = klam_SII_6731 * ebv_map

    # ----------------------------------- Get average extinction values in the different regions and vel comps ----------------------------------- #
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
    a_halpha_map_bridge = ma.array(a_Halpha_map, mask=bridge_mask)
    a_halpha_map_north = ma.array(a_Halpha_map, mask=north_mask)
    a_halpha_map_south = ma.array(a_Halpha_map, mask=south_mask)

    print "\n", "Mean Halhpa extinction in north galaxy:", np.nanmean(a_halpha_map_north)
    print "Mean Halhpa extinction in south galaxy:", np.nanmean(a_halpha_map_south)
    print "Mean Halhpa extinction in bridge:", np.nanmean(a_halpha_map_bridge)

    # Intrinsic Halpha fluxes using results from other code (halpha_frac_from_sf.py)
    aha_n = np.nanmean(a_halpha_map_north)
    aha_s = np.nanmean(a_halpha_map_south)
    aha_br = np.nanmean(a_halpha_map_bridge)

    # Luminostiy from other code
    lha_low_obs_n = 2.823e+40
    lha_low_obs_s = 6.021e+40
    lha_low_obs_br = 1.907e+40

    lha_high_obs_n = 3.264e+40
    lha_high_obs_s = 2.289e+40
    lha_high_obs_br = 2.131e+40

    lum_ha_low_int = lha_low_obs_n * 10**(0.4*aha_n) + lha_low_obs_s * 10**(0.4*aha_s) + lha_low_obs_br * 10**(0.4*aha_br)
    lum_ha_high_int = lha_high_obs_n * 10**(0.4*aha_n) + lha_high_obs_s * 10**(0.4*aha_s) + lha_high_obs_br * 10**(0.4*aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity in the low vel comp:", "{:.3e}".format(lum_ha_low_int)
    print "Intrinsic H-alpha luminosity in the high vel comp:", "{:.3e}".format(lum_ha_high_int)
    lum_ha_from_sf_int = 0.62*lum_ha_low_int + 0.48*lum_ha_high_int
    print "Total intrinsic H-alpha luminosity from star-formation:", "{:.3e}".format(lum_ha_from_sf_int)

    # ------- by component ------- #
    lha_low_int_n = lha_low_obs_n * 10**(0.4 * aha_n)
    lha_low_int_s = lha_low_obs_s * 10**(0.4 * aha_s)
    lha_low_int_br = lha_low_obs_br * 10**(0.4 * aha_br)

    lha_high_int_n = lha_high_obs_n * 10**(0.4 * aha_n)
    lha_high_int_s = lha_high_obs_s * 10**(0.4 * aha_s)
    lha_high_int_br = lha_high_obs_br * 10**(0.4 * aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity in north galaxy in comp1:", lha_low_int_n
    print "Intrinsic H-alpha luminosity in south galaxy in comp1:", lha_low_int_s
    print "Intrinsic H-alpha luminosity in bridge galaxy in comp1:", lha_low_int_br

    print "\n"
    print "Intrinsic H-alpha luminosity in north galaxy in comp2:", lha_high_int_n
    print "Intrinsic H-alpha luminosity in south galaxy in comp2:", lha_high_int_s
    print "Intrinsic H-alpha luminosity in bridge galaxy in comp2:", lha_high_int_br

    #plt.imshow(av_map_bridge, vmin=0, vmax=4, origin='lower', interpolation='None')
    #plt.show()

    # read in i band SDSS image
    #sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    #lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # ----------------------------------- Get intrinsic fluxes for all lines ----------------------------------- #
    # Factor for converting flux to luminosity for Taffy
    flux_to_lum = 4 * np.pi * (63.2 * 1e6 * 3.08567758128e18)**2  # lum dist to Taffy assumed to be 63.2 Mpc

    # Apply only region masks # NO checkerboard mask # See written notes for details
    # Re-reading these in because I modifie the earlier ones with the checkerboard mask
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    for i in range(len(line_name_list)):

        # Get observed flux for line
        # Read in line map from LZIFU
        line_map = fits.open(taffy_extdir + 'stitched_cube_' + line_name_list[i] + '.fits')
        line_total = line_map[0].data[0]

        # Only consider valid and infinite indices
        i0 = np.where(line_total != -9999.0)
        val_line_total = line_total[i0]  # line_total[i0] gives a 1D array # Use line_total[i0[0]][i0[1]] to get the proper 2D array
        print "\n", "Total", line_name_list[i], "observed luminosity:", "{:.3e}".format(np.sum(val_line_total[np.isfinite(val_line_total)]) * 1e-18 * flux_to_lum)

        # Fluxes from regions
        # --------- total --------- #
        line_total_br = ma.array(line_total, mask=bridge_mask)
        line_total_n = ma.array(line_total, mask=north_mask)
        line_total_s = ma.array(line_total, mask=south_mask)

        i0_br = np.where(line_total_br != -9999.0)
        val_line_total_br = line_total_br[i0_br]
        line_br = np.sum(val_line_total_br[np.isfinite(val_line_total_br)]) * 1e-18 * flux_to_lum
        print "Total", line_name_list[i], "observed luminosity from bridge:", "{:.3e}".format(line_br)

        i0_n = np.where(line_total_n != -9999.0)
        val_line_total_n = line_total_n[i0_n]
        line_n = np.sum(val_line_total_n[np.isfinite(val_line_total_n)]) * 1e-18 * flux_to_lum
        print "Total", line_name_list[i], "observed luminosity from north:", "{:.3e}".format(line_n)

        i0_s = np.where(line_total_s != -9999.0)
        val_line_total_s = line_total_s[i0_s]
        line_s = np.sum(val_line_total_s[np.isfinite(val_line_total_s)]) * 1e-18 * flux_to_lum
        print "Total", line_name_list[i], "observed luminosity from south:", "{:.3e}".format(line_s)

        print "Sum of observed luminosity from regions:", "{:.3e}".format(line_br + line_n + line_s)
        # I think that when these are summed they are just less than the total printed above
        # because the regions I've defined don't exactly cover everything. I do beleive that 
        # the sum as printed here should be more reliable than the simplistic total above.

    sys.exit(0)

    # ----------------------------------- Plotting ----------------------------------- #
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