from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

def correct_spaxel_wise(all_extinction_maps, line_name_list):

    # Apply only region masks # NO checkerboard mask # See written notes for details
    # Re-reading these in because I modifie the earlier ones with the checkerboard mask
    bridge_mask_noCB = vcm.get_region_mask('bridge_bpt_new')
    north_mask_noCB = vcm.get_region_mask('north_galaxy_bpt')
    south_mask_noCB = vcm.get_region_mask('south_galaxy_bpt')

    for i in range(len(line_name_list)):

        current_ext_map = all_extinction_maps[i]

        # Read in maps from LZIFU
        line_map = fits.open(taffy_extdir + 'stitched_cube_' + line_name_list[i] + '.fits')
        line_map_err = fits.open(taffy_extdir + 'stitched_cube_' + line_name_list[i] + '_ERR.fits')
        line_total = line_map[0].data[0]
        line_total_err = line_map_err[0].data[0]

        # Remove the -9999.0 values before you do the correction 
        # because they won't remain -9999.0 after teh correction is multiplied.
        line_total[np.where(line_total == -9999.0)] = np.nan
        line_total_err[np.where(line_total_err == -9999.0)] = np.nan

        # Get ext corr line map and error map
        ext_corr_line_map = line_total * 10**(0.4 * current_ext_map)
        ext_corr_line_err = line_total_err * 10**(0.4 * current_ext_map)

        ext_corr_line_map *= 1e-5
        ext_corr_line_err *= 1e-5

        # Fluxes from regions
        # --------- total --------- #
        line_total_br = ma.array(ext_corr_line_map, mask=bridge_mask_noCB)
        line_total_n = ma.array(ext_corr_line_map, mask=north_mask_noCB)
        line_total_s = ma.array(ext_corr_line_map, mask=south_mask_noCB)

        # --------- total errors --------- #
        line_total_br_err = ma.array(ext_corr_line_err, mask=bridge_mask_noCB)
        line_total_n_err = ma.array(ext_corr_line_err, mask=north_mask_noCB)
        line_total_s_err = ma.array(ext_corr_line_err, mask=south_mask_noCB)

        # -------------------------------------------
        # ------- North
        i0_n = np.where(line_total_n != -9999.0)
        val_line_total_n = line_total_n[i0_n]
        line_n = np.sum(val_line_total_n[np.isfinite(val_line_total_n)])

        # ------- North Errors
        i0_n = np.where(line_total_n_err != -9999.0)
        val_line_total_n_err = line_total_n_err[i0_n]
        line_n_err = np.sum(val_line_total_n_err[np.isfinite(val_line_total_n_err)])
        print "\n"
        print "Total", line_name_list[i], "intrinsic flux from north:", "{:.2f}".format(line_n), "+-", "{:.2f}".format(line_n_err)

        # ------- South
        i0_s = np.where(line_total_s != -9999.0)
        val_line_total_s = line_total_s[i0_s]
        line_s = np.sum(val_line_total_s[np.isfinite(val_line_total_s)])

        # ------- South Errors
        i0_s = np.where(line_total_s_err != -9999.0)
        val_line_total_s_err = line_total_s_err[i0_s]
        line_s_err = np.sum(val_line_total_s_err[np.isfinite(val_line_total_s_err)])
        print "Total", line_name_list[i], "intrinsic flux from south:", "{:.2f}".format(line_s), "+-", "{:.2f}".format(line_s_err)

        # ------- Bridge
        i0_br = np.where(line_total_br != -9999.0)
        val_line_total_br = line_total_br[i0_br]
        line_br = np.sum(val_line_total_br[np.isfinite(val_line_total_br)])

        # ------- Bridge Error
        i0_br = np.where(line_total_br_err != -9999.0)
        val_line_total_br_err = line_total_br_err[i0_br]
        line_br_err = np.sum(val_line_total_br_err[np.isfinite(val_line_total_br_err)])
        print "Total", line_name_list[i], "intrinsic flux from bridge:", "{:.2f}".format(line_br), "+-", "{:.2f}".format(line_br_err)

    return None

def get_halpha_lum(stitched_cube, a_Halpha_map):

    # Factor for converting flux to luminosity for Taffy
    flux_to_lum = 4 * np.pi * (63.2 * 1e6 * 3.08567758128e18)**2  # lum dist to Taffy assumed to be 63.2 Mpc

    # -----------------
    # Read in total and indivdual vel comp
    halpha = stitched_cube['HALPHA'].data[0]
    halpha_comp1 = stitched_cube['HALPHA'].data[1]
    halpha_comp2 = stitched_cube['HALPHA'].data[2]
    
    # Read in total and indivdual vel comp ERRORS
    halpha_err = stitched_cube['HALPHA_ERR'].data[0]
    halpha_comp1_err = stitched_cube['HALPHA_ERR'].data[1]
    halpha_comp2_err = stitched_cube['HALPHA_ERR'].data[2]

    # -----------------
    # Remove the -9999.0 values before you do the correction 
    # because they won't remain -9999.0 after teh correction is multiplied.
    halpha[np.where(halpha == -9999.0)] = np.nan
    halpha_comp1[np.where(halpha_comp1 == -9999.0)] = np.nan
    halpha_comp2[np.where(halpha_comp2 == -9999.0)] = np.nan
    halpha_err[np.where(halpha_err == -9999.0)] = np.nan
    halpha_comp1_err[np.where(halpha_comp1_err == -9999.0)] = np.nan
    halpha_comp2_err[np.where(halpha_comp2_err == -9999.0)] = np.nan

    # -----------------
    # Get ext corr line map and error map
    ext_corr_halpha = halpha * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum
    ext_corr_halpha_comp1 = halpha_comp1 * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum
    ext_corr_halpha_comp2 = halpha_comp2 * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum

    ext_corr_halpha_err = halpha_err * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum
    ext_corr_halpha_comp1_err = halpha_comp1_err * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum
    ext_corr_halpha_comp2_err = halpha_comp2_err * 10**(0.4 * a_Halpha_map) * 1e-18 * flux_to_lum

    # -----------------
    # Print total luminosity after getting valid and finite indices
    val_halpha = ext_corr_halpha[np.where(ext_corr_halpha != -9999.0)]
    val_halpha_err = ext_corr_halpha_err[np.where(ext_corr_halpha_err != -9999.0)]
    print "Total H-alpha luminosity:", "{:.3e}".format(np.sum(val_halpha[np.isfinite(val_halpha)])), "+-",
    print "{:.3e}".format(np.sum(val_halpha_err[np.isfinite(val_halpha_err)]))

    val_halpha_comp1 = ext_corr_halpha_comp1[np.where(ext_corr_halpha_comp1 != -9999.0)]
    val_halpha_comp1_err = ext_corr_halpha_comp1_err[np.where(ext_corr_halpha_comp1_err != -9999.0)]
    print "Total H-alpha comp1 luminosity:", "{:.3e}".format(np.sum(val_halpha_comp1[np.isfinite(val_halpha_comp1)])), "+-",
    print "{:.3e}".format(np.sum(val_halpha_comp1_err[np.isfinite(val_halpha_comp1_err)]))
    lum_ha_low_int = np.sum(val_halpha_comp1[np.isfinite(val_halpha_comp1)])
    lum_ha_low_int_err = np.sum(val_halpha_comp1_err[np.isfinite(val_halpha_comp1_err)])

    val_halpha_comp2 = ext_corr_halpha_comp2[np.where(ext_corr_halpha_comp2 != -9999.0)]
    val_halpha_comp2_err = ext_corr_halpha_comp2_err[np.where(ext_corr_halpha_comp2_err != -9999.0)]
    print "Total H-alpha comp2 luminosity:", "{:.3e}".format(np.sum(val_halpha_comp2[np.isfinite(val_halpha_comp2)])), "+-",
    print "{:.3e}".format(np.sum(val_halpha_comp2_err[np.isfinite(val_halpha_comp2_err)]))
    lum_ha_high_int = np.sum(val_halpha_comp2[np.isfinite(val_halpha_comp2)])
    lum_ha_high_int_err = np.sum(val_halpha_comp2_err[np.isfinite(val_halpha_comp2_err)])

    # -----------------
    # SF Luminosity
    lum_ha_from_sf_int = 0.63*lum_ha_low_int + 0.45*lum_ha_high_int
    lum_ha_from_sf_int_err = 0.63*lum_ha_low_int_err + 0.45*lum_ha_high_int_err
    print lum_ha_from_sf_int, "+-", lum_ha_from_sf_int_err

    # From the bridge
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')

    ext_corr_halpha_comp1_br = ma.array(ext_corr_halpha_comp1, mask=bridge_mask)
    ext_corr_halpha_comp2_br = ma.array(ext_corr_halpha_comp2, mask=bridge_mask)
    ext_corr_halpha_comp1_br_err = ma.array(ext_corr_halpha_comp1_err, mask=bridge_mask)
    ext_corr_halpha_comp2_br_err = ma.array(ext_corr_halpha_comp2_err, mask=bridge_mask)

    val_halpha_br_comp1 = ext_corr_halpha_comp1_br[np.where(ext_corr_halpha_comp1_br != -9999.0)]
    val_halpha_br_comp1_err = ext_corr_halpha_comp1_br_err[np.where(ext_corr_halpha_comp1_br_err != -9999.0)]
    print "Total H-alpha comp1 luminosity in bridge:", "{:.3e}".format(np.sum(val_halpha_br_comp1[np.isfinite(val_halpha_br_comp1)])), "+-",
    print "{:.3e}".format(np.sum(val_halpha_br_comp1_err[np.isfinite(val_halpha_br_comp1_err)]))
    lum_ha_br_low_int = np.sum(val_halpha_br_comp1[np.isfinite(val_halpha_br_comp1)])
    lum_ha_br_low_int_err = np.sum(val_halpha_br_comp1_err[np.isfinite(val_halpha_br_comp1_err)])

    val_halpha_br_comp2 = ext_corr_halpha_comp2_br[np.where(ext_corr_halpha_comp2_br != -9999.0)]
    val_halpha_br_comp2_err = ext_corr_halpha_comp2_br_err[np.where(ext_corr_halpha_comp2_br_err != -9999.0)]
    print "Total H-alpha comp2 luminosity in bridge:", "{:.3e}".format(np.sum(val_halpha_br_comp2[np.isfinite(val_halpha_br_comp2)])), "+-",
    print "{:.3e}".format(np.sum(val_halpha_br_comp2_err[np.isfinite(val_halpha_br_comp2_err)]))
    lum_ha_br_high_int = np.sum(val_halpha_br_comp2[np.isfinite(val_halpha_br_comp2)])
    lum_ha_br_high_int_err = np.sum(val_halpha_br_comp2_err[np.isfinite(val_halpha_br_comp2_err)])

    # SF Luminosity
    lum_ha_br_from_sf_int = 0.63*lum_ha_br_low_int + 0.45*lum_ha_br_high_int
    lum_ha_br_from_sf_int_err = 0.63*lum_ha_br_low_int_err + 0.45*lum_ha_br_high_int_err
    print lum_ha_br_from_sf_int, "+-", lum_ha_br_from_sf_int_err

    return None

def get_wav_arr():

    # create wavelength array
    # I read these data from the corresponding headers
    delt_b = 0.3  # i.e. the wav axis is sampled at 0.3A
    blue_wav_start = 4662.0
    total_blue_res_elem = 2227
    blue_res = 1.6

    blue_wav_arr = [blue_wav_start + delt_b*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350
    red_res = 1.5

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    return blue_wav_arr, red_wav_arr

def get_pseudo_cont_level(obs_data, line_idx, linepad_left, linepad_right):

    pseudo_cont_arr_left = obs_data[line_idx-10-(linepad_left-10):line_idx-10, i, j]
    pseudo_cont_arr_right = obs_data[line_idx+10:line_idx+10+(linepad_right-10), i, j]
    cont_level = np.nanmean(np.concatenate((pseudo_cont_arr_left, pseudo_cont_arr_right)))

    return cont_level

if __name__ == '__main__':
    """
    To get extinction corrected emission line fluxes.

    Getting the emission line fluxes by averaging A_lambda over regions 
    and then applying that average to luminosity averaged over regions
    versus the spaxel-by-spaxel way will give slightly different results.

    For the former method:
    The A_lambda average is computed from a map that HAS the checkerboard mask 
    applied. This A_lambda average is then used with a line_map that DOES NOT
    have the checkerboard (CB) mask applied. Since we confirmed that the original
    spatial smoothing conserved the flux, I don't think we should be applying
    the CB mask anywhere in this correction.

    For the latter method:
    The luminosities are corrected spaxel-by-spaxel with no regard to any mask.
    The luminosities from different regions are then computed by applying the
    individual region masks with NO CB mask applied.
    This method is used in the correct_spaxel_wise() function above.

    It will also be different simply because you're applying an average value
    of dust extinction to another value averaged over the region which ignores
    any small scale over/under-densities and only gives an average result.

    ----------------------
    Now to find the ratio and get dust extinction
    assumes that all the halpha and hbeta comes from
    case B which is not the case. We are explicitly saying 
    htat some of the ionized gas emission comes from shocks.
    you will have to break down the halpha and Hbeta emission 
    into the fractions coming from shocks and from SF separately.

    I do not think there is a way around this problem. To get the 
    luminosity of any line that comes just from star-formation 
    we need the extinction corrected values which we cannot do
    without first doing this calculation. I do think the problem
    is mitigated because this is a ratio but even so we should
    mention the caveat (perhaps after the referee report) that 
    our extinction values might be slightly overestimated.

    loop over all pixels and get a ebv value for each
    ebv[i,j] = 1.97 * np.log10(halpha[i,j]/hbeta[i,j] / 2.86)
    """

    # ----------------------------------- Read in arrays ----------------------------------- #
    # Read in stitched cube
    #stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

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

    # ------- Read lines to plot when checking individual fits --------- #
    b_line = fits.open(taffy_extdir + 'stitched_cube_B_LINE.fits')
    r_line = fits.open(taffy_extdir + 'stitched_cube_R_LINE.fits')
    # read in observed data
    obs_b = fits.open(taffy_data + 'Taffy_B.fits')
    obs_r = fits.open(taffy_data + 'Taffy_R.fits')
    obs_data_b = obs_b[0].data
    obs_data_r = obs_r[0].data

    # find line index in wavelength array
    redshift = 0.0145  # average z
    hbeta_air_wav = 4861.363
    halpha_air_wav = 6562.80

    blue_wav_arr, red_wav_arr = get_wav_arr()

    hbeta_wav = hbeta_air_wav*(1+redshift)
    hbeta_idx = np.argmin(abs(blue_wav_arr - hbeta_wav))

    halpha_wav = halpha_air_wav*(1+redshift)
    halpha_idx = np.argmin(abs(red_wav_arr - halpha_wav))

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
    ebv_map = np.zeros((58,58))

    # conv ds9 coords to array coords 
    # to be able to check with ds9
    pix_x = 39
    pix_y = 19
    arr_x = pix_y - 1
    arr_y = pix_x - 1
    box_size = 5

    # If you want to analyze a block enter the ds9 pix coords of the low 
    # left corner above. Otherwise you need to loop both i and j over range(58)
    # Also uncomment plotting below if you want to check individual fits
    for i in range(arr_x, arr_x + box_size):
        for j in range(arr_y, arr_y + box_size):

            ebv_map[i,j] = 1.97 * np.log10(halpha_total[i,j]/hbeta_total[i,j] / 2.86)

            print j+1, i+1, "{:.2f}".format(halpha_total[i,j]), "{:.2f}".format(hbeta_total[i,j]), \
            "{:.2f}".format(ebv_map[i,j]), "{:.2f}".format(4.05 * ebv_map[i,j])

            # PLot to check
            fig = plt.figure()
            gs = gridspec.GridSpec(5,10)
            gs.update(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=1.0, hspace=0)

            ax1 = fig.add_subplot(gs[:,:5])
            ax2 = fig.add_subplot(gs[:,5:])

            linepad_left = 70
            linepad_right = 70
            
            hbeta_y_arr_data = obs_data_b[hbeta_idx-linepad_left:hbeta_idx+linepad_right, i, j]
            hbeta_x_arr_data = np.linspace(hbeta_idx-linepad_left, hbeta_idx+linepad_right, len(hbeta_y_arr_data))

            halpha_y_arr_data = obs_data_r[halpha_idx-linepad_left:halpha_idx+linepad_right, i, j]
            halpha_x_arr_data = np.linspace(halpha_idx-linepad_left, halpha_idx+linepad_right, len(halpha_y_arr_data))

            # find pseudo continuum and subtract
            hbeta_cont_level= get_pseudo_cont_level(obs_data_b, hbeta_idx, linepad_left, linepad_right)
            halpha_cont_level= get_pseudo_cont_level(obs_data_r, halpha_idx, linepad_left, linepad_right)
            hbeta_y_arr_data -= hbeta_cont_level
            halpha_y_arr_data -= halpha_cont_level

            ax1.plot(hbeta_x_arr_data, hbeta_y_arr_data, color='midnightblue')
            ax2.plot(halpha_x_arr_data, halpha_y_arr_data, color='maroon')

            plt.show()
            sys.exit()

    sys.exit(0)

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

    all_extinction_maps = [a_Hbeta_map, a_OIII_5007_map, a_OI_6300_map, a_Halpha_map, a_NII_6583_map, a_SII_6716_map, a_SII_6731_map]

    correct_spaxel_wise(all_extinction_maps, line_name_list)
    get_halpha_lum(stitched_cube, a_Halpha_map)
    sys.exit(0)

    # --------------------------------- Get average extinction values in the different regions and vel comps --------------------------------- #
    # get the region masks first
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # Combine all masks with checker board mask
    # See checkerboard mask comment in BPT velo comp code
    """
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
    """

    # Now apply masks
    av_map_bridge = ma.array(av_map, mask=bridge_mask)
    av_map_north = ma.array(av_map, mask=north_mask)
    av_map_south = ma.array(av_map, mask=south_mask)

    print "\n", "Mean visual extinction in north galaxy:", np.nanmean(av_map_north)
    print "Mean visual extinction in south galaxy:", np.nanmean(av_map_south)
    print "Mean visual extinction in bridge:", np.nanmean(av_map_bridge)

    # ------------- Extinction at Halpha -------------- #
    """
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

    # Errors:
    lha_low_obs_n_err = 1.097e+39
    lha_low_obs_s_err = 3.091e+38
    lha_low_obs_br_err = 1.168e+39

    lha_high_obs_n_err = 1.096e+39
    lha_high_obs_s_err = 3.145e+38
    lha_high_obs_br_err = 1.168e+39

    lum_ha_low_int = lha_low_obs_n * 10**(0.4*aha_n) + lha_low_obs_s * 10**(0.4*aha_s) + lha_low_obs_br * 10**(0.4*aha_br)
    lum_ha_high_int = lha_high_obs_n * 10**(0.4*aha_n) + lha_high_obs_s * 10**(0.4*aha_s) + lha_high_obs_br * 10**(0.4*aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity in the low vel comp:", "{:.3e}".format(lum_ha_low_int)
    print "Intrinsic H-alpha luminosity in the high vel comp:", "{:.3e}".format(lum_ha_high_int)
    lum_ha_from_sf_int = 0.63*lum_ha_low_int + 0.45*lum_ha_high_int
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

    # ----------------- Errors ----------------- #
    lum_ha_low_int_err = lha_low_obs_n_err * 10**(0.4*aha_n) + lha_low_obs_s_err * 10**(0.4*aha_s) + lha_low_obs_br_err * 10**(0.4*aha_br)
    lum_ha_high_int_err = lha_high_obs_n_err * 10**(0.4*aha_n) + lha_high_obs_s_err * 10**(0.4*aha_s) + lha_high_obs_br_err * 10**(0.4*aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity error in the low vel comp:", "{:.3e}".format(lum_ha_low_int_err)
    print "Intrinsic H-alpha luminosity error in the high vel comp:", "{:.3e}".format(lum_ha_high_int_err)
    lum_ha_from_sf_int_err = 0.63*lum_ha_low_int_err + 0.45*lum_ha_high_int_err
    print "Total intrinsic H-alpha luminosity error from star-formation:", "{:.3e}".format(lum_ha_from_sf_int_err)

    # ---------------------------- components
    lha_low_int_n_err = lha_low_obs_n_err * 10**(0.4 * aha_n)
    lha_low_int_s_err = lha_low_obs_s_err * 10**(0.4 * aha_s)
    lha_low_int_br_err = lha_low_obs_br_err * 10**(0.4 * aha_br)

    lha_high_int_n_err = lha_high_obs_n_err * 10**(0.4 * aha_n)
    lha_high_int_s_err = lha_high_obs_s_err * 10**(0.4 * aha_s)
    lha_high_int_br_err = lha_high_obs_br_err * 10**(0.4 * aha_br)

    print "\n"
    print "Intrinsic H-alpha luminosity error in north galaxy in comp1:", lha_low_int_n_err
    print "Intrinsic H-alpha luminosity error in south galaxy in comp1:", lha_low_int_s_err
    print "Intrinsic H-alpha luminosity error in bridge galaxy in comp1:", lha_low_int_br_err

    print "\n"
    print "Intrinsic H-alpha luminosity error in north galaxy in comp2:", lha_high_int_n_err
    print "Intrinsic H-alpha luminosity error in south galaxy in comp2:", lha_high_int_s_err
    print "Intrinsic H-alpha luminosity error in bridge galaxy in comp2:", lha_high_int_br_err

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
    bridge_mask_noCB = vcm.get_region_mask('bridge_bpt_new')
    north_mask_noCB = vcm.get_region_mask('north_galaxy_bpt')
    south_mask_noCB = vcm.get_region_mask('south_galaxy_bpt')

    for i in range(len(line_name_list)):

        a_line_map = all_extinction_maps[i]

        # Get extinction values for line for the different regions # Uses CB mask
        a_line_map_bridge = ma.array(a_line_map, mask=bridge_mask)
        a_line_map_north = ma.array(a_line_map, mask=north_mask)
        a_line_map_south = ma.array(a_line_map, mask=south_mask)

        line_ext_n = np.nanmean(a_line_map_north)
        line_ext_s = np.nanmean(a_line_map_south)
        line_ext_br = np.nanmean(a_line_map_bridge)

        # Get observed flux and flux errors for line
        # Read in maps from LZIFU
        line_map = fits.open(taffy_extdir + 'stitched_cube_' + line_name_list[i] + '.fits')
        line_map_err = fits.open(taffy_extdir + 'stitched_cube_' + line_name_list[i] + '_ERR.fits')
        line_total = line_map[0].data[0]
        line_total_err = line_map_err[0].data[0]

        # Only consider valid and infinite indices
        i0 = np.where(line_total != -9999.0)
        val_line_total = line_total[i0]  # line_total[i0] gives a 1D array # Use line_total[i0[0]][i0[1]] to get the proper 2D array
        print "\n", "Total", line_name_list[i], "observed luminosity (no region mask applied):", \
        "{:.3e}".format(np.sum(val_line_total[np.isfinite(val_line_total)]) * 1e-18 * flux_to_lum)

        # Fluxes from regions
        # --------- total --------- #
        line_total_br = ma.array(line_total, mask=bridge_mask_noCB)
        line_total_n = ma.array(line_total, mask=north_mask_noCB)
        line_total_s = ma.array(line_total, mask=south_mask_noCB)

        # --------- total errors --------- #
        line_total_br_err = ma.array(line_total_err, mask=bridge_mask_noCB)
        line_total_n_err = ma.array(line_total_err, mask=north_mask_noCB)
        line_total_s_err = ma.array(line_total_err, mask=south_mask_noCB)

        # --------------------
        # ------- Bridge
        i0_br = np.where(line_total_br != -9999.0)
        val_line_total_br = line_total_br[i0_br]
        line_br = np.sum(val_line_total_br[np.isfinite(val_line_total_br)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_br_int = line_br * 10**(0.4 * line_ext_br)
        print "Total", line_name_list[i], "observed and intrinsic luminosities from bridge:", \
        "{:.3e}".format(line_br), "{:.3e}".format(line_br_int)

        # ------- Bridge Error
        i0_br = np.where(line_total_br_err != -9999.0)
        val_line_total_br_err = line_total_br_err[i0_br]
        line_br_err = np.sum(val_line_total_br_err[np.isfinite(val_line_total_br_err)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_br_int_err = line_br_err * 10**(0.4 * line_ext_br)
        print "Total", line_name_list[i], "observed and intrinsic luminosity errors from bridge:", \
        "{:.3e}".format(line_br_err), "{:.3e}".format(line_br_int_err)

        # --------------------
        # ------- North
        i0_n = np.where(line_total_n != -9999.0)
        val_line_total_n = line_total_n[i0_n]
        line_n = np.sum(val_line_total_n[np.isfinite(val_line_total_n)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_n_int = line_n * 10**(0.4 * line_ext_n)
        print "Total", line_name_list[i], "observed and intrinsic luminosities from north:", "{:.3e}".format(line_n), "{:.3e}".format(line_n_int)

        # ------- North Errors
        i0_n = np.where(line_total_n_err != -9999.0)
        val_line_total_n_err = line_total_n_err[i0_n]
        line_n_err = np.sum(val_line_total_n_err[np.isfinite(val_line_total_n_err)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_n_int_err = line_n_err * 10**(0.4 * line_ext_n)
        print "Total", line_name_list[i], "observed and intrinsic luminosity errors from north:", \
        "{:.3e}".format(line_n_err), "{:.3e}".format(line_n_int_err)

        # --------------------
        # ------- South
        i0_s = np.where(line_total_s != -9999.0)
        val_line_total_s = line_total_s[i0_s]
        line_s = np.sum(val_line_total_s[np.isfinite(val_line_total_s)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_s_int = line_s * 10**(0.4 * line_ext_s)
        print "Total", line_name_list[i], "observed and intrinsic luminosities from south:", "{:.3e}".format(line_s), "{:.3e}".format(line_s_int)

        # ------- South Errors
        i0_s = np.where(line_total_s_err != -9999.0)
        val_line_total_s_err = line_total_s_err[i0_s]
        line_s_err = np.sum(val_line_total_s_err[np.isfinite(val_line_total_s_err)]) * 1e-18 * flux_to_lum
        # Intrinsic
        line_s_int_err = line_s_err * 10**(0.4 * line_ext_s)
        print "Total", line_name_list[i], "observed and intrinsic luminosity errors from south:", \
        "{:.3e}".format(line_s_err), "{:.3e}".format(line_s_int_err)

        # --------------------
        print "Sum of observed and intrinsic luminosities from regions:", \
        "{:.3e}".format(line_br + line_n + line_s), "{:.3e}".format(line_br_int + line_n_int + line_s_int)
        # I think that when these are summed they are just less than the total printed above
        # because the regions I've defined don't exactly cover everything. I do beleive that 
        # the sum as printed here should be more reliable than the simplistic total above.

        # Output to copy-paste in tex file
        print "TeX for line (these are now fluxes [W/m^2] not luminosities)",
        print line_name_list[i], ":  ",
        print "{:.3}".format(line_n_int * 1e-3 / (flux_to_lum * 1e-16)), "$\pm$", "{:.2}".format(line_n_int_err * 1e-3 / (flux_to_lum * 1e-16)),
        print "&", "{:.3}".format(line_s_int * 1e-3 / (flux_to_lum * 1e-16)), "$\pm$", "{:.2}".format(line_s_int_err * 1e-3 / (flux_to_lum * 1e-16)),
        print "&", "{:.3}".format(line_br_int * 1e-3 / (flux_to_lum * 1e-16)), "$\pm$", "{:.2}".format(line_br_int_err * 1e-3 / (flux_to_lum * 1e-16)), 
        print "\\\\"

        # Note on unit conversion above:
        # I divided by the 4*pi*d_L^2 to go from luminosity to flux.
        # Then the 1e-3 converts that to W m^-2 and the 1e-16 further will write it in units of [W/m^2]x1e-16.

    sys.exit(0)
    """

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