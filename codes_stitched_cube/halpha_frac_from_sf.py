from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import bpt_plots as bpt

def plot_bpt_with_hii_shaded(plottype, vel_comp, xarr_br, xarr_n, xarr_s, yarr_br, yarr_n, yarr_s, \
    xarr_err_br, xarr_err_n, xarr_err_s, yarr_err_br, yarr_err_n, yarr_err_s, \
    xarr_snuc, yarr_snuc, xarr_nw, yarr_nw, xarr_nb, yarr_nb, \
    valid_indices, figdir, yerr_avg):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)
    ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)

    # read in Mappings III models and overplot
    mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
    mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
    mappings_oi_halpha_v250, mappings_oi_halpha_v300, mappings_oi_halpha_v500,\
    mappings_oi_halpha_v800,\
    mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
    mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
    mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
    mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
    mappings_nii_halpha_v800,\
    mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
    mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
    mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
    mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
    mappings_oiii_hbeta_v800,\
    mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
    mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
    mappings_sii_halpha_v250, mappings_sii_halpha_v300, mappings_sii_halpha_v500,\
    mappings_sii_halpha_v800 = bpt.mappings_oi_nii_sii()

    ax.errorbar(xarr_br[valid_indices], yarr_br[valid_indices], \
        xerr=xarr_err_br[valid_indices], yerr=yarr_err_br[valid_indices], \
        color='maroon', markersize=6, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.2)
    ax.errorbar(xarr_n[valid_indices], yarr_n[valid_indices], \
        xerr=xarr_err_n[valid_indices], yerr=yarr_err_n[valid_indices], \
        color='darkgreen', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.2)
    ax.errorbar(xarr_s[valid_indices], yarr_s[valid_indices], \
        xerr=xarr_err_s[valid_indices], yerr=yarr_err_s[valid_indices], \
        color='midnightblue', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.2)

    # Circle interesting regions
    ax.scatter(xarr_snuc[valid_indices], yarr_snuc[valid_indices], s=30, marker='d', edgecolors='midnightblue', facecolors='midnightblue')
    ax.scatter(xarr_nw[valid_indices], yarr_nw[valid_indices], s=50, edgecolors='darkorchid', facecolors='none', zorder=5)
    ax.scatter(xarr_nb[valid_indices], yarr_nb[valid_indices], s=50, lw=1.5, edgecolors='darkorange', facecolors='none', zorder=5)

    """
    All of the BPT classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    """
    # plot classification lines
    # HII line
    x_arr = np.arange(-1.5, 0.0, 0.01)
    y_arr = 1.3 + 0.61 / (x_arr - 0.05)

    ax.plot(x_arr, y_arr, '-', color='k')

    # AGN/Seyfert line
    y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)
    ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')
    
    # Plot grey padded region for getting H-alpha SF excited fraction
    y_upper_hii_line = y_arr + 2*yerr_avg
    ax.fill_between(x_arr, -2.0, y_upper_hii_line, color='lightgray', alpha=0.5)

    # Plot shock models
    ax.plot(mappings_nii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s')
    ax.plot(mappings_nii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s')
    ax.plot(mappings_nii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s')
    ax.plot(mappings_nii_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s')
    ax.plot(mappings_nii_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s')
    ax.plot(mappings_nii_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s')

    ax.set_xlim(-1,0.3)
    ax.set_ylim(-1,1)

    # labels
    agnbox = TextArea('AGN', textprops=dict(color='k', size=16))
    anc_agnbox = AnchoredOffsetbox(loc=2, child=agnbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.57, 0.93),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_agnbox) 

    compbox = TextArea('HII-AGN Composite', textprops=dict(color='k', size=16))
    anc_compbox = AnchoredOffsetbox(loc=2, child=compbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.56, 0.1),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_compbox) 

    hiibox = TextArea('HII', textprops=dict(color='k', size=16))
    anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.32, 0.3),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_hiibox)

    # Other auxilliary plot commands
    ax.legend(loc='center left', prop={'size':10})

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    # save and close
    fig.savefig(figdir + 'bpt_' + plottype + '_comp' + vel_comp + '_errbar_hii_shaded.png', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    return None

def calc_lower_lim(nii_halpha_withcut, xarr_br, xarr_n, xarr_s, yarr_br, yarr_n, yarr_s, \
    xarr_err_br, xarr_err_n, xarr_err_s, yarr_err_br, yarr_err_n, yarr_err_s, \
    halpha_withcut_bridge, halpha_withcut_north, halpha_withcut_south):

    # ---------------------- Now get the padding area around the HII line ---------------------- #
    # these classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
    y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

    # get non-zero indices
    nii_nonzero = np.nonzero(nii_halpha_withcut)

    # Get average values of the errors to get padding area 
    xerr_avg_arr = ma.concatenate((xarr_err_br[nii_nonzero], xarr_err_n[nii_nonzero], xarr_err_s[nii_nonzero]))
    yerr_avg_arr = ma.concatenate((yarr_err_br[nii_nonzero], yarr_err_n[nii_nonzero], yarr_err_s[nii_nonzero]))
    
    # Nan out values = -9999.0
    invalid_idx_xerr = np.where(xerr_avg_arr == -9999.0)
    invalid_idx_yerr = np.where(yerr_avg_arr == -9999.0)

    xerr_avg_arr[invalid_idx_xerr] = np.nan
    yerr_avg_arr[invalid_idx_yerr] = np.nan

    xerr_avg = np.nanmean(xerr_avg_arr)
    yerr_avg = np.nanmean(yerr_avg_arr)

    print "Average error on log([NII]/Ha):", xerr_avg
    print "Average error on log([OIII]/Hb):", yerr_avg

    # Find total h-alpha coming from spaxels within this grey region
    # loop over all points and check where they are
    hii_x = []
    hii_y = []
    full_x_arr = ma.concatenate((xarr_br[nii_nonzero], xarr_n[nii_nonzero], xarr_s[nii_nonzero]))
    full_y_arr = ma.concatenate((yarr_br[nii_nonzero], yarr_n[nii_nonzero], yarr_s[nii_nonzero]))

    full_halpha_arr = ma.concatenate((halpha_withcut_bridge[nii_nonzero], halpha_withcut_north[nii_nonzero], \
        halpha_withcut_south[nii_nonzero]))
    halpha_sum = 0
    for i in range(len(full_x_arr)):

        current_x = full_x_arr[i]
        current_y = full_y_arr[i]

        # Find if the current spaxel is within the HII region 
        y_on_grayline = 1.3 + 0.61 / (current_x - 0.05)  # HII classification line

        if current_y <= y_on_grayline + 2*yerr_avg:
            hii_x.append(current_x)
            hii_y.append(current_y)

            # Now sum up halpha in that spaxel
            if full_halpha_arr[i] != -9999.0:
                halpha_sum += float(full_halpha_arr[i])

    # ax.scatter(hii_x, hii_y, s=5, color='green', zorder=20)
    # Do not remove above line. Useful for checking method to tag 
    # HII excited spaxels.

    val_idx = np.where(full_halpha_arr != -9999.0)[0]
    total_halpha = np.sum(full_halpha_arr[val_idx])
    print "Total H-alpha:", total_halpha
    print "Lower limit to H-alpha from SF:", halpha_sum
    print "Fractional lower limit to H-alpha from SF:", halpha_sum / total_halpha

    return xerr_avg, yerr_avg

def prep_total(stitched_cube):

    # assign line arrays
    # Even though I only need the [NII] arrays 
    # I have to read all of these in so that I
    # can use the functions written in the 
    # bpt_plots.py code.
    halpha = stitched_cube['HALPHA'].data[0]
    hbeta = stitched_cube['HBETA'].data[0]
    nii6583 = stitched_cube['NII6583'].data[0]
    oiii5007 = stitched_cube['OIII5007'].data[0]
    oi6300 = stitched_cube['OI6300'].data[0]
    oi6364 = stitched_cube['OI6364'].data[0]
    sii6716 = stitched_cube['SII6716'].data[0]
    sii6731 = stitched_cube['SII6731'].data[0]

    halpha_err = stitched_cube['HALPHA_ERR'].data[0]
    hbeta_err = stitched_cube['HBETA_ERR'].data[0]
    nii6583_err = stitched_cube['NII6583_ERR'].data[0]
    oiii5007_err = stitched_cube['OIII5007_ERR'].data[0]
    oi6300_err = stitched_cube['OI6300_ERR'].data[0]
    oi6364_err = stitched_cube['OI6364_ERR'].data[0]
    sii6716_err = stitched_cube['SII6716_ERR'].data[0]
    sii6731_err = stitched_cube['SII6731_ERR'].data[0]

    # add lines which are doublets
    sii = sii6716 + sii6731
    sii_err = np.sqrt((sii6716_err)**2 + (sii6731_err)**2)

    # get arrays with significance cut applied
    nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut = \
    bpt.get_arr_withsigcut(3, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err, \
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha.shape)

    # get the region masks
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

    # ------------ apply bridge mask ------------ #
    nii_halpha_withcut_bridge = ma.array(nii_halpha_withcut, mask=bridge_mask)
    oiii_hbeta_for_nii_withcut_bridge = ma.array(oiii_hbeta_for_nii_withcut, mask=bridge_mask)

    # on errors
    nii_halpha_err_withcut_bridge = ma.array(nii_halpha_err_withcut, mask=bridge_mask)
    oiii_hbeta_for_nii_err_withcut_bridge = ma.array(oiii_hbeta_for_nii_err_withcut, mask=bridge_mask)

    # ------------ apply north mask ------------ #
    nii_halpha_withcut_north = ma.array(nii_halpha_withcut, mask=north_mask)
    oiii_hbeta_for_nii_withcut_north = ma.array(oiii_hbeta_for_nii_withcut, mask=north_mask)

    # on errors
    nii_halpha_err_withcut_north = ma.array(nii_halpha_err_withcut, mask=north_mask)
    oiii_hbeta_for_nii_err_withcut_north = ma.array(oiii_hbeta_for_nii_err_withcut, mask=north_mask)

    # ------------ apply south mask ------------ #
    nii_halpha_withcut_south = ma.array(nii_halpha_withcut, mask=south_mask)
    oiii_hbeta_for_nii_withcut_south = ma.array(oiii_hbeta_for_nii_withcut, mask=south_mask)

    # on errors
    nii_halpha_err_withcut_south = ma.array(nii_halpha_err_withcut, mask=south_mask)
    oiii_hbeta_for_nii_err_withcut_south = ma.array(oiii_hbeta_for_nii_err_withcut, mask=south_mask)

    # ----------- only for halpha map ----------- #
    halpha_withcut_bridge = ma.array(halpha_withcut, mask=bridge_mask)
    halpha_withcut_north = ma.array(halpha_withcut, mask=north_mask)
    halpha_withcut_south = ma.array(halpha_withcut, mask=south_mask)

    # read in Mappings III models and overplot
    mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
    mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
    mappings_oi_halpha_v250, mappings_oi_halpha_v300, mappings_oi_halpha_v500,\
    mappings_oi_halpha_v800,\
    mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
    mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
    mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
    mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
    mappings_nii_halpha_v800,\
    mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
    mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
    mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
    mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
    mappings_oiii_hbeta_v800,\
    mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
    mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
    mappings_sii_halpha_v250, mappings_sii_halpha_v300, mappings_sii_halpha_v500,\
    mappings_sii_halpha_v800 = bpt.mappings_oi_nii_sii()

    xerr_avg_total, yerr_avg_total = calc_lower_lim(nii_halpha_withcut, nii_halpha_withcut_bridge, nii_halpha_withcut_north, nii_halpha_withcut_south, \
        oiii_hbeta_for_nii_withcut_bridge, oiii_hbeta_for_nii_withcut_north, oiii_hbeta_for_nii_withcut_south, \
        nii_halpha_err_withcut_bridge, nii_halpha_err_withcut_north, nii_halpha_err_withcut_south, \
        oiii_hbeta_for_nii_err_withcut_bridge, oiii_hbeta_for_nii_err_withcut_north, oiii_hbeta_for_nii_err_withcut_south, \
        halpha_withcut_bridge, halpha_withcut_north, halpha_withcut_south)

    # ------------------------------------------------------------------- #
    # Do not delete. Block useful for debugging.
    """
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.imshow(halpha, origin='lower', vmin=0, vmax=1000)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.imshow(halpha_withcut, origin='lower', vmin=0, vmax=1000)

    #plt.show()

    print '\n'

    print "[NII]/Ha shapes:", nii_halpha_withcut.shape, nii_halpha_withcut_bridge.shape, nii_halpha_withcut_north.shape, nii_halpha_withcut_south.shape
    print "[OIII]/Hb shapes:", oiii_hbeta_for_nii_withcut_bridge.shape, oiii_hbeta_for_nii_withcut_north.shape, oiii_hbeta_for_nii_withcut_south.shape
    print "Ha shapes:", halpha_withcut_bridge.shape, halpha_withcut_north.shape, halpha_withcut_south.shape
    print np.nonzero(nii_halpha_withcut)

    # if I apply the 2comp mask
    halpha_withcut = ma.array(halpha_withcut, mask=comb_mask)
    halpha_withcut_bridge = ma.array(halpha_withcut, mask=bridge_mask)
    halpha_withcut_north = ma.array(halpha_withcut, mask=north_mask)
    halpha_withcut_south = ma.array(halpha_withcut, mask=south_mask)

    val_idx2 = np.where(halpha_withcut != -9999.0)

    print "Sum of orig halpha:", np.sum(halpha[np.isfinite(halpha)])
    print "Sum of halpha with sig cut:", np.sum(halpha_withcut[val_idx2])

    val_br = np.where(halpha_withcut_bridge != -9999.0)
    val_n = np.where(halpha_withcut_north != -9999.0)
    val_s = np.where(halpha_withcut_south != -9999.0)

    print "Sum of halpha in bridge:", np.sum(halpha_withcut_bridge[val_br])
    print "Sum of halpha in north:", np.sum(halpha_withcut_north[val_n])
    print "Sum of halpha in south:", np.sum(halpha_withcut_south[val_s])

    # Shape after concatenating
    nonzero_idx = np.nonzero(nii_halpha_withcut)
    temp_h = ma.concatenate((halpha_withcut_bridge[nonzero_idx], halpha_withcut_north[nonzero_idx], halpha_withcut_south[nonzero_idx]))
    vi = np.where(temp_h != -9999.0)[0]
    print "Sum if ma.concatenate is used instead of np.concatenate:", np.sum(temp_h[vi])
    """
    # ------------------------------------------------------------------- #

    return nii_halpha_withcut_bridge, nii_halpha_withcut_north, nii_halpha_withcut_south, \
    oiii_hbeta_for_nii_withcut_bridge, oiii_hbeta_for_nii_withcut_north, oiii_hbeta_for_nii_withcut_south, \
    nii_halpha_err_withcut_bridge, nii_halpha_err_withcut_north, nii_halpha_err_withcut_south, \
    oiii_hbeta_for_nii_err_withcut_bridge, oiii_hbeta_for_nii_err_withcut_north, oiii_hbeta_for_nii_err_withcut_south, \
    nii_halpha_withcut, xerr_avg_total, yerr_avg_total

def get_total_fluxes(stitched_cube):

    # Factor for converting flux to luminosity for Taffy
    flux_to_lum = 4 * np.pi * (63.2 * 1e6 * 3.08567758128e18)**2  # lum dist to Taffy assumed to be 63.2 Mpc

    # Total and indivdual vel comp
    halpha = stitched_cube['HALPHA'].data[0]
    halpha_comp1 = stitched_cube['HALPHA'].data[1]
    halpha_comp2 = stitched_cube['HALPHA'].data[2]

    # Print total flux after getting valid and finite indices
    i0 = np.where(halpha != -9999.0)
    val_halpha = halpha[i0]  # halpha[i0] gives a 1D array # Use halpha[i0[0]][i0[1]] to get the proper 2D array
    print "\n", "Total H-alpha luminosity:", "{:.3e}".format(np.sum(val_halpha[np.isfinite(val_halpha)]) * 1e-18 * flux_to_lum)

    i1 = np.where(halpha_comp1 != -9999.0)
    val_halpha1 = halpha_comp1[i1]
    print "H-alpha luminosity in comp 1:", "{:.3e}".format(np.sum(val_halpha1[np.isfinite(val_halpha1)]) * 1e-18 * flux_to_lum)

    i2 = np.where(halpha_comp2 != -9999.0)
    val_halpha2 = halpha_comp2[i2]
    print "H-alpha luminosity in comp 2:", "{:.3e}".format(np.sum(val_halpha2[np.isfinite(val_halpha2)]) * 1e-18 * flux_to_lum)

    # Apply only region masks # NO checkerboard mask # See written notes for details
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # Fluxes from regions
    # --------- total --------- #
    halpha_br = ma.array(halpha, mask=bridge_mask)
    halpha_n = ma.array(halpha, mask=north_mask)
    halpha_s = ma.array(halpha, mask=south_mask)

    i0_br = np.where(halpha_br != -9999.0)
    val_halpha_br = halpha_br[i0_br]
    print "\n", "Total H-alpha luminosity from bridge:", "{:.3e}".format(np.sum(val_halpha_br[np.isfinite(val_halpha_br)]) * 1e-18 * flux_to_lum)
    i0_n = np.where(halpha_n != -9999.0)
    val_halpha_n = halpha_n[i0_n]
    print "Total H-alpha luminosity from north:", "{:.3e}".format(np.sum(val_halpha_n[np.isfinite(val_halpha_n)]) * 1e-18 * flux_to_lum)
    i0_s = np.where(halpha_s != -9999.0)
    val_halpha_s = halpha_s[i0_s]
    print "Total H-alpha luminosity from south:", "{:.3e}".format(np.sum(val_halpha_s[np.isfinite(val_halpha_s)]) * 1e-18 * flux_to_lum)

    # --------- comp1 --------- #
    halpha_br1 = ma.array(halpha_comp1, mask=bridge_mask)
    halpha_n1 = ma.array(halpha_comp1, mask=north_mask)
    halpha_s1 = ma.array(halpha_comp1, mask=south_mask)

    i0_br1 = np.where(halpha_br1 != -9999.0)
    val_halpha_br1 = halpha_br1[i0_br1]
    print "\n", "Comp1 H-alpha luminosity from bridge:", "{:.3e}".format(np.sum(val_halpha_br1[np.isfinite(val_halpha_br1)]) * 1e-18 * flux_to_lum)
    i0_n1 = np.where(halpha_n1 != -9999.0)
    val_halpha_n1 = halpha_n1[i0_n1]
    print "Comp1 H-alpha luminosity from north:", "{:.3e}".format(np.sum(val_halpha_n1[np.isfinite(val_halpha_n1)]) * 1e-18 * flux_to_lum)
    i0_s1 = np.where(halpha_s1 != -9999.0)
    val_halpha_s1 = halpha_s1[i0_s1]
    print "Comp1 H-alpha luminosity from south:", "{:.3e}".format(np.sum(val_halpha_s1[np.isfinite(val_halpha_s1)]) * 1e-18 * flux_to_lum)

    # --------- comp2 --------- #
    halpha_br2 = ma.array(halpha_comp2, mask=bridge_mask)
    halpha_n2 = ma.array(halpha_comp2, mask=north_mask)
    halpha_s2 = ma.array(halpha_comp2, mask=south_mask)

    i0_br2 = np.where(halpha_br2 != -9999.0)
    val_halpha_br2 = halpha_br2[i0_br2]
    print "\n", "Comp2 H-alpha luminosity from bridge:", "{:.3e}".format(np.sum(val_halpha_br2[np.isfinite(val_halpha_br2)]) * 1e-18 * flux_to_lum)
    i0_n2 = np.where(halpha_n2 != -9999.0)
    val_halpha_n2 = halpha_n2[i0_n2]
    print "Comp2 H-alpha luminosity from north:", "{:.3e}".format(np.sum(val_halpha_n2[np.isfinite(val_halpha_n2)]) * 1e-18 * flux_to_lum)
    i0_s2 = np.where(halpha_s2 != -9999.0)
    val_halpha_s2 = halpha_s2[i0_s2]
    print "Comp2 H-alpha luminosity from south:", "{:.3e}".format(np.sum(val_halpha_s2[np.isfinite(val_halpha_s2)]) * 1e-18 * flux_to_lum)

    return None

def get_total_errors(stitched_cube):

    print "----------------------------------"
    print "Printing errors in fluxes"

    # Factor for converting flux to luminosity for Taffy
    flux_to_lum = 4 * np.pi * (63.2 * 1e6 * 3.08567758128e18)**2  # lum dist to Taffy assumed to be 63.2 Mpc

    # Total and indivdual vel comp ERRORS
    halpha_err = stitched_cube['HALPHA_ERR'].data[0]
    halpha_comp1_err = stitched_cube['HALPHA_ERR'].data[1]
    halpha_comp2_err = stitched_cube['HALPHA_ERR'].data[2]
    # ------------------------------------

    # Print total flux error after getting valid and finite indices
    # Flux error has to be summed in quadrature 
    i0 = np.where(halpha_err != -9999.0)
    val_halpha_err = halpha_err[i0]  # halpha[i0] gives a 1D array # Use halpha[i0[0]][i0[1]] to get the proper 2D array
    print "\n", "Total H-alpha luminosity error:", \
    "{:.3e}".format(np.sqrt(np.sum( val_halpha_err[np.isfinite(val_halpha_err)]**2 )) * 1e-18 * flux_to_lum)

    i1 = np.where(halpha_comp1_err != -9999.0)
    val_halpha1_err = halpha_comp1_err[i1]
    print "H-alpha luminosity error in comp 1:", \
    "{:.3e}".format(np.sqrt(np.sum( val_halpha1_err[np.isfinite(val_halpha1_err)]**2 )) * 1e-18 * flux_to_lum)

    i2 = np.where(halpha_comp2_err != -9999.0)
    val_halpha2_err = halpha_comp2_err[i2]
    print "H-alpha luminosity error in comp 2:", \
    "{:.3e}".format(np.sqrt(np.sum( val_halpha2_err[np.isfinite(val_halpha2_err)]**2 )) * 1e-18 * flux_to_lum)

    # ------------------------------------
    # Apply only region masks # NO checkerboard mask # See written notes for details
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')
    # ------------------------------------

    # Flux Errors from regions
    # --------- total --------- #
    halpha_br_err = ma.array(halpha_err, mask=bridge_mask)
    halpha_n_err = ma.array(halpha_err, mask=north_mask)
    halpha_s_err = ma.array(halpha_err, mask=south_mask)

    i0_br = np.where(halpha_br_err != -9999.0)
    val_halpha_br_err = halpha_br_err[i0_br]
    print "\n", "Total H-alpha luminosity error from bridge:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_br_err[np.isfinite(val_halpha_br_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_n = np.where(halpha_n_err != -9999.0)
    val_halpha_n_err = halpha_n_err[i0_n]
    print "Total H-alpha luminosity error from north:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_n_err[np.isfinite(val_halpha_n_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_s = np.where(halpha_s_err != -9999.0)
    val_halpha_s_err = halpha_s_err[i0_s]
    print "Total H-alpha luminosity error from south:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_s_err[np.isfinite(val_halpha_s_err)]**2 )) * 1e-18 * flux_to_lum)

    # --------- comp1 --------- #
    halpha_br1_err = ma.array(halpha_comp1_err, mask=bridge_mask)
    halpha_n1_err = ma.array(halpha_comp1_err, mask=north_mask)
    halpha_s1_err = ma.array(halpha_comp1_err, mask=south_mask)

    i0_br1 = np.where(halpha_br1_err != -9999.0)
    val_halpha_br1_err = halpha_br1_err[i0_br1]
    print "\n", "Comp1 H-alpha luminosity error from bridge:",\
     "{:.3e}".format(np.sqrt(np.sum( val_halpha_br1_err[np.isfinite(val_halpha_br1_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_n1 = np.where(halpha_n1_err != -9999.0)
    val_halpha_n1_err = halpha_n1_err[i0_n1]
    print "Comp1 H-alpha luminosity error from north:",\
     "{:.3e}".format(np.sqrt(np.sum( val_halpha_n1_err[np.isfinite(val_halpha_n1_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_s1 = np.where(halpha_s1_err != -9999.0)
    val_halpha_s1_err = halpha_s1_err[i0_s1]
    print "Comp1 H-alpha luminosity error from south:",\
     "{:.3e}".format(np.sqrt(np.sum( val_halpha_s1_err[np.isfinite(val_halpha_s1_err)]**2 )) * 1e-18 * flux_to_lum)

    # --------- comp2 --------- #
    halpha_br2_err = ma.array(halpha_comp2_err, mask=bridge_mask)
    halpha_n2_err = ma.array(halpha_comp2_err, mask=north_mask)
    halpha_s2_err = ma.array(halpha_comp2_err, mask=south_mask)

    i0_br2 = np.where(halpha_br2_err != -9999.0)
    val_halpha_br2_err = halpha_br2_err[i0_br2]
    print "\n", "Comp2 H-alpha luminosity error from bridge:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_br2_err[np.isfinite(val_halpha_br2_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_n2 = np.where(halpha_n2_err != -9999.0)
    val_halpha_n2_err = halpha_n2_err[i0_n2]
    print "Comp2 H-alpha luminosity error from north:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_n2_err[np.isfinite(val_halpha_n2_err)]**2 )) * 1e-18 * flux_to_lum)

    i0_s2 = np.where(halpha_s2_err != -9999.0)
    val_halpha_s2_err = halpha_s2_err[i0_s2]
    print "Comp2 H-alpha luminosity error from south:",\
     "{:.3e}".format(np.sqrt( np.sum(val_halpha_s2_err[np.isfinite(val_halpha_s2_err)]**2 )) * 1e-18 * flux_to_lum)

    return None

def apply_mask_to_nii(mask_to_apply):

    # Comp1
    nii_halpha_withcut_comp1_withmask = ma.array(nii_halpha_withcut_comp1, mask=mask_to_apply)
    oiii_hbeta_for_nii_withcut_comp1_withmask = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=mask_to_apply)

    # Comp2
    nii_halpha_withcut_comp2_withmask = ma.array(nii_halpha_withcut_comp2, mask=mask_to_apply)
    oiii_hbeta_for_nii_withcut_comp2_withmask = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=mask_to_apply)

    return nii_halpha_withcut_comp1_withmask, oiii_hbeta_for_nii_withcut_comp1_withmask, \
    nii_halpha_withcut_comp2_withmask, oiii_hbeta_for_nii_withcut_comp2_withmask

if __name__ == '__main__':
    """
    This code only needs the [NII] stuff. See the notes called taffy_notes.txt
    for how this code works.
    """

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # Read in stitched cube
    stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

    get_total_fluxes(stitched_cube)
    get_total_errors(stitched_cube)

    # read in masks for single and two comp fit
    all_cases = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = all_cases['COMP1_INV'].data.astype(bool)
    comp2_inv_idx = all_cases['COMP2_INV'].data.astype(bool)
    single_idx = all_cases['SINGLE_IDX'].data.astype(bool)

    # also get mask for all possible not nan spaxels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # make masks for single comp and two comp
    comp_inv_mask = np.ma.mask_or(comp1_inv_idx, comp2_inv_idx)
    onecomp_mask = np.ma.mask_or(single_idx, comp_inv_mask)
    onecomp_mask = np.ma.masked_array(onecomp_mask, mask=all_mask)
    onecomp_mask = np.logical_not(onecomp_mask)  
    # this logical not above here is to make sure that the mask has the
    # correct boolean scheme. In the numpy scheme, 0 is ok and 1 is to
    # be masked. When I initially made the masks I made them such that 
    # 0 is where the condition failed and 1 is where it was satisfied.
    # So it is the reverse of what numpy wants. Therefore, I had to 
    # use logical_not here.

    twocomp_mask = np.logical_not(onecomp_mask)
    twocomp_mask = np.ma.masked_array(twocomp_mask, mask=all_mask)

    # assign line arrays
    # Even though I only need the [NII] arrays 
    # I have to read all of these in so that I
    # can use the functions written in the 
    # bpt_plots.py code.
    # -------------- component 1 -------------- #
    halpha_comp1 = stitched_cube['HALPHA'].data[1]
    hbeta_comp1 = stitched_cube['HBETA'].data[1]
    nii6583_comp1 = stitched_cube['NII6583'].data[1]
    oiii5007_comp1 = stitched_cube['OIII5007'].data[1]
    oi6300_comp1 = stitched_cube['OI6300'].data[1]
    oi6364_comp1 = stitched_cube['OI6364'].data[1]
    sii6716_comp1 = stitched_cube['SII6716'].data[1]
    sii6731_comp1 = stitched_cube['SII6731'].data[1]

    halpha_err_comp1 = stitched_cube['HALPHA_ERR'].data[1]
    hbeta_err_comp1 = stitched_cube['HBETA_ERR'].data[1]
    nii6583_err_comp1 = stitched_cube['NII6583_ERR'].data[1]
    oiii5007_err_comp1 = stitched_cube['OIII5007_ERR'].data[1]
    oi6300_err_comp1 = stitched_cube['OI6300_ERR'].data[1]
    oi6364_err_comp1 = stitched_cube['OI6364_ERR'].data[1]
    sii6716_err_comp1 = stitched_cube['SII6716_ERR'].data[1]
    sii6731_err_comp1 = stitched_cube['SII6731_ERR'].data[1]

    # -------------- component 2 -------------- #
    halpha_comp2 = stitched_cube['HALPHA'].data[2]
    hbeta_comp2 = stitched_cube['HBETA'].data[2]
    nii6583_comp2 = stitched_cube['NII6583'].data[2]
    oiii5007_comp2 = stitched_cube['OIII5007'].data[2]
    oi6300_comp2 = stitched_cube['OI6300'].data[2]
    oi6364_comp2 = stitched_cube['OI6364'].data[2]
    sii6716_comp2 = stitched_cube['SII6716'].data[2]
    sii6731_comp2 = stitched_cube['SII6731'].data[2]

    halpha_err_comp2 = stitched_cube['HALPHA_ERR'].data[2]
    hbeta_err_comp2 = stitched_cube['HBETA_ERR'].data[2]
    nii6583_err_comp2 = stitched_cube['NII6583_ERR'].data[2]
    oiii5007_err_comp2 = stitched_cube['OIII5007_ERR'].data[2]
    oi6300_err_comp2 = stitched_cube['OI6300_ERR'].data[2]
    oi6364_err_comp2 = stitched_cube['OI6364_ERR'].data[2]
    sii6716_err_comp2 = stitched_cube['SII6716_ERR'].data[2]
    sii6731_err_comp2 = stitched_cube['SII6731_ERR'].data[2]

    # -------------------------------------------------------------------
    # Simple check to see if the compoennts add up as they should
    # Do not delete code block.
    """
    halpha = stitched_cube['HALPHA'].data[0]
    hbeta = stitched_cube['HBETA'].data[0]
    nii6583 = stitched_cube['NII6583'].data[0]
    oiii5007 = stitched_cube['OIII5007'].data[0]
    oi6300 = stitched_cube['OI6300'].data[0]
    oi6364 = stitched_cube['OI6364'].data[0]
    sii6716 = stitched_cube['SII6716'].data[0]
    sii6731 = stitched_cube['SII6731'].data[0]

    halpha_err = stitched_cube['HALPHA_ERR'].data[0]
    hbeta_err = stitched_cube['HBETA_ERR'].data[0]
    nii6583_err = stitched_cube['NII6583_ERR'].data[0]
    oiii5007_err = stitched_cube['OIII5007_ERR'].data[0]
    oi6300_err = stitched_cube['OI6300_ERR'].data[0]
    oi6364_err = stitched_cube['OI6364_ERR'].data[0]
    sii6716_err = stitched_cube['SII6716_ERR'].data[0]
    sii6731_err = stitched_cube['SII6731_ERR'].data[0]

    # add lines which are doublets
    sii = sii6716 + sii6731
    sii_err = np.sqrt((sii6716_err)**2 + (sii6731_err)**2)

    # get arrays with significance cut applied
    nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut = \
    bpt.get_arr_withsigcut(3, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err, \
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha.shape)

    for i in range(58):
        for j in range(58):

            if (halpha_comp1[i,j] == -9999.0):
                halpha_comp1[i,j] = 0.0

            if (halpha_comp2[i,j] == -9999.0):
                halpha_comp2[i,j] = 0.0

            if np.isnan(halpha[i,j]) or np.isnan(halpha_comp1[i,j]) or np.isnan(halpha_comp2[i,j]):
                continue

            if halpha[i,j] != halpha_comp1[i,j]+halpha_comp2[i,j]:
                print "DS9 pixel:", j+1, i+1
                print halpha_comp1[i,j], halpha_comp2[i,j], halpha[i,j], halpha_comp1[i,j]+halpha_comp2[i,j]

    sys.exit(0)
    """
    # -----------

    # velocity dispersions of each component
    #vdisp_line1 = hdu_vdisp[0].data[1]
    #vdisp_line2 = hdu_vdisp[0].data[2]

    # add lines which are doublets
    sii_comp1 = sii6716_comp1 + sii6731_comp1
    sii_err_comp1 = np.sqrt((sii6716_err_comp1)**2 + (sii6731_err_comp1)**2)

    sii_comp2 = sii6716_comp2 + sii6731_comp2
    sii_err_comp2 = np.sqrt((sii6716_err_comp2)**2 + (sii6731_err_comp2)**2)

    # apply sig and baseline cuts
    # sig cut for comp 1
    nii_halpha_withcut_comp1, oi_halpha_withcut_comp1, sii_halpha_withcut_comp1, \
    nii_halpha_err_withcut_comp1, oi_halpha_err_withcut_comp1, sii_halpha_err_withcut_comp1, \
    halpha_withcut_comp1, hbeta_withcut_comp1, oiii5007_withcut_comp1, oi6300_withcut_comp1, \
    nii6583_withcut_comp1, sii_withcut_comp1, \
    oiii_hbeta_for_nii_withcut_comp1, oiii_hbeta_for_oi_withcut_comp1, oiii_hbeta_for_sii_withcut_comp1, \
    oiii_hbeta_for_nii_err_withcut_comp1, oiii_hbeta_for_oi_err_withcut_comp1, oiii_hbeta_for_sii_err_withcut_comp1 = \
    bpt.get_arr_withsigcut(3.0, halpha_comp1, halpha_err_comp1, hbeta_comp1, hbeta_err_comp1, oiii5007_comp1, oiii5007_err_comp1,\
    nii6583_comp1, nii6583_err_comp1, oi6300_comp1, oi6300_err_comp1, sii_comp1, sii_err_comp1, (58,58))
    
    # sig cut for comp 2
    nii_halpha_withcut_comp2, oi_halpha_withcut_comp2, sii_halpha_withcut_comp2, \
    nii_halpha_err_withcut_comp2, oi_halpha_err_withcut_comp2, sii_halpha_err_withcut_comp2, \
    halpha_withcut_comp2, hbeta_withcut_comp2, oiii5007_withcut_comp2, oi6300_withcut_comp2, \
    nii6583_withcut_comp2, sii_withcut_comp2, \
    oiii_hbeta_for_nii_withcut_comp2, oiii_hbeta_for_oi_withcut_comp2, oiii_hbeta_for_sii_withcut_comp2, \
    oiii_hbeta_for_nii_err_withcut_comp2, oiii_hbeta_for_oi_err_withcut_comp2, oiii_hbeta_for_sii_err_withcut_comp2 = \
    bpt.get_arr_withsigcut(3.0, halpha_comp2, halpha_err_comp2, hbeta_comp2, hbeta_err_comp2, oiii5007_comp2, oiii5007_err_comp2,\
    nii6583_comp2, nii6583_err_comp2, oi6300_comp2, oi6300_err_comp2, sii_comp2, sii_err_comp2, (58,58))

    # get region mask for region defined first in ds9
    # see process to do this detailed in the comments in the bpt_plots.py code.
    # get the region masks
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # ------------- Use the saved specific region masks to circle interesting points in the BPT ------------- #
    # i.e. like the south nucleus and the north west regions
    # Read masks in
    south_nuc_mask = vcm.get_region_mask('south_galaxy_nuc_bpt')
    north_west_mask = vcm.get_region_mask('north_galaxy_west_bpt')
    north_bridge_mask = vcm.get_region_mask('north_bridge_bpt')

    # Apply checker board pattern mask to only select independent spaxels
    # Very simple and elegant solution came from here:
    # https://www.w3resource.com/python-exercises/numpy/python-numpy-exercise-10.php
    checkerboard_mask = np.zeros((58,58), dtype=np.int)
    checkerboard_mask[::2, 1::2] = 1
    checkerboard_mask[1::2, ::2] = 1

    # Combine checkerboard and two comp masks and then apply
    checkerboard_mask = checkerboard_mask.astype(bool)
    comb_mask = np.ma.mask_or(checkerboard_mask, twocomp_mask)

    # -------------------------- apply two comp mask -------------------------- #
    # comp1
    nii_halpha_withcut_comp1 = ma.array(nii_halpha_withcut_comp1, mask=comb_mask)
    oiii_hbeta_for_nii_withcut_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=comb_mask)

    # errors
    nii_halpha_err_withcut_comp1 = ma.array(nii_halpha_err_withcut_comp1, mask=comb_mask)
    oiii_hbeta_for_nii_err_withcut_comp1 = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=comb_mask)

    # comp2
    nii_halpha_withcut_comp2 = ma.array(nii_halpha_withcut_comp2, mask=comb_mask)
    oiii_hbeta_for_nii_withcut_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=comb_mask)

    # errors
    nii_halpha_err_withcut_comp2 = ma.array(nii_halpha_err_withcut_comp2, mask=comb_mask)
    oiii_hbeta_for_nii_err_withcut_comp2 = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=comb_mask)

    # -------------------------- apply region masks -------------------------- #
    # ----------------- bridge ----------------- #
    # Comp1
    nii_halpha_withcut_bridge_comp1 = ma.array(nii_halpha_withcut_comp1, mask=bridge_mask)
    oiii_hbeta_for_nii_withcut_bridge_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=bridge_mask)

    # Errors
    nii_halpha_err_withcut_bridge_comp1 = ma.array(nii_halpha_err_withcut_comp1, mask=bridge_mask)
    oiii_hbeta_for_nii_err_withcut_bridge_comp1 = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=bridge_mask)

    # Comp2
    nii_halpha_withcut_bridge_comp2 = ma.array(nii_halpha_withcut_comp2, mask=bridge_mask)
    oiii_hbeta_for_nii_withcut_bridge_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=bridge_mask)

    # Errors
    nii_halpha_err_withcut_bridge_comp2 = ma.array(nii_halpha_err_withcut_comp2, mask=bridge_mask)
    oiii_hbeta_for_nii_err_withcut_bridge_comp2 = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=bridge_mask)

    # ----------------- north ----------------- #
    # Comp1
    nii_halpha_withcut_north_comp1 = ma.array(nii_halpha_withcut_comp1, mask=north_mask)
    oiii_hbeta_for_nii_withcut_north_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=north_mask)

    # Errors
    nii_halpha_err_withcut_north_comp1 = ma.array(nii_halpha_err_withcut_comp1, mask=north_mask)
    oiii_hbeta_for_nii_err_withcut_north_comp1 = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=north_mask)

    # Comp2
    nii_halpha_withcut_north_comp2 = ma.array(nii_halpha_withcut_comp2, mask=north_mask)
    oiii_hbeta_for_nii_withcut_north_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=north_mask)

    # Errors
    nii_halpha_err_withcut_north_comp2 = ma.array(nii_halpha_err_withcut_comp2, mask=north_mask)
    oiii_hbeta_for_nii_err_withcut_north_comp2 = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=north_mask)

    # ----------------- south ----------------- #
    # Comp1
    nii_halpha_withcut_south_comp1 = ma.array(nii_halpha_withcut_comp1, mask=south_mask)
    oiii_hbeta_for_nii_withcut_south_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=south_mask)

    # Errors
    nii_halpha_err_withcut_south_comp1 = ma.array(nii_halpha_err_withcut_comp1, mask=south_mask)
    oiii_hbeta_for_nii_err_withcut_south_comp1 = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=south_mask)

    # Comp2
    nii_halpha_withcut_south_comp2 = ma.array(nii_halpha_withcut_comp2, mask=south_mask)
    oiii_hbeta_for_nii_withcut_south_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=south_mask)

    # Errors
    nii_halpha_err_withcut_south_comp2 = ma.array(nii_halpha_err_withcut_comp2, mask=south_mask)
    oiii_hbeta_for_nii_err_withcut_south_comp2 = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=south_mask)

    # ----------------- mask for h-alpha ----------------- #
    halpha_withcut_comp1 = ma.array(halpha_withcut_comp1, mask=comb_mask)
    halpha_withcut_bridge_comp1 = ma.array(halpha_withcut_comp1, mask=bridge_mask)
    halpha_withcut_north_comp1 = ma.array(halpha_withcut_comp1, mask=north_mask)
    halpha_withcut_south_comp1 = ma.array(halpha_withcut_comp1, mask=south_mask)

    halpha_withcut_comp2 = ma.array(halpha_withcut_comp2, mask=comb_mask)
    halpha_withcut_bridge_comp2 = ma.array(halpha_withcut_comp2, mask=bridge_mask)
    halpha_withcut_north_comp2 = ma.array(halpha_withcut_comp2, mask=north_mask)
    halpha_withcut_south_comp2 = ma.array(halpha_withcut_comp2, mask=south_mask)

    # apply SNuc and NW masks
    nii_halpha_withcut_southnuc_comp1, oiii_hbeta_for_nii_withcut_southnuc_comp1, \
    nii_halpha_withcut_southnuc_comp2, oiii_hbeta_for_nii_withcut_southnuc_comp2 = apply_mask_to_nii(south_nuc_mask)

    nii_halpha_withcut_nw_comp1, oiii_hbeta_for_nii_withcut_nw_comp1, \
    nii_halpha_withcut_nw_comp2, oiii_hbeta_for_nii_withcut_nw_comp2 = apply_mask_to_nii(north_west_mask)

    nii_halpha_withcut_nb_comp1, oiii_hbeta_for_nii_withcut_nb_comp1, \
    nii_halpha_withcut_nb_comp2, oiii_hbeta_for_nii_withcut_nb_comp2 = apply_mask_to_nii(north_bridge_mask)

    # print info
    print '\n', "COMPONENT 1"
    val_idx1 = np.where(halpha_withcut_comp1 != -9999.0)
    val_idx0 = np.where(halpha_comp1 != -9999.0)

    print "Sum of orig halpha:", np.nansum(halpha_comp1[val_idx0])
    print "Sum of halpha with sig cut:", np.sum(halpha_withcut_comp1[val_idx1])

    val_br1 = np.where(halpha_withcut_bridge_comp1 != -9999.0)
    val_n1 = np.where(halpha_withcut_north_comp1 != -9999.0)
    val_s1 = np.where(halpha_withcut_south_comp1 != -9999.0)

    print "Sum of halpha in bridge:", np.sum(halpha_withcut_bridge_comp1[val_br1])
    print "Sum of halpha in north:", np.sum(halpha_withcut_north_comp1[val_n1])
    print "Sum of halpha in south:", np.sum(halpha_withcut_south_comp1[val_s1])

    # --------------
    print '\n', "COMPONENT 2"
    val_idx2 = np.where(halpha_withcut_comp2 != -9999.0)
    val_idx00 = np.where(halpha_comp2 != -9999.0)

    print "Sum of orig halpha:", np.nansum(halpha_comp2[val_idx00])
    print "Sum of halpha with sig cut:", np.sum(halpha_withcut_comp2[val_idx2])

    val_br2 = np.where(halpha_withcut_bridge_comp2 != -9999.0)
    val_n2 = np.where(halpha_withcut_north_comp2 != -9999.0)
    val_s2 = np.where(halpha_withcut_south_comp2 != -9999.0)

    print "Sum of halpha in bridge:", np.sum(halpha_withcut_bridge_comp2[val_br2])
    print "Sum of halpha in north:", np.sum(halpha_withcut_north_comp2[val_n2])
    print "Sum of halpha in south:", np.sum(halpha_withcut_south_comp2[val_s2])

    # ----------------------------------------------------------------------------------- #
    # ----------------------------------------------------------------------------------- #
    print '\n', 'Component 1 (Low velocity):'
    xerr_avg_comp1, yerr_avg_comp1 = calc_lower_lim(nii_halpha_withcut_comp1, nii_halpha_withcut_bridge_comp1, nii_halpha_withcut_north_comp1, \
        nii_halpha_withcut_south_comp1, \
        oiii_hbeta_for_nii_withcut_bridge_comp1, oiii_hbeta_for_nii_withcut_north_comp1, oiii_hbeta_for_nii_withcut_south_comp1, \
        nii_halpha_err_withcut_bridge_comp1, nii_halpha_err_withcut_north_comp1, nii_halpha_err_withcut_south_comp1, \
        oiii_hbeta_for_nii_err_withcut_bridge_comp1, oiii_hbeta_for_nii_err_withcut_north_comp1, \
        oiii_hbeta_for_nii_err_withcut_south_comp1, \
        halpha_withcut_bridge_comp1, halpha_withcut_north_comp1, halpha_withcut_south_comp1)

    plot_bpt_with_hii_shaded('nii', '1', nii_halpha_withcut_bridge_comp1, nii_halpha_withcut_north_comp1, \
        nii_halpha_withcut_south_comp1, \
        oiii_hbeta_for_nii_withcut_bridge_comp1, oiii_hbeta_for_nii_withcut_north_comp1, oiii_hbeta_for_nii_withcut_south_comp1, \
        nii_halpha_err_withcut_bridge_comp1, nii_halpha_err_withcut_north_comp1, nii_halpha_err_withcut_south_comp1, \
        oiii_hbeta_for_nii_err_withcut_bridge_comp1, oiii_hbeta_for_nii_err_withcut_north_comp1, \
        oiii_hbeta_for_nii_err_withcut_south_comp1, \
        nii_halpha_withcut_southnuc_comp1, oiii_hbeta_for_nii_withcut_southnuc_comp1, \
        nii_halpha_withcut_nw_comp1, oiii_hbeta_for_nii_withcut_nw_comp1, \
        nii_halpha_withcut_nb_comp1, oiii_hbeta_for_nii_withcut_nb_comp1, \
        np.nonzero(nii_halpha_withcut_comp1), ipac_taffy_figdir, yerr_avg_comp1)

    print '\n', 'Component 2 (High velocity):'
    xerr_avg_comp2, yerr_avg_comp2 = calc_lower_lim(nii_halpha_withcut_comp2, nii_halpha_withcut_bridge_comp2, nii_halpha_withcut_north_comp2, \
        nii_halpha_withcut_south_comp2, \
        oiii_hbeta_for_nii_withcut_bridge_comp2, oiii_hbeta_for_nii_withcut_north_comp2, oiii_hbeta_for_nii_withcut_south_comp2, \
        nii_halpha_err_withcut_bridge_comp2, nii_halpha_err_withcut_north_comp2, nii_halpha_err_withcut_south_comp2, \
        oiii_hbeta_for_nii_err_withcut_bridge_comp2, oiii_hbeta_for_nii_err_withcut_north_comp2, \
        oiii_hbeta_for_nii_err_withcut_south_comp2, \
        halpha_withcut_bridge_comp2, halpha_withcut_north_comp2, halpha_withcut_south_comp2)

    plot_bpt_with_hii_shaded('nii', '2', nii_halpha_withcut_bridge_comp2, nii_halpha_withcut_north_comp2, \
        nii_halpha_withcut_south_comp2, \
        oiii_hbeta_for_nii_withcut_bridge_comp2, oiii_hbeta_for_nii_withcut_north_comp2, oiii_hbeta_for_nii_withcut_south_comp2, \
        nii_halpha_err_withcut_bridge_comp2, nii_halpha_err_withcut_north_comp2, nii_halpha_err_withcut_south_comp2, \
        oiii_hbeta_for_nii_err_withcut_bridge_comp2, oiii_hbeta_for_nii_err_withcut_north_comp2, \
        oiii_hbeta_for_nii_err_withcut_south_comp2, \
        nii_halpha_withcut_southnuc_comp2, oiii_hbeta_for_nii_withcut_southnuc_comp2, \
        nii_halpha_withcut_nw_comp2, oiii_hbeta_for_nii_withcut_nw_comp2, \
        nii_halpha_withcut_nb_comp2, oiii_hbeta_for_nii_withcut_nb_comp2, \
        np.nonzero(nii_halpha_withcut_comp2), ipac_taffy_figdir, yerr_avg_comp2)

    # Close HDU
    stitched_cube.close()

    sys.exit(0)

    # ------------ Total ----------- #
    print '\n', 'Total (Both velocity components summed):'
    nii_halpha_withcut_bridge, nii_halpha_withcut_north, nii_halpha_withcut_south, \
    oiii_hbeta_for_nii_withcut_bridge, oiii_hbeta_for_nii_withcut_north, oiii_hbeta_for_nii_withcut_south, \
    nii_halpha_err_withcut_bridge, nii_halpha_err_withcut_north, nii_halpha_err_withcut_south, \
    oiii_hbeta_for_nii_err_withcut_bridge, oiii_hbeta_for_nii_err_withcut_north, oiii_hbeta_for_nii_err_withcut_south, \
    nii_halpha_withcut, xerr_avg_total, yerr_avg_total = prep_total(stitched_cube)
    plot_bpt_with_hii_shaded('nii', 'total', nii_halpha_withcut_bridge, nii_halpha_withcut_north, nii_halpha_withcut_south, \
        oiii_hbeta_for_nii_withcut_bridge, oiii_hbeta_for_nii_withcut_north, oiii_hbeta_for_nii_withcut_south, \
        nii_halpha_err_withcut_bridge, nii_halpha_err_withcut_north, nii_halpha_err_withcut_south, \
        oiii_hbeta_for_nii_err_withcut_bridge, oiii_hbeta_for_nii_err_withcut_north, oiii_hbeta_for_nii_err_withcut_south, \
        np.nonzero(nii_halpha_withcut), ipac_taffy_figdir, yerr_avg_total)

    # Close HDU
    stitched_cube.close()

    sys.exit(0)
