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
        color='maroon', markersize=8, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(xarr_n[valid_indices], yarr_n[valid_indices], \
        xerr=xarr_err_n[valid_indices], yerr=yarr_err_n[valid_indices], \
        color='goldenrod', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(xarr_s[valid_indices], yarr_s[valid_indices], \
        xerr=xarr_err_s[valid_indices], yerr=yarr_err_s[valid_indices], \
        color='midnightblue', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

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
    ax.fill_between(x_arr, y_arr + 2*yerr_avg, y_arr - 2*yerr_avg, color='lightgray')

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
    ax.legend(loc=0, prop={'size':10})

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    # save and close
    fig.savefig(figdir + 'bpt_' + plottype + '_comp' + vel_comp + '_errbar_hii_shaded.eps', dpi=300, bbox_inches='tight')

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
    xerr_avg_arr = np.concatenate((xarr_err_br[nii_nonzero], xarr_err_n[nii_nonzero], xarr_err_s[nii_nonzero]))
    yerr_avg_arr = np.concatenate((yarr_err_br[nii_nonzero], yarr_err_n[nii_nonzero], yarr_err_s[nii_nonzero]))
    
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
    full_x_arr = np.concatenate((xarr_br[nii_nonzero], xarr_n[nii_nonzero], xarr_s[nii_nonzero]))
    full_y_arr = np.concatenate((yarr_br[nii_nonzero], yarr_n[nii_nonzero], yarr_s[nii_nonzero]))

    full_halpha_arr = np.concatenate((halpha_withcut_bridge[nii_nonzero], halpha_withcut_north[nii_nonzero], \
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

    return nii_halpha_withcut_bridge, nii_halpha_withcut_north, nii_halpha_withcut_south, \
    oiii_hbeta_for_nii_withcut_bridge, oiii_hbeta_for_nii_withcut_north, oiii_hbeta_for_nii_withcut_south, \
    nii_halpha_err_withcut_bridge, nii_halpha_err_withcut_north, nii_halpha_err_withcut_south, \
    oiii_hbeta_for_nii_err_withcut_bridge, oiii_hbeta_for_nii_err_withcut_north, oiii_hbeta_for_nii_err_withcut_south, \
    nii_halpha_withcut, xerr_avg_total, yerr_avg_total

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
    bpt.get_arr_withsigcut(2.5, halpha_comp1, halpha_err_comp1, hbeta_comp1, hbeta_err_comp1, oiii5007_comp1, oiii5007_err_comp1,\
    nii6583_comp1, nii6583_err_comp1, oi6300_comp1, oi6300_err_comp1, sii_comp1, sii_err_comp1, (58,58))
    
    # sig cut for comp 2
    nii_halpha_withcut_comp2, oi_halpha_withcut_comp2, sii_halpha_withcut_comp2, \
    nii_halpha_err_withcut_comp2, oi_halpha_err_withcut_comp2, sii_halpha_err_withcut_comp2, \
    halpha_withcut_comp2, hbeta_withcut_comp2, oiii5007_withcut_comp2, oi6300_withcut_comp2, \
    nii6583_withcut_comp2, sii_withcut_comp2, \
    oiii_hbeta_for_nii_withcut_comp2, oiii_hbeta_for_oi_withcut_comp2, oiii_hbeta_for_sii_withcut_comp2, \
    oiii_hbeta_for_nii_err_withcut_comp2, oiii_hbeta_for_oi_err_withcut_comp2, oiii_hbeta_for_sii_err_withcut_comp2 = \
    bpt.get_arr_withsigcut(2.5, halpha_comp2, halpha_err_comp2, hbeta_comp2, hbeta_err_comp2, oiii5007_comp2, oiii5007_err_comp2,\
    nii6583_comp2, nii6583_err_comp2, oi6300_comp2, oi6300_err_comp2, sii_comp2, sii_err_comp2, (58,58))

    # get region mask for region defined first in ds9
    # see process to do this detailed in the comments in the bpt_plots.py code.
    # get the region masks
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

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
        np.nonzero(nii_halpha_withcut_comp2), ipac_taffy_figdir, yerr_avg_comp2)

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

    sys.exit(0)