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

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import bpt_plots as bpt

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

    # ---------------------- Now get the padding area around the HII line ---------------------- #
    # these classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
    y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

    # get non-zero indices
    nii_nonzero = np.nonzero(nii_halpha_withcut)

    # Get average values of the errors to get padding area 
    xerr_avg_arr = np.concatenate((nii_halpha_err_withcut_bridge[nii_nonzero], \
        nii_halpha_err_withcut_north[nii_nonzero], \
        nii_halpha_err_withcut_south[nii_nonzero]))

    yerr_avg_arr = np.concatenate((oiii_hbeta_for_nii_err_withcut_bridge[nii_nonzero], \
        oiii_hbeta_for_nii_err_withcut_north[nii_nonzero], \
        oiii_hbeta_for_nii_err_withcut_south[nii_nonzero]))
    
    # Nan out values = -9999.0
    invalid_idx_xerr = np.where(xerr_avg_arr == -9999.0)
    invalid_idx_yerr = np.where(yerr_avg_arr == -9999.0)

    xerr_avg_arr[invalid_idx_xerr] = np.nan
    yerr_avg_arr[invalid_idx_yerr] = np.nan

    xerr_avg = np.nanmean(xerr_avg_arr)
    yerr_avg = np.nanmean(yerr_avg_arr)

    print "Average error on log([NII]/Ha):", xerr_avg
    print "Average error on log([OIII]/Hb):", yerr_avg

    # BPT with [NII]
    # This figure plotting code is only slightly different from the 
    # one in the bpt_plots.py code.
    # I basically copy-pasted it and modified a small portion of it.
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # labels
    ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)
    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    # plot points
    ax.errorbar(nii_halpha_withcut_bridge[nii_nonzero], oiii_hbeta_for_nii_withcut_bridge[nii_nonzero], \
        xerr=nii_halpha_err_withcut_bridge[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_bridge[nii_nonzero], \
        color='maroon', markersize=8, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(nii_halpha_withcut_north[nii_nonzero], oiii_hbeta_for_nii_withcut_north[nii_nonzero], \
        xerr=nii_halpha_err_withcut_north[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_north[nii_nonzero], \
        color='goldenrod', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(nii_halpha_withcut_south[nii_nonzero], oiii_hbeta_for_nii_withcut_south[nii_nonzero], \
        xerr=nii_halpha_err_withcut_south[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_south[nii_nonzero], \
        color='midnightblue', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

    # plot classification lines
    x_arr = np.arange(-1.5, 0.0, 0.01)
    y_arr = 1.3 + 0.61 / (x_arr - 0.05)

    ax.plot(x_arr, y_arr, '-', color='k')
    ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')
    
    # Plot grey padded region for getting H-alpha SF excited fraction
    ax.fill_between(x_arr, y_arr + 2*yerr_avg, y_arr - 2*yerr_avg, color='lightgray')

    # Find total h-alpha coming from spaxels within this grey region
    # loop over all points and check where they are
    hii_x = []
    hii_y = []
    full_x_arr = np.concatenate((nii_halpha_withcut_bridge[nii_nonzero], nii_halpha_withcut_north[nii_nonzero], \
        nii_halpha_withcut_south[nii_nonzero]))
    full_y_arr = np.concatenate((oiii_hbeta_for_nii_withcut_bridge[nii_nonzero], oiii_hbeta_for_nii_withcut_north[nii_nonzero], \
        oiii_hbeta_for_nii_withcut_south[nii_nonzero]))
    for i in range(len(full_x_arr)):

        current_x = full_x_arr[i]
        current_y = full_y_arr[i]

        y_on_grayline = 1.3 + 0.61 / (current_x - 0.05)

        if current_y <= y_on_grayline + 2*yerr_avg:
            hii_x.append(current_x)
            hii_y.append(current_y)

    # ax.scatter(hii_x, hii_y, s=5, color='green', zorder=20)
    # Do not remove above line. Useful for checking method to tag 
    # HII excited spaxels.

    # plot shock models
    ax.plot(mappings_nii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s', zorder=10)

    # other plot requirements
    ax.legend(loc=0, prop={'size':10})

    ax.set_xlim(-1,0.3)
    ax.set_ylim(-1,1)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    # region labels
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

    # save
    fig.savefig(ipac_taffy_figdir + 'bpt_nii_errbar_hii_shaded.eps', dpi=300, bbox_inches='tight')

    plt.show()
    plt.clf()
    plt.cla()
    plt.close()

    sys.exit(0)