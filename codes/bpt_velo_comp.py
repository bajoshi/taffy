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

import bpt_plots as bpt
import vel_channel_map as vm

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
ipac_taffy_figdir = home + "/Desktop/ipac/taffy/figures/"

def plotbpt(plottype, vel_comp, xarr_br, yarr_br, xarr_n, yarr_n, xarr_s, yarr_s, valid_indices, figdir):
    """
    All of the BPT classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    # read in Mappings III models and overplot
    mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
    mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
    mappings_oi_halpha_v250, mappings_oi_halpha_v300,\
    mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
    mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
    mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
    mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
        mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
        mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
        mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
        mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
        mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
        mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
        mappings_sii_halpha_v250, mappings_sii_halpha_v300 = bpt.mappings_oi_nii_sii()

    if plottype == 'nii':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
        y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

        ax.plot(xarr_br[valid_indices], yarr_br[valid_indices], 'x', color='r', markersize=6, markeredgecolor='r')
        ax.plot(xarr_n[valid_indices], yarr_n[valid_indices], 'o', color='g', markersize=2, markeredgecolor='g')
        ax.plot(xarr_s[valid_indices], yarr_s[valid_indices], 'o', color='b', markersize=2, markeredgecolor='b')

        ax.plot(np.arange(-1, 0, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.plot(mappings_nii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=1.25, label='125 km/s')
        ax.plot(mappings_nii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=1.25, label='175 km/s')
        ax.plot(mappings_nii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=1.25, label='200 km/s')
        ax.plot(mappings_nii_halpha_v250, mappings_oiii_hbeta_v250, '.-', lw=1.25, label='250 km/s')

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
                                             bbox_to_anchor=(0.55, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_compbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.22, 0.3),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    elif plottype == 'oi':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[OI]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.33 + 0.73 / (np.arange(-2.5, -0.8, 0.01) + 0.59)
        y_liner_seyfert_line = 1.30 + 1.18 * np.arange(-1.1, 0, 0.01)

        ax.plot(xarr_br[valid_indices], yarr_br[valid_indices], 'x', color='r', markersize=6, markeredgecolor='r')
        ax.plot(xarr_n[valid_indices], yarr_n[valid_indices], 'o', color='g', markersize=2, markeredgecolor='g')
        ax.plot(xarr_s[valid_indices], yarr_s[valid_indices], 'o', color='b', markersize=2, markeredgecolor='b')

        ax.plot(np.arange(-2.5, -0.8, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1.1, 0, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.plot(mappings_oi_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=1.25, label='125 km/s')
        ax.plot(mappings_oi_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=1.25, label='175 km/s')
        ax.plot(mappings_oi_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=1.25, label='200 km/s')
        ax.plot(mappings_oi_halpha_v250, mappings_oiii_hbeta_v250, '.-', lw=1.25, label='250 km/s')

        ax.set_xlim(-2.0,0)
        ax.set_ylim(-1,1)

        # labels
        seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
        anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.35, 0.93),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_seyfertbox) 

        linerbox = TextArea('LINER', textprops=dict(color='k', size=16))
        anc_linerbox = AnchoredOffsetbox(loc=2, child=linerbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.8, 0.45),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_linerbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.2, 0.2),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    elif plottype == 'sii':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[SII]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.3 + 0.72 / (np.arange(-1, 0.1, 0.01) - 0.32)
        y_liner_seyfert_line = 0.76 + 1.89 * np.arange(-0.3, 1, 0.01)

        ax.plot(xarr_br[valid_indices], yarr_br[valid_indices], 'x', color='r', markersize=6, markeredgecolor='r')
        ax.plot(xarr_n[valid_indices], yarr_n[valid_indices], 'o', color='g', markersize=2, markeredgecolor='g')
        ax.plot(xarr_s[valid_indices], yarr_s[valid_indices], 'o', color='b', markersize=2, markeredgecolor='b')

        ax.plot(np.arange(-1, 0.1, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-0.3, 1, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.plot(mappings_sii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=1.25, label='125 km/s')
        ax.plot(mappings_sii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=1.25, label='175 km/s')
        ax.plot(mappings_sii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=1.25, label='200 km/s')
        ax.plot(mappings_sii_halpha_v250, mappings_oiii_hbeta_v250, '.-', lw=1.25, label='250 km/s')

        ax.set_xlim(-1,0.5)
        ax.set_ylim(-1,1)

        # labels
        seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
        anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.35, 0.93),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_seyfertbox) 

        linerbox = TextArea('LINER', textprops=dict(color='k', size=16))
        anc_linerbox = AnchoredOffsetbox(loc=2, child=linerbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.75, 0.45),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_linerbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.22, 0.3),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    ax.legend(loc=0, prop={'size':10})

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(figdir + 'bpt_' + plottype + '_comp' + vel_comp + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()
    
    return None

def bpt_range_to_spatial(xarr, yarr, xran, yran):

    xlow = xran[0]
    xhigh = xran[1]

    ylow = yran[0]
    yhigh = yran[1]

    idx = np.where((xarr >= xlow) & (xarr <= xhigh) & (yarr >= ylow) & (yarr <= yhigh))

    spatial_mask = np.zeros(xarr.shape)
    spatial_mask[idx] = 1.0

    return spatial_mask, idx

def overlay_spatial_mask_on_sdss(spatial_mask_idx):

    # read in i band SDSS image
    sdss_i, wcs_sdss = vm.get_sdss('i')

    # only need the IFU WCS here
    h, wcs_lzifu = vm.get_lzifu_products()

    # plot sdss image
    fig, ax = vm.plot_sdss_image(sdss_i, wcs_sdss)

    im = ax.scatter(spatial_mask_idx[1], spatial_mask_idx[0], s=34, c='r', marker='s',\
     alpha=0.3, edgecolors='none', transform=ax.get_transform(wcs_lzifu))
    # had to use scatter instead of using another imshow on the same axes
    # it was ignoring the transform on the second imshow

    plt.show()

    return None

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort.fits')
    hdu_vdisp = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort_VDISP.fits')

    # assign line arrays for each component and their errors
    # -------------- component 1 -------------- #
    halpha_comp1 = h['HALPHA'].data[1]
    hbeta_comp1 = h['HBETA'].data[1]
    nii6583_comp1 = h['NII6583'].data[1]
    oiii5007_comp1 = h['OIII5007'].data[1]
    oi6300_comp1 = h['OI6300'].data[1]
    oi6364_comp1 = h['OI6364'].data[1]
    sii6716_comp1 = h['SII6716'].data[1]
    sii6731_comp1 = h['SII6731'].data[1]

    halpha_err_comp1 = h['HALPHA_ERR'].data[1]
    hbeta_err_comp1 = h['HBETA_ERR'].data[1]
    nii6583_err_comp1 = h['NII6583_ERR'].data[1]
    oiii5007_err_comp1 = h['OIII5007_ERR'].data[1]
    oi6300_err_comp1 = h['OI6300_ERR'].data[1]
    oi6364_err_comp1 = h['OI6364_ERR'].data[1]
    sii6716_err_comp1 = h['SII6716_ERR'].data[1]
    sii6731_err_comp1 = h['SII6731_ERR'].data[1]

    # -------------- component 2 -------------- #
    halpha_comp2 = h['HALPHA'].data[2]
    hbeta_comp2 = h['HBETA'].data[2]
    nii6583_comp2 = h['NII6583'].data[2]
    oiii5007_comp2 = h['OIII5007'].data[2]
    oi6300_comp2 = h['OI6300'].data[2]
    oi6364_comp2 = h['OI6364'].data[2]
    sii6716_comp2 = h['SII6716'].data[2]
    sii6731_comp2 = h['SII6731'].data[2]

    halpha_err_comp2 = h['HALPHA_ERR'].data[2]
    hbeta_err_comp2 = h['HBETA_ERR'].data[2]
    nii6583_err_comp2 = h['NII6583_ERR'].data[2]
    oiii5007_err_comp2 = h['OIII5007_ERR'].data[2]
    oi6300_err_comp2 = h['OI6300_ERR'].data[2]
    oi6364_err_comp2 = h['OI6364_ERR'].data[2]
    sii6716_err_comp2 = h['SII6716_ERR'].data[2]
    sii6731_err_comp2 = h['SII6731_ERR'].data[2]

    # velocity dispersions of each component
    vdisp_line1 = hdu_vdisp[0].data[1]
    vdisp_line2 = hdu_vdisp[0].data[2]

    # add lines which are doublets
    sii_comp1 = sii6716_comp1 + sii6731_comp1
    sii_err_comp1 = np.sqrt((sii6716_err_comp1)**2 + (sii6731_err_comp1)**2)

    sii_comp2 = sii6716_comp2 + sii6731_comp2
    sii_err_comp2 = np.sqrt((sii6716_err_comp2)**2 + (sii6731_err_comp2)**2)

    # apply sig and baseline cuts
    # sig cut for comp 1
    nii_halpha_withcut_comp1, oi_halpha_withcut_comp1, sii_halpha_withcut_comp1, \
    halpha_withcut_comp1, hbeta_withcut_comp1, oiii5007_withcut_comp1, oi6300_withcut_comp1, nii6583_withcut_comp1, sii_withcut_comp1, \
    oiii_hbeta_for_nii_withcut_comp1, oiii_hbeta_for_oi_withcut_comp1, oiii_hbeta_for_sii_withcut_comp1 = \
    bpt.get_arr_withsigcut(2.5, halpha_comp1, halpha_err_comp1, hbeta_comp1, hbeta_err_comp1, oiii5007_comp1, oiii5007_err_comp1,\
    nii6583_comp1, nii6583_err_comp1, oi6300_comp1, oi6300_err_comp1, sii_comp1, sii_err_comp1, (58,58))
    
    # sig cut for comp 2
    nii_halpha_withcut_comp2, oi_halpha_withcut_comp2, sii_halpha_withcut_comp2, \
    halpha_withcut_comp2, hbeta_withcut_comp2, oiii5007_withcut_comp2, oi6300_withcut_comp2, nii6583_withcut_comp2, sii_withcut_comp2, \
    oiii_hbeta_for_nii_withcut_comp2, oiii_hbeta_for_oi_withcut_comp2, oiii_hbeta_for_sii_withcut_comp2 = \
    bpt.get_arr_withsigcut(2.5, halpha_comp2, halpha_err_comp2, hbeta_comp2, hbeta_err_comp2, oiii5007_comp2, oiii5007_err_comp2,\
    nii6583_comp2, nii6583_err_comp2, oi6300_comp2, oi6300_err_comp2, sii_comp2, sii_err_comp2, (58,58))

    # get region mask for region defined first in ds9
    # see process to do this detailed in the comments in the bpt_plots.py code.
    region_file = open(taffy_extdir + 'bridge_bpt.reg')
    region_list = np.array(region_file.readlines()[-1].split('(')[1].split(')')[0].split(','))
    region_list = region_list.astype(np.float64)
    region_file.close()

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    region_pn = pg.Polygon(pn_list)
    bridge_mask = bpt.getregionmask(region_pn, (58,58), "bridge region.")

    # north galaxy mask
    region_file = open(taffy_extdir + 'north_galaxy_bpt.reg')
    region_list = np.array(region_file.readlines()[-1].split('(')[1].split(')')[0].split(','))
    region_list = region_list.astype(np.float64)
    region_file.close()

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    region_pn = pg.Polygon(pn_list)
    north_mask = bpt.getregionmask(region_pn, (58,58), "north galaxy region.")

    # south galaxy mask
    region_file = open(taffy_extdir + 'south_galaxy_bpt.reg')
    region_list = np.array(region_file.readlines()[-1].split('(')[1].split(')')[0].split(','))
    region_list = region_list.astype(np.float64)
    region_file.close()

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    region_pn = pg.Polygon(pn_list)
    south_mask = bpt.getregionmask(region_pn, (58,58), "south galaxy region.")

    # apply mask
    # ----------------- bridge ----------------- #
    nii_halpha_withcut_bridge_comp1 = ma.array(nii_halpha_withcut_comp1, mask=bridge_mask)
    oi_halpha_withcut_bridge_comp1 = ma.array(oi_halpha_withcut_comp1, mask=bridge_mask)
    sii_halpha_withcut_bridge_comp1 = ma.array(sii_halpha_withcut_comp1, mask=bridge_mask)

    oiii_hbeta_for_nii_withcut_bridge_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=bridge_mask)
    oiii_hbeta_for_oi_withcut_bridge_comp1 = ma.array(oiii_hbeta_for_oi_withcut_comp1, mask=bridge_mask)
    oiii_hbeta_for_sii_withcut_bridge_comp1 = ma.array(oiii_hbeta_for_sii_withcut_comp1, mask=bridge_mask)

    #halpha_withcut_bridge_comp1 = ma.array(halpha_withcut_comp1, mask=bridge_mask)
    #hbeta_withcut_bridge_comp1 = ma.array(hbeta_withcut_comp1, mask=bridge_mask)
    #oiii5007_withcut_bridge_comp1 = ma.array(oiii5007_withcut_comp1, mask=bridge_mask)
    #oi6300_withcut_bridge_comp1 = ma.array(oi6300_withcut_comp1, mask=bridge_mask)
    #nii6583_withcut_bridge_comp1 = ma.array(nii6583_withcut_comp1, mask=bridge_mask)
    #sii_withcut_bridge_comp1 = ma.array(sii_withcut_comp1, mask=bridge_mask)

    nii_halpha_withcut_bridge_comp2 = ma.array(nii_halpha_withcut_comp2, mask=bridge_mask)
    oi_halpha_withcut_bridge_comp2 = ma.array(oi_halpha_withcut_comp2, mask=bridge_mask)
    sii_halpha_withcut_bridge_comp2 = ma.array(sii_halpha_withcut_comp2, mask=bridge_mask)

    oiii_hbeta_for_nii_withcut_bridge_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=bridge_mask)
    oiii_hbeta_for_oi_withcut_bridge_comp2 = ma.array(oiii_hbeta_for_oi_withcut_comp2, mask=bridge_mask)
    oiii_hbeta_for_sii_withcut_bridge_comp2 = ma.array(oiii_hbeta_for_sii_withcut_comp2, mask=bridge_mask)

    # ----------------- north galaxy ----------------- #
    nii_halpha_withcut_north_comp1 = ma.array(nii_halpha_withcut_comp1, mask=north_mask)
    oi_halpha_withcut_north_comp1 = ma.array(oi_halpha_withcut_comp1, mask=north_mask)
    sii_halpha_withcut_north_comp1 = ma.array(sii_halpha_withcut_comp1, mask=north_mask)

    oiii_hbeta_for_nii_withcut_north_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=north_mask)
    oiii_hbeta_for_oi_withcut_north_comp1 = ma.array(oiii_hbeta_for_oi_withcut_comp1, mask=north_mask)
    oiii_hbeta_for_sii_withcut_north_comp1 = ma.array(oiii_hbeta_for_sii_withcut_comp1, mask=north_mask)

    nii_halpha_withcut_north_comp2 = ma.array(nii_halpha_withcut_comp2, mask=north_mask)
    oi_halpha_withcut_north_comp2 = ma.array(oi_halpha_withcut_comp2, mask=north_mask)
    sii_halpha_withcut_north_comp2 = ma.array(sii_halpha_withcut_comp2, mask=north_mask)

    oiii_hbeta_for_nii_withcut_north_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=north_mask)
    oiii_hbeta_for_oi_withcut_north_comp2 = ma.array(oiii_hbeta_for_oi_withcut_comp2, mask=north_mask)
    oiii_hbeta_for_sii_withcut_north_comp2 = ma.array(oiii_hbeta_for_sii_withcut_comp2, mask=north_mask)

    # ----------------- south galaxy ----------------- #
    nii_halpha_withcut_south_comp1 = ma.array(nii_halpha_withcut_comp1, mask=south_mask)
    oi_halpha_withcut_south_comp1 = ma.array(oi_halpha_withcut_comp1, mask=south_mask)
    sii_halpha_withcut_south_comp1 = ma.array(sii_halpha_withcut_comp1, mask=south_mask)

    oiii_hbeta_for_nii_withcut_south_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=south_mask)
    oiii_hbeta_for_oi_withcut_south_comp1 = ma.array(oiii_hbeta_for_oi_withcut_comp1, mask=south_mask)
    oiii_hbeta_for_sii_withcut_south_comp1 = ma.array(oiii_hbeta_for_sii_withcut_comp1, mask=south_mask)

    nii_halpha_withcut_south_comp2 = ma.array(nii_halpha_withcut_comp2, mask=south_mask)
    oi_halpha_withcut_south_comp2 = ma.array(oi_halpha_withcut_comp2, mask=south_mask)
    sii_halpha_withcut_south_comp2 = ma.array(sii_halpha_withcut_comp2, mask=south_mask)

    oiii_hbeta_for_nii_withcut_south_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=south_mask)
    oiii_hbeta_for_oi_withcut_south_comp2 = ma.array(oiii_hbeta_for_oi_withcut_comp2, mask=south_mask)
    oiii_hbeta_for_sii_withcut_south_comp2 = ma.array(oiii_hbeta_for_sii_withcut_comp2, mask=south_mask)

    # spatial mapping 
    #spatial_mask, spatial_mask_idx = \
    #bpt_range_to_spatial(nii_halpha_withcut_comp2, oiii_hbeta_for_nii_withcut_bridge_comp2, [-0.3,-0.05], [-0.2,0.0])
    #plt.imshow(spatial_mask, origin='lower')
    #plt.show()
    #overlay_spatial_mask_on_sdss(spatial_mask_idx)
    #sys.exit(0)

    # plot bpt diagrams
    # BPT with [NII]
    # -------------- component 1 -------------- #
    plotbpt('nii', '1', nii_halpha_withcut_bridge_comp1, oiii_hbeta_for_nii_withcut_bridge_comp1, nii_halpha_withcut_north_comp1,\
    oiii_hbeta_for_nii_withcut_north_comp1, nii_halpha_withcut_south_comp1, oiii_hbeta_for_nii_withcut_south_comp1, np.nonzero(nii_halpha_withcut_comp1))
    # -------------- component 2 -------------- #
    plotbpt('nii', '2', nii_halpha_withcut_bridge_comp2, oiii_hbeta_for_nii_withcut_bridge_comp2, nii_halpha_withcut_north_comp2,\
    oiii_hbeta_for_nii_withcut_north_comp2, nii_halpha_withcut_south_comp2, oiii_hbeta_for_nii_withcut_south_comp2, np.nonzero(nii_halpha_withcut_comp2))

    # BPT with [OI]
    # -------------- component 1 -------------- #
    plotbpt('oi', '1', oi_halpha_withcut_bridge_comp1, oiii_hbeta_for_oi_withcut_bridge_comp1, oi_halpha_withcut_north_comp1,\
    oiii_hbeta_for_oi_withcut_north_comp1, oi_halpha_withcut_south_comp1, oiii_hbeta_for_oi_withcut_south_comp1, np.nonzero(oi_halpha_withcut_comp1))
    # -------------- component 2 -------------- #
    plotbpt('oi', '2', oi_halpha_withcut_bridge_comp2, oiii_hbeta_for_oi_withcut_bridge_comp2, oi_halpha_withcut_north_comp2,\
    oiii_hbeta_for_oi_withcut_north_comp2, oi_halpha_withcut_south_comp2, oiii_hbeta_for_oi_withcut_south_comp2, np.nonzero(oi_halpha_withcut_comp2))

    # BPT with [SII]
    # -------------- component 1 -------------- #
    plotbpt('sii', '1', sii_halpha_withcut_bridge_comp1, oiii_hbeta_for_sii_withcut_bridge_comp1, sii_halpha_withcut_north_comp1,\
    oiii_hbeta_for_sii_withcut_north_comp1, sii_halpha_withcut_south_comp1, oiii_hbeta_for_sii_withcut_south_comp1, np.nonzero(sii_halpha_withcut_comp1))
    # -------------- component 2 -------------- #
    plotbpt('sii', '2', sii_halpha_withcut_bridge_comp2, oiii_hbeta_for_sii_withcut_bridge_comp2, sii_halpha_withcut_north_comp2,\
    oiii_hbeta_for_sii_withcut_north_comp2, sii_halpha_withcut_south_comp2, oiii_hbeta_for_sii_withcut_south_comp2, np.nonzero(sii_halpha_withcut_comp2))

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)