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
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

sys.path.append(taffydir + 'codes/')
import bpt_plots as bpt
import vel_channel_map as vcm

def plotbpt(plottype, vel_comp, xarr_br, yarr_br, xarr_n, yarr_n, xarr_s, yarr_s, \
    xarr_err_br, yarr_err_br, xarr_err_n, yarr_err_n, xarr_err_s, yarr_err_s, \
    xarr_snuc, yarr_snuc, xarr_err_snuc, yarr_err_snuc, \
    xarr_nw, yarr_nw, xarr_err_nw, yarr_err_nw, \
    xarr_nb, yarr_nb, xarr_err_nb, yarr_err_nb, \
    xarr_snucm, yarr_snucm, xarr_err_snucm, yarr_err_snucm, \
    valid_indices, figdir):
    """
    All of the BPT classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

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

    # -------------- Actual plotting -------------- #
    ax.errorbar(xarr_br[valid_indices], yarr_br[valid_indices], \
        xerr=xarr_err_br[valid_indices], yerr=yarr_err_br[valid_indices], \
        color='maroon', markersize=6, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(xarr_n[valid_indices], yarr_n[valid_indices], \
        xerr=xarr_err_n[valid_indices], yerr=yarr_err_n[valid_indices], \
        color='darkgreen', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(xarr_s[valid_indices], yarr_s[valid_indices], \
        xerr=xarr_err_s[valid_indices], yerr=yarr_err_s[valid_indices], \
        color='midnightblue', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

    ax.errorbar(xarr_snuc[valid_indices], yarr_snuc[valid_indices], \
        xerr=xarr_err_snuc[valid_indices], yerr=yarr_err_snuc[valid_indices], \
        color='midnightblue', markersize=4.5, markeredgecolor='midnightblue', fmt='d', zorder=5, capsize=0, elinewidth=0.2)
    #ax.errorbar(xarr_snucm[valid_indices], yarr_snucm[valid_indices], \
    #    xerr=xarr_err_snucm[valid_indices], yerr=yarr_err_snucm[valid_indices], \
    #    color='None', markersize=5, markeredgecolor='limegreen', fmt='o', zorder=5, capsize=0, elinewidth=0.2)
    ax.errorbar(xarr_nw[valid_indices], yarr_nw[valid_indices], \
        xerr=xarr_err_nw[valid_indices], yerr=yarr_err_nw[valid_indices], \
        color='darkgreen', markersize=6, markeredgecolor='darkgreen', fmt='+', zorder=5, capsize=0, elinewidth=0.2)
    ax.errorbar(xarr_nb[valid_indices], yarr_nb[valid_indices], \
        xerr=xarr_err_nb[valid_indices], yerr=yarr_err_nb[valid_indices], \
        color='darkorange', markersize=4, markeredgecolor='darkorange', fmt='o', zorder=5, capsize=0, elinewidth=0.2)

    # Try plotting the entire north bridge region as a single point too
    # Make sure you only consider the valid indices
    # I'm not plotting hte errors here. Too much effort and I'm not sure I'll learn anything new.
    valid_idx1 = np.where(xarr_nb[valid_indices] != -9999.0)[0]
    valid_idx2 = np.where(yarr_nb[valid_indices] != -9999.0)[0]
    valid_idx = reduce(np.intersect1d, (valid_idx1, valid_idx2))

    x_nb = np.mean(xarr_nb[valid_indices][valid_idx])
    y_nb = np.mean(yarr_nb[valid_indices][valid_idx])

    #if (type(x_nb) is np.float64) and (type(y_nb) is np.float64):
    #    ax.scatter(x_nb, y_nb, s=50, edgecolors='darkorange', facecolors='darkorange', zorder=5)

    if plottype == 'nii':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
        y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

        ax.plot(np.arange(-1, 0, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')

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
                                             bbox_to_anchor=(0.55, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_compbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.32, 0.3),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    elif plottype == 'oi':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[OI]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.33 + 0.73 / (np.arange(-2.5, -0.8, 0.01) + 0.59)
        y_liner_seyfert_line = 1.30 + 1.18 * np.arange(-1.1, 0, 0.01)

        ax.plot(np.arange(-2.5, -0.8, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1.1, 0, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.plot(mappings_oi_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s')
        ax.plot(mappings_oi_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s')
        ax.plot(mappings_oi_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s')
        ax.plot(mappings_oi_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s')
        ax.plot(mappings_oi_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s')
        ax.plot(mappings_oi_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s')

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

        ax.plot(np.arange(-1, 0.1, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-0.3, 1, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.plot(mappings_sii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s')
        ax.plot(mappings_sii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s')
        ax.plot(mappings_sii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s')
        ax.plot(mappings_sii_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s')
        ax.plot(mappings_sii_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s')
        ax.plot(mappings_sii_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s')

        ax.set_xlim(-1,0.5)
        ax.set_ylim(-1,1)

        # labels
        seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
        anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.3, 0.86),\
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

    # Text indicating panel # also make it bold
    f = FontProperties()
    f.set_weight('bold')

    if vel_comp == '1':
        ax.text(0.03, 0.97, '(a) Low velocity component', verticalalignment='top', horizontalalignment='left', \
            transform=ax.transAxes, color='k', fontproperties=f, size=14)
    elif vel_comp == '2':
        ax.text(0.03, 0.97, '(b) High velocity component', verticalalignment='top', horizontalalignment='left', \
            transform=ax.transAxes, color='k', fontproperties=f, size=14)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    fig.savefig(figdir + 'bpt_' + plottype + '_comp' + vel_comp + '_full_errbar_indep.png', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()
    
    return None

def bpt_range_to_spatial(xarr, yarr, xran, yran):
    # This funciton will be passed a masked array
    # i want all masked values to be NaN'd out before
    # using np.where() because it seems to be including
    # those as well.

    # Get mask
    mx = ma.getmask(xarr)
    my = ma.getmask(yarr)

    # Set all masked values to NaN
    xarr[mx] = np.nan
    yarr[my] = np.nan

    # now search spatially 
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
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # only need the IFU WCS here
    h, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    im = ax.scatter(spatial_mask_idx[1], spatial_mask_idx[0], s=34, c='r', marker='s',\
     alpha=0.3, edgecolors='none', transform=ax.get_transform(wcs_lzifu))
    # had to use scatter instead of using another imshow on the same axes
    # it was ignoring the transform on the second imshow

    plt.show()

    return None

def apply_mask(mask_to_apply):

    # comp1
    nii_halpha_withcut_comp1_withmask = ma.array(nii_halpha_withcut_comp1, mask=mask_to_apply)
    oi_halpha_withcut_comp1_withmask = ma.array(oi_halpha_withcut_comp1, mask=mask_to_apply)
    sii_halpha_withcut_comp1_withmask = ma.array(sii_halpha_withcut_comp1, mask=mask_to_apply)

    oiii_hbeta_for_nii_withcut_comp1_withmask = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=mask_to_apply)
    oiii_hbeta_for_oi_withcut_comp1_withmask = ma.array(oiii_hbeta_for_oi_withcut_comp1, mask=mask_to_apply)
    oiii_hbeta_for_sii_withcut_comp1_withmask = ma.array(oiii_hbeta_for_sii_withcut_comp1, mask=mask_to_apply)

    # Errors
    nii_halpha_err_withcut_comp1_withmask = ma.array(nii_halpha_err_withcut_comp1, mask=mask_to_apply)
    oi_halpha_err_withcut_comp1_withmask = ma.array(oi_halpha_err_withcut_comp1, mask=mask_to_apply)
    sii_halpha_err_withcut_comp1_withmask = ma.array(sii_halpha_err_withcut_comp1, mask=mask_to_apply)

    oiii_hbeta_for_nii_err_withcut_comp1_withmask = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=mask_to_apply)
    oiii_hbeta_for_oi_err_withcut_comp1_withmask = ma.array(oiii_hbeta_for_oi_err_withcut_comp1, mask=mask_to_apply)
    oiii_hbeta_for_sii_err_withcut_comp1_withmask = ma.array(oiii_hbeta_for_sii_err_withcut_comp1, mask=mask_to_apply)

    # comp2
    nii_halpha_withcut_comp2_withmask = ma.array(nii_halpha_withcut_comp2, mask=mask_to_apply)
    oi_halpha_withcut_comp2_withmask = ma.array(oi_halpha_withcut_comp2, mask=mask_to_apply)
    sii_halpha_withcut_comp2_withmask = ma.array(sii_halpha_withcut_comp2, mask=mask_to_apply)

    oiii_hbeta_for_nii_withcut_comp2_withmask = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=mask_to_apply)
    oiii_hbeta_for_oi_withcut_comp2_withmask = ma.array(oiii_hbeta_for_oi_withcut_comp2, mask=mask_to_apply)
    oiii_hbeta_for_sii_withcut_comp2_withmask = ma.array(oiii_hbeta_for_sii_withcut_comp2, mask=mask_to_apply)

    # Errors
    nii_halpha_err_withcut_comp2_withmask = ma.array(nii_halpha_err_withcut_comp2, mask=mask_to_apply)
    oi_halpha_err_withcut_comp2_withmask = ma.array(oi_halpha_err_withcut_comp2, mask=mask_to_apply)
    sii_halpha_err_withcut_comp2_withmask = ma.array(sii_halpha_err_withcut_comp2, mask=mask_to_apply)

    oiii_hbeta_for_nii_err_withcut_comp2_withmask = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=mask_to_apply)
    oiii_hbeta_for_oi_err_withcut_comp2_withmask = ma.array(oiii_hbeta_for_oi_err_withcut_comp2, mask=mask_to_apply)
    oiii_hbeta_for_sii_err_withcut_comp2_withmask = ma.array(oiii_hbeta_for_sii_err_withcut_comp2, mask=mask_to_apply)

    return nii_halpha_withcut_comp1_withmask, oi_halpha_withcut_comp1_withmask, sii_halpha_withcut_comp1_withmask, \
    oiii_hbeta_for_nii_withcut_comp1_withmask, oiii_hbeta_for_oi_withcut_comp1_withmask, oiii_hbeta_for_sii_withcut_comp1_withmask, \
    nii_halpha_err_withcut_comp1_withmask, oi_halpha_err_withcut_comp1_withmask, sii_halpha_err_withcut_comp1_withmask, \
    oiii_hbeta_for_nii_err_withcut_comp1_withmask, oiii_hbeta_for_oi_err_withcut_comp1_withmask, oiii_hbeta_for_sii_err_withcut_comp1_withmask, \
    nii_halpha_withcut_comp2_withmask, oi_halpha_withcut_comp2_withmask, sii_halpha_withcut_comp2_withmask, \
    oiii_hbeta_for_nii_withcut_comp2_withmask, oiii_hbeta_for_oi_withcut_comp2_withmask, oiii_hbeta_for_sii_withcut_comp2_withmask, \
    nii_halpha_err_withcut_comp2_withmask, oi_halpha_err_withcut_comp2_withmask, sii_halpha_err_withcut_comp2_withmask, \
    oiii_hbeta_for_nii_err_withcut_comp2_withmask, oiii_hbeta_for_oi_err_withcut_comp2_withmask, oiii_hbeta_for_sii_err_withcut_comp2_withmask

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in lzifu output file
    stitched = True
    if stitched:
        savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'
        h = fits.open(taffy_extdir + 'stitched_cube.fits')
        ipac_taffy_figdir = home + "/Desktop/ipac/taffy_lzifu/figures_stitched_cube/"

        # read in masks for single and two comp fit
        all_cases = fits.open(savedir + 'all_cases_indices.fits')

        comp1_inv_idx = all_cases['COMP1_INV'].data.astype(bool)
        comp2_inv_idx = all_cases['COMP2_INV'].data.astype(bool)
        single_idx = all_cases['SINGLE_IDX'].data.astype(bool)
        diffmean_idx = all_cases['DIFFMEAN_IDX'].data.astype(bool)
        diffstd_idx = all_cases['DIFFSTD_IDX'].data.astype(bool)
        diffboth_idx = all_cases['DIFFBOTH_IDX'].data.astype(bool)

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

    else:
        h = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort.fits')
        hdu_vdisp = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort_VDISP.fits')
        ipac_taffy_figdir = home + "/Desktop/ipac/taffy/figures/"

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
    halpha_withcut_comp1, hbeta_withcut_comp1, oiii5007_withcut_comp1, oi6300_withcut_comp1, nii6583_withcut_comp1, sii_withcut_comp1, \
    oiii_hbeta_for_nii_withcut_comp1, oiii_hbeta_for_oi_withcut_comp1, oiii_hbeta_for_sii_withcut_comp1, \
    oiii_hbeta_for_nii_err_withcut_comp1, oiii_hbeta_for_oi_err_withcut_comp1, oiii_hbeta_for_sii_err_withcut_comp1 = \
    bpt.get_arr_withsigcut(3.0, halpha_comp1, halpha_err_comp1, hbeta_comp1, hbeta_err_comp1, oiii5007_comp1, oiii5007_err_comp1,\
    nii6583_comp1, nii6583_err_comp1, oi6300_comp1, oi6300_err_comp1, sii_comp1, sii_err_comp1, (58,58))
    
    # sig cut for comp 2
    nii_halpha_withcut_comp2, oi_halpha_withcut_comp2, sii_halpha_withcut_comp2, \
    nii_halpha_err_withcut_comp2, oi_halpha_err_withcut_comp2, sii_halpha_err_withcut_comp2, \
    halpha_withcut_comp2, hbeta_withcut_comp2, oiii5007_withcut_comp2, oi6300_withcut_comp2, nii6583_withcut_comp2, sii_withcut_comp2, \
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
    snuc_minoraxis_mask = vcm.get_region_mask('snuc_minorax_bpt')
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

    # apply mask
    if stitched:
        # apply two comp mask
        # comp1
        nii_halpha_withcut_comp1 = ma.array(nii_halpha_withcut_comp1, mask=comb_mask)
        oi_halpha_withcut_comp1 = ma.array(oi_halpha_withcut_comp1, mask=comb_mask)
        sii_halpha_withcut_comp1 = ma.array(sii_halpha_withcut_comp1, mask=comb_mask)

        oiii_hbeta_for_nii_withcut_comp1 = ma.array(oiii_hbeta_for_nii_withcut_comp1, mask=comb_mask)
        oiii_hbeta_for_oi_withcut_comp1 = ma.array(oiii_hbeta_for_oi_withcut_comp1, mask=comb_mask)
        oiii_hbeta_for_sii_withcut_comp1 = ma.array(oiii_hbeta_for_sii_withcut_comp1, mask=comb_mask)

        # errors
        nii_halpha_err_withcut_comp1 = ma.array(nii_halpha_err_withcut_comp1, mask=comb_mask)
        oi_halpha_err_withcut_comp1 = ma.array(oi_halpha_err_withcut_comp1, mask=comb_mask)
        sii_halpha_err_withcut_comp1 = ma.array(sii_halpha_err_withcut_comp1, mask=comb_mask)

        oiii_hbeta_for_nii_err_withcut_comp1 = ma.array(oiii_hbeta_for_nii_err_withcut_comp1, mask=comb_mask)
        oiii_hbeta_for_oi_err_withcut_comp1 = ma.array(oiii_hbeta_for_oi_err_withcut_comp1, mask=comb_mask)
        oiii_hbeta_for_sii_err_withcut_comp1 = ma.array(oiii_hbeta_for_sii_err_withcut_comp1, mask=comb_mask)

        # comp2
        nii_halpha_withcut_comp2 = ma.array(nii_halpha_withcut_comp2, mask=comb_mask)
        oi_halpha_withcut_comp2 = ma.array(oi_halpha_withcut_comp2, mask=comb_mask)
        sii_halpha_withcut_comp2 = ma.array(sii_halpha_withcut_comp2, mask=comb_mask)

        oiii_hbeta_for_nii_withcut_comp2 = ma.array(oiii_hbeta_for_nii_withcut_comp2, mask=comb_mask)
        oiii_hbeta_for_oi_withcut_comp2 = ma.array(oiii_hbeta_for_oi_withcut_comp2, mask=comb_mask)
        oiii_hbeta_for_sii_withcut_comp2 = ma.array(oiii_hbeta_for_sii_withcut_comp2, mask=comb_mask)

        # errors
        nii_halpha_err_withcut_comp2 = ma.array(nii_halpha_err_withcut_comp2, mask=comb_mask)
        oi_halpha_err_withcut_comp2 = ma.array(oi_halpha_err_withcut_comp2, mask=comb_mask)
        sii_halpha_err_withcut_comp2 = ma.array(sii_halpha_err_withcut_comp2, mask=comb_mask)

        oiii_hbeta_for_nii_err_withcut_comp2 = ma.array(oiii_hbeta_for_nii_err_withcut_comp2, mask=comb_mask)
        oiii_hbeta_for_oi_err_withcut_comp2 = ma.array(oiii_hbeta_for_oi_err_withcut_comp2, mask=comb_mask)
        oiii_hbeta_for_sii_err_withcut_comp2 = ma.array(oiii_hbeta_for_sii_err_withcut_comp2, mask=comb_mask)

    # spatial mapping 
    """
    spatial_mask, spatial_mask_idx = \
    bpt_range_to_spatial(nii_halpha_withcut_north_comp1, \
        oiii_hbeta_for_nii_withcut_north_comp1, [-0.35,0.0], [-0.75,0.1])
    overlay_spatial_mask_on_sdss(spatial_mask_idx)
    sys.exit(0)
    """

    nii_halpha_withcut_north_comp1, oi_halpha_withcut_north_comp1, sii_halpha_withcut_north_comp1, \
    oiii_hbeta_for_nii_withcut_north_comp1, oiii_hbeta_for_oi_withcut_north_comp1, oiii_hbeta_for_sii_withcut_north_comp1, \
    nii_halpha_err_withcut_north_comp1, oi_halpha_err_withcut_north_comp1, sii_halpha_err_withcut_north_comp1, \
    oiii_hbeta_for_nii_err_withcut_north_comp1, oiii_hbeta_for_oi_err_withcut_north_comp1, oiii_hbeta_for_sii_err_withcut_north_comp1, \
    nii_halpha_withcut_north_comp2, oi_halpha_withcut_north_comp2, sii_halpha_withcut_north_comp2, \
    oiii_hbeta_for_nii_withcut_north_comp2, oiii_hbeta_for_oi_withcut_north_comp2, oiii_hbeta_for_sii_withcut_north_comp2, \
    nii_halpha_err_withcut_north_comp2, oi_halpha_err_withcut_north_comp2, sii_halpha_err_withcut_north_comp2, \
    oiii_hbeta_for_nii_err_withcut_north_comp2, oiii_hbeta_for_oi_err_withcut_north_comp2, oiii_hbeta_for_sii_err_withcut_north_comp2 \
    = apply_mask(north_mask)

    nii_halpha_withcut_south_comp1, oi_halpha_withcut_south_comp1, sii_halpha_withcut_south_comp1, \
    oiii_hbeta_for_nii_withcut_south_comp1, oiii_hbeta_for_oi_withcut_south_comp1, oiii_hbeta_for_sii_withcut_south_comp1, \
    nii_halpha_err_withcut_south_comp1, oi_halpha_err_withcut_south_comp1, sii_halpha_err_withcut_south_comp1, \
    oiii_hbeta_for_nii_err_withcut_south_comp1, oiii_hbeta_for_oi_err_withcut_south_comp1, oiii_hbeta_for_sii_err_withcut_south_comp1, \
    nii_halpha_withcut_south_comp2, oi_halpha_withcut_south_comp2, sii_halpha_withcut_south_comp2, \
    oiii_hbeta_for_nii_withcut_south_comp2, oiii_hbeta_for_oi_withcut_south_comp2, oiii_hbeta_for_sii_withcut_south_comp2, \
    nii_halpha_err_withcut_south_comp2, oi_halpha_err_withcut_south_comp2, sii_halpha_err_withcut_south_comp2, \
    oiii_hbeta_for_nii_err_withcut_south_comp2, oiii_hbeta_for_oi_err_withcut_south_comp2, oiii_hbeta_for_sii_err_withcut_south_comp2 \
    = apply_mask(south_mask)

    nii_halpha_withcut_bridge_comp1, oi_halpha_withcut_bridge_comp1, sii_halpha_withcut_bridge_comp1, \
    oiii_hbeta_for_nii_withcut_bridge_comp1, oiii_hbeta_for_oi_withcut_bridge_comp1, oiii_hbeta_for_sii_withcut_bridge_comp1, \
    nii_halpha_err_withcut_bridge_comp1, oi_halpha_err_withcut_bridge_comp1, sii_halpha_err_withcut_bridge_comp1, \
    oiii_hbeta_for_nii_err_withcut_bridge_comp1, oiii_hbeta_for_oi_err_withcut_bridge_comp1, oiii_hbeta_for_sii_err_withcut_bridge_comp1, \
    nii_halpha_withcut_bridge_comp2, oi_halpha_withcut_bridge_comp2, sii_halpha_withcut_bridge_comp2, \
    oiii_hbeta_for_nii_withcut_bridge_comp2, oiii_hbeta_for_oi_withcut_bridge_comp2, oiii_hbeta_for_sii_withcut_bridge_comp2, \
    nii_halpha_err_withcut_bridge_comp2, oi_halpha_err_withcut_bridge_comp2, sii_halpha_err_withcut_bridge_comp2, \
    oiii_hbeta_for_nii_err_withcut_bridge_comp2, oiii_hbeta_for_oi_err_withcut_bridge_comp2, oiii_hbeta_for_sii_err_withcut_bridge_comp2 \
    = apply_mask(bridge_mask)

    nii_halpha_withcut_southnuc_comp1, oi_halpha_withcut_southnuc_comp1, sii_halpha_withcut_southnuc_comp1, \
    oiii_hbeta_for_nii_withcut_southnuc_comp1, oiii_hbeta_for_oi_withcut_southnuc_comp1, oiii_hbeta_for_sii_withcut_southnuc_comp1, \
    nii_halpha_err_withcut_southnuc_comp1, oi_halpha_err_withcut_southnuc_comp1, sii_halpha_err_withcut_southnuc_comp1, \
    oiii_hbeta_for_nii_err_withcut_southnuc_comp1, oiii_hbeta_for_oi_err_withcut_southnuc_comp1, oiii_hbeta_for_sii_err_withcut_southnuc_comp1, \
    nii_halpha_withcut_southnuc_comp2, oi_halpha_withcut_southnuc_comp2, sii_halpha_withcut_southnuc_comp2, \
    oiii_hbeta_for_nii_withcut_southnuc_comp2, oiii_hbeta_for_oi_withcut_southnuc_comp2, oiii_hbeta_for_sii_withcut_southnuc_comp2, \
    nii_halpha_err_withcut_southnuc_comp2, oi_halpha_err_withcut_southnuc_comp2, sii_halpha_err_withcut_southnuc_comp2, \
    oiii_hbeta_for_nii_err_withcut_southnuc_comp2, oiii_hbeta_for_oi_err_withcut_southnuc_comp2, oiii_hbeta_for_sii_err_withcut_southnuc_comp2 \
    = apply_mask(south_nuc_mask)

    nii_halpha_withcut_southnucm_comp1, oi_halpha_withcut_southnucm_comp1, sii_halpha_withcut_southnucm_comp1, \
    oiii_hbeta_for_nii_withcut_southnucm_comp1, oiii_hbeta_for_oi_withcut_southnucm_comp1, oiii_hbeta_for_sii_withcut_southnucm_comp1, \
    nii_halpha_err_withcut_southnucm_comp1, oi_halpha_err_withcut_southnucm_comp1, sii_halpha_err_withcut_southnucm_comp1, \
    oiii_hbeta_for_nii_err_withcut_southnucm_comp1, oiii_hbeta_for_oi_err_withcut_southnucm_comp1, oiii_hbeta_for_sii_err_withcut_southnucm_comp1, \
    nii_halpha_withcut_southnucm_comp2, oi_halpha_withcut_southnucm_comp2, sii_halpha_withcut_southnucm_comp2, \
    oiii_hbeta_for_nii_withcut_southnucm_comp2, oiii_hbeta_for_oi_withcut_southnucm_comp2, oiii_hbeta_for_sii_withcut_southnucm_comp2, \
    nii_halpha_err_withcut_southnucm_comp2, oi_halpha_err_withcut_southnucm_comp2, sii_halpha_err_withcut_southnucm_comp2, \
    oiii_hbeta_for_nii_err_withcut_southnucm_comp2, oiii_hbeta_for_oi_err_withcut_southnucm_comp2, oiii_hbeta_for_sii_err_withcut_southnucm_comp2 \
    = apply_mask(snuc_minoraxis_mask)

    nii_halpha_withcut_nw_comp1, oi_halpha_withcut_nw_comp1, sii_halpha_withcut_nw_comp1, \
    oiii_hbeta_for_nii_withcut_nw_comp1, oiii_hbeta_for_oi_withcut_nw_comp1, oiii_hbeta_for_sii_withcut_nw_comp1, \
    nii_halpha_err_withcut_nw_comp1, oi_halpha_err_withcut_nw_comp1, sii_halpha_err_withcut_nw_comp1, \
    oiii_hbeta_for_nii_err_withcut_nw_comp1, oiii_hbeta_for_oi_err_withcut_nw_comp1, oiii_hbeta_for_sii_err_withcut_nw_comp1, \
    nii_halpha_withcut_nw_comp2, oi_halpha_withcut_nw_comp2, sii_halpha_withcut_nw_comp2, \
    oiii_hbeta_for_nii_withcut_nw_comp2, oiii_hbeta_for_oi_withcut_nw_comp2, oiii_hbeta_for_sii_withcut_nw_comp2, \
    nii_halpha_err_withcut_nw_comp2, oi_halpha_err_withcut_nw_comp2, sii_halpha_err_withcut_nw_comp2, \
    oiii_hbeta_for_nii_err_withcut_nw_comp2, oiii_hbeta_for_oi_err_withcut_nw_comp2, oiii_hbeta_for_sii_err_withcut_nw_comp2 \
    = apply_mask(north_west_mask)

    nii_halpha_withcut_nb_comp1, oi_halpha_withcut_nb_comp1, sii_halpha_withcut_nb_comp1, \
    oiii_hbeta_for_nii_withcut_nb_comp1, oiii_hbeta_for_oi_withcut_nb_comp1, oiii_hbeta_for_sii_withcut_nb_comp1, \
    nii_halpha_err_withcut_nb_comp1, oi_halpha_err_withcut_nb_comp1, sii_halpha_err_withcut_nb_comp1, \
    oiii_hbeta_for_nii_err_withcut_nb_comp1, oiii_hbeta_for_oi_err_withcut_nb_comp1, oiii_hbeta_for_sii_err_withcut_nb_comp1, \
    nii_halpha_withcut_nb_comp2, oi_halpha_withcut_nb_comp2, sii_halpha_withcut_nb_comp2, \
    oiii_hbeta_for_nii_withcut_nb_comp2, oiii_hbeta_for_oi_withcut_nb_comp2, oiii_hbeta_for_sii_withcut_nb_comp2, \
    nii_halpha_err_withcut_nb_comp2, oi_halpha_err_withcut_nb_comp2, sii_halpha_err_withcut_nb_comp2, \
    oiii_hbeta_for_nii_err_withcut_nb_comp2, oiii_hbeta_for_oi_err_withcut_nb_comp2, oiii_hbeta_for_sii_err_withcut_nb_comp2 \
    = apply_mask(north_bridge_mask)

    # plot bpt diagrams
    # BPT with [NII]
    # -------------- component 1 -------------- #
    plotbpt('nii', '1', nii_halpha_withcut_bridge_comp1, oiii_hbeta_for_nii_withcut_bridge_comp1, nii_halpha_withcut_north_comp1, \
    oiii_hbeta_for_nii_withcut_north_comp1, nii_halpha_withcut_south_comp1, oiii_hbeta_for_nii_withcut_south_comp1, \
    nii_halpha_err_withcut_bridge_comp1, oiii_hbeta_for_nii_err_withcut_bridge_comp1, nii_halpha_err_withcut_north_comp1, \
    oiii_hbeta_for_nii_err_withcut_north_comp1, nii_halpha_err_withcut_south_comp1, oiii_hbeta_for_nii_err_withcut_south_comp1, \
    nii_halpha_withcut_southnuc_comp1, oiii_hbeta_for_nii_withcut_southnuc_comp1, \
    nii_halpha_err_withcut_southnuc_comp1, oiii_hbeta_for_nii_err_withcut_southnuc_comp1, \
    nii_halpha_withcut_nw_comp1, oiii_hbeta_for_nii_withcut_nw_comp1,
    nii_halpha_err_withcut_nw_comp1, oiii_hbeta_for_nii_err_withcut_nw_comp1,
    nii_halpha_withcut_nb_comp1, oiii_hbeta_for_nii_withcut_nb_comp1,
    nii_halpha_err_withcut_nb_comp1, oiii_hbeta_for_nii_err_withcut_nb_comp1,
    nii_halpha_withcut_southnucm_comp1, oiii_hbeta_for_nii_withcut_southnucm_comp1,
    nii_halpha_err_withcut_southnucm_comp1, oiii_hbeta_for_nii_err_withcut_southnucm_comp1,
    np.nonzero(nii_halpha_withcut_comp1), ipac_taffy_figdir)
    # -------------- component 2 -------------- #
    plotbpt('nii', '2', nii_halpha_withcut_bridge_comp2, oiii_hbeta_for_nii_withcut_bridge_comp2, nii_halpha_withcut_north_comp2, \
    oiii_hbeta_for_nii_withcut_north_comp2, nii_halpha_withcut_south_comp2, oiii_hbeta_for_nii_withcut_south_comp2, \
    nii_halpha_err_withcut_bridge_comp2, oiii_hbeta_for_nii_err_withcut_bridge_comp2, nii_halpha_err_withcut_north_comp2, \
    oiii_hbeta_for_nii_err_withcut_north_comp2, nii_halpha_err_withcut_south_comp2, oiii_hbeta_for_nii_err_withcut_south_comp2, \
    nii_halpha_withcut_southnuc_comp2, oiii_hbeta_for_nii_withcut_southnuc_comp2, \
    nii_halpha_err_withcut_southnuc_comp2, oiii_hbeta_for_nii_err_withcut_southnuc_comp2, \
    nii_halpha_withcut_nw_comp2, oiii_hbeta_for_nii_withcut_nw_comp2,
    nii_halpha_err_withcut_nw_comp2, oiii_hbeta_for_nii_err_withcut_nw_comp2,
    nii_halpha_withcut_nb_comp2, oiii_hbeta_for_nii_withcut_nb_comp2,
    nii_halpha_err_withcut_nb_comp2, oiii_hbeta_for_nii_err_withcut_nb_comp2,
    nii_halpha_withcut_southnucm_comp2, oiii_hbeta_for_nii_withcut_southnucm_comp2,
    nii_halpha_err_withcut_southnucm_comp2, oiii_hbeta_for_nii_err_withcut_southnucm_comp2,
    np.nonzero(nii_halpha_withcut_comp2), ipac_taffy_figdir)

    # BPT with [OI]
    # -------------- component 1 -------------- #
    plotbpt('oi', '1', oi_halpha_withcut_bridge_comp1, oiii_hbeta_for_oi_withcut_bridge_comp1, oi_halpha_withcut_north_comp1, \
    oiii_hbeta_for_oi_withcut_north_comp1, oi_halpha_withcut_south_comp1, oiii_hbeta_for_oi_withcut_south_comp1, \
    oi_halpha_err_withcut_bridge_comp1, oiii_hbeta_for_oi_err_withcut_bridge_comp1, oi_halpha_err_withcut_north_comp1, \
    oiii_hbeta_for_oi_err_withcut_north_comp1, oi_halpha_err_withcut_south_comp1, oiii_hbeta_for_oi_err_withcut_south_comp1, \
    oi_halpha_withcut_southnuc_comp1, oiii_hbeta_for_oi_withcut_southnuc_comp1, \
    oi_halpha_err_withcut_southnuc_comp1, oiii_hbeta_for_oi_err_withcut_southnuc_comp1, \
    oi_halpha_withcut_nw_comp1, oiii_hbeta_for_oi_withcut_nw_comp1,
    oi_halpha_err_withcut_nw_comp1, oiii_hbeta_for_oi_err_withcut_nw_comp1,
    oi_halpha_withcut_nb_comp1, oiii_hbeta_for_oi_withcut_nb_comp1,
    oi_halpha_err_withcut_nb_comp1, oiii_hbeta_for_oi_err_withcut_nb_comp1,
    oi_halpha_withcut_southnucm_comp1, oiii_hbeta_for_oi_withcut_southnucm_comp1,
    oi_halpha_err_withcut_southnucm_comp1, oiii_hbeta_for_oi_err_withcut_southnucm_comp1,
    np.nonzero(oi_halpha_withcut_comp1), ipac_taffy_figdir)
    # -------------- component 2 -------------- #
    plotbpt('oi', '2', oi_halpha_withcut_bridge_comp2, oiii_hbeta_for_oi_withcut_bridge_comp2, oi_halpha_withcut_north_comp2, \
    oiii_hbeta_for_oi_withcut_north_comp2, oi_halpha_withcut_south_comp2, oiii_hbeta_for_oi_withcut_south_comp2, \
    oi_halpha_err_withcut_bridge_comp2, oiii_hbeta_for_oi_err_withcut_bridge_comp2, oi_halpha_err_withcut_north_comp2, \
    oiii_hbeta_for_oi_err_withcut_north_comp2, oi_halpha_err_withcut_south_comp2, oiii_hbeta_for_oi_err_withcut_south_comp2, \
    oi_halpha_withcut_southnuc_comp2, oiii_hbeta_for_oi_withcut_southnuc_comp2, \
    oi_halpha_err_withcut_southnuc_comp2, oiii_hbeta_for_oi_err_withcut_southnuc_comp2, \
    oi_halpha_withcut_nw_comp2, oiii_hbeta_for_oi_withcut_nw_comp2,
    oi_halpha_err_withcut_nw_comp2, oiii_hbeta_for_oi_err_withcut_nw_comp2,
    oi_halpha_withcut_nb_comp2, oiii_hbeta_for_oi_withcut_nb_comp2,
    oi_halpha_err_withcut_nb_comp2, oiii_hbeta_for_oi_err_withcut_nb_comp2,
    oi_halpha_withcut_southnucm_comp2, oiii_hbeta_for_oi_withcut_southnucm_comp2,
    oi_halpha_err_withcut_southnucm_comp2, oiii_hbeta_for_oi_err_withcut_southnucm_comp2,
    np.nonzero(oi_halpha_withcut_comp2), ipac_taffy_figdir)

    # BPT with [SII]
    # -------------- component 1 -------------- #
    plotbpt('sii', '1', sii_halpha_withcut_bridge_comp1, oiii_hbeta_for_sii_withcut_bridge_comp1, sii_halpha_withcut_north_comp1, \
    oiii_hbeta_for_sii_withcut_north_comp1, sii_halpha_withcut_south_comp1, oiii_hbeta_for_sii_withcut_south_comp1, \
    sii_halpha_err_withcut_bridge_comp1, oiii_hbeta_for_sii_err_withcut_bridge_comp1, sii_halpha_err_withcut_north_comp1, \
    oiii_hbeta_for_sii_err_withcut_north_comp1, sii_halpha_err_withcut_south_comp1, oiii_hbeta_for_sii_err_withcut_south_comp1, \
    sii_halpha_withcut_southnuc_comp1, oiii_hbeta_for_sii_withcut_southnuc_comp1, \
    sii_halpha_err_withcut_southnuc_comp1, oiii_hbeta_for_sii_err_withcut_southnuc_comp1, \
    sii_halpha_withcut_nw_comp1, oiii_hbeta_for_sii_withcut_nw_comp1,
    sii_halpha_err_withcut_nw_comp1, oiii_hbeta_for_sii_err_withcut_nw_comp1,
    sii_halpha_withcut_nb_comp1, oiii_hbeta_for_sii_withcut_nb_comp1,
    sii_halpha_err_withcut_nb_comp1, oiii_hbeta_for_sii_err_withcut_nb_comp1,
    sii_halpha_withcut_southnucm_comp1, oiii_hbeta_for_sii_withcut_southnucm_comp1,
    sii_halpha_err_withcut_southnucm_comp1, oiii_hbeta_for_sii_err_withcut_southnucm_comp1,
    np.nonzero(sii_halpha_withcut_comp1), ipac_taffy_figdir)
    # -------------- component 2 -------------- #
    plotbpt('sii', '2', sii_halpha_withcut_bridge_comp2, oiii_hbeta_for_sii_withcut_bridge_comp2, sii_halpha_withcut_north_comp2, \
    oiii_hbeta_for_sii_withcut_north_comp2, sii_halpha_withcut_south_comp2, oiii_hbeta_for_sii_withcut_south_comp2, \
    sii_halpha_err_withcut_bridge_comp2, oiii_hbeta_for_sii_err_withcut_bridge_comp2, sii_halpha_err_withcut_north_comp2, \
    oiii_hbeta_for_sii_err_withcut_north_comp2, sii_halpha_err_withcut_south_comp2, oiii_hbeta_for_sii_err_withcut_south_comp2, \
    sii_halpha_withcut_southnuc_comp2, oiii_hbeta_for_sii_withcut_southnuc_comp2, \
    sii_halpha_err_withcut_southnuc_comp2, oiii_hbeta_for_sii_err_withcut_southnuc_comp2, \
    sii_halpha_withcut_nw_comp2, oiii_hbeta_for_sii_withcut_nw_comp2,
    sii_halpha_err_withcut_nw_comp2, oiii_hbeta_for_sii_err_withcut_nw_comp2,
    sii_halpha_withcut_nb_comp2, oiii_hbeta_for_sii_withcut_nb_comp2,
    sii_halpha_err_withcut_nb_comp2, oiii_hbeta_for_sii_err_withcut_nb_comp2,
    sii_halpha_withcut_southnucm_comp2, oiii_hbeta_for_sii_withcut_southnucm_comp2,
    sii_halpha_err_withcut_southnucm_comp2, oiii_hbeta_for_sii_err_withcut_southnucm_comp2,
    np.nonzero(sii_halpha_withcut_comp2), ipac_taffy_figdir)

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)