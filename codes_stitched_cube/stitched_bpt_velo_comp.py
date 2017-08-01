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

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
ipac_taffy_figdir = home + "/Desktop/ipac/taffy/figures_stitched_cube/"
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import bpt_plots as bpt
import bpt_velo_comp as bptv
import vel_channel_map as vcm

def save_comp_masks(onecomp_mask, twocomp_mask):

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.set_title('onecomp_mask 0=ok 1=mask')
    cax1 = ax1.imshow(onecomp_mask, origin='lower', cmap='Greys')
    fig1.colorbar(cax1)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.set_title('twocomp_mask 0=ok 1=mask')
    cax2 = ax2.imshow(twocomp_mask, origin='lower', cmap='Greys')
    fig2.colorbar(cax2)

    fig1.savefig(ipac_taffy_figdir + 'onecomp_mask.eps', dpi=150, bbox_inches='tight')
    fig2.savefig(ipac_taffy_figdir + 'twocomp_mask.eps', dpi=150, bbox_inches='tight')
    plt.show()

    return None

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # Read in stitched cube
    stitched_cube = fits.open(savedir + 'stitched_cube.fits')

    # read in masks for single and two comp fit
    all_cases = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = all_cases['COMP1_INV'].data.astype(bool)
    comp2_inv_idx = all_cases['COMP2_INV'].data.astype(bool)
    single_idx = all_cases['SINGLE_IDX'].data.astype(bool)
    diffmean_idx = all_cases['DIFFMEAN_IDX'].data.astype(bool)
    diffstd_idx = all_cases['DIFFSTD_IDX'].data.astype(bool)
    diffboth_idx = all_cases['DIFFBOTH_IDX'].data.astype(bool)

    # also get mask for all possible not nan spaxels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels')

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

    #save_comp_masks(onecomp_mask, twocomp_mask)

    # assign line arrays
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
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut = \
    bpt.get_arr_withsigcut(2.5, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha.shape)

    # get the region masks
    #bridge_mask, north_mask, south_mask = bpt.getallmasks(halpha.shape)
    bridge_mask = vcm.get_region_mask('bridge_bpt')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # apply masks
    # one comp 
    nii_halpha_onecomp = ma.array(nii_halpha_withcut, mask=onecomp_mask)
    oi_halpha_onecomp = ma.array(oi_halpha_withcut, mask=onecomp_mask)
    sii_halpha_onecomp = ma.array(sii_halpha_withcut, mask=onecomp_mask)

    oiii_hbeta_for_nii_onecomp = ma.array(oiii_hbeta_for_nii_withcut, mask=onecomp_mask)
    oiii_hbeta_for_oi_onecomp = ma.array(oiii_hbeta_for_oi_withcut, mask=onecomp_mask)
    oiii_hbeta_for_sii_onecomp = ma.array(oiii_hbeta_for_sii_withcut, mask=onecomp_mask)

    # apply bridge mask
    nii_halpha_onecomp_bridge = ma.array(nii_halpha_onecomp, mask=bridge_mask)
    oi_halpha_onecomp_bridge = ma.array(oi_halpha_onecomp, mask=bridge_mask)
    sii_halpha_onecomp_bridge = ma.array(sii_halpha_onecomp, mask=bridge_mask)

    oiii_hbeta_for_nii_onecomp_bridge = ma.array(oiii_hbeta_for_nii_onecomp, mask=bridge_mask)
    oiii_hbeta_for_oi_onecomp_bridge = ma.array(oiii_hbeta_for_oi_onecomp, mask=bridge_mask)
    oiii_hbeta_for_sii_onecomp_bridge = ma.array(oiii_hbeta_for_sii_onecomp, mask=bridge_mask)

    # apply north mask
    nii_halpha_onecomp_north = ma.array(nii_halpha_onecomp, mask=north_mask)
    oi_halpha_onecomp_north = ma.array(oi_halpha_onecomp, mask=north_mask)
    sii_halpha_onecomp_north = ma.array(sii_halpha_onecomp, mask=north_mask)

    oiii_hbeta_for_nii_onecomp_north = ma.array(oiii_hbeta_for_nii_onecomp, mask=north_mask)
    oiii_hbeta_for_oi_onecomp_north = ma.array(oiii_hbeta_for_oi_onecomp, mask=north_mask)
    oiii_hbeta_for_sii_onecomp_north = ma.array(oiii_hbeta_for_sii_onecomp, mask=north_mask)

    # apply south mask
    nii_halpha_onecomp_south = ma.array(nii_halpha_onecomp, mask=south_mask)
    oi_halpha_onecomp_south = ma.array(oi_halpha_onecomp, mask=south_mask)
    sii_halpha_onecomp_south = ma.array(sii_halpha_onecomp, mask=south_mask)

    oiii_hbeta_for_nii_onecomp_south = ma.array(oiii_hbeta_for_nii_onecomp, mask=south_mask)
    oiii_hbeta_for_oi_onecomp_south = ma.array(oiii_hbeta_for_oi_onecomp, mask=south_mask)
    oiii_hbeta_for_sii_onecomp_south = ma.array(oiii_hbeta_for_sii_onecomp, mask=south_mask)

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

    # plot bpt diagrams for one comp fits
    # BPT with [NII]
    valid_idx = np.nonzero(nii_halpha_withcut)
    xarr_br = nii_halpha_onecomp_bridge
    yarr_br = oiii_hbeta_for_nii_onecomp_bridge
    xarr_n = nii_halpha_onecomp_north
    yarr_n = oiii_hbeta_for_nii_onecomp_north
    xarr_s = nii_halpha_onecomp_south
    yarr_s = oiii_hbeta_for_nii_onecomp_south
    bptv.plotbpt('nii', 'single', xarr_br, yarr_br, xarr_n, yarr_n, xarr_s, yarr_s, valid_idx, ipac_taffy_figdir)

    # BPT with [OI]
    valid_idx = np.nonzero(oi_halpha_withcut)
    xarr_br = oi_halpha_onecomp_bridge
    yarr_br = oiii_hbeta_for_oi_onecomp_bridge
    xarr_n = oi_halpha_onecomp_north
    yarr_n = oiii_hbeta_for_oi_onecomp_north
    xarr_s = oi_halpha_onecomp_south
    yarr_s = oiii_hbeta_for_oi_onecomp_south
    bptv.plotbpt('oi', 'single', xarr_br, yarr_br, xarr_n, yarr_n, xarr_s, yarr_s, valid_idx, ipac_taffy_figdir)

    # BPT with [SII]
    valid_idx = np.nonzero(sii_halpha_withcut)
    xarr_br = sii_halpha_onecomp_bridge
    yarr_br = oiii_hbeta_for_sii_onecomp_bridge
    xarr_n = sii_halpha_onecomp_north
    yarr_n = oiii_hbeta_for_sii_onecomp_north
    xarr_s = sii_halpha_onecomp_south
    yarr_s = oiii_hbeta_for_sii_onecomp_south
    bptv.plotbpt('sii', 'single', xarr_br, yarr_br, xarr_n, yarr_n, xarr_s, yarr_s, valid_idx, ipac_taffy_figdir)

    # plot bpt diagrams for two comp fits


    stitched_cube.close()
    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)


