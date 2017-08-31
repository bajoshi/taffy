from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from astropy.wcs import WCS
from astropy.visualization import ManualInterval, ZScaleInterval, LogStretch, ImageNormalize

import os
import sys
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import bpt_plots as bpt
import bpt_velo_comp as bptv
import vel_channel_map as vcm

def overplot_mask(mask_idx):

    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # only need the IFU WCS here
    h, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    im = ax.scatter(mask_idx[1], mask_idx[0], s=34, c='r', marker='s',\
     alpha=0.3, edgecolors='none', transform=ax.get_transform(wcs_lzifu))
    # had to use scatter instead of using another imshow on the same axes
    # it was ignoring the transform on the second imshow

    return fig, ax

def save_npy_to_fits(fullpath, clobber=False):

    # read in arr
    npy_arr = np.load(fullpath)

    # get dir name and filename
    filename = os.path.basename(fullpath)
    filename = filename.split('.')[0]
    basedir = os.path.dirname(fullpath) + '/'

    # get shape
    shape = npy_arr.shape

    # save as fits
    hdu = fits.PrimaryHDU(data=npy_arr)
    hdu.writeto(basedir + filename + '.fits', clobber=clobber)

    return None

def get_two_comp_mask(halpha_comp1, halpha_comp2, halpha_err_comp1, halpha_err_comp2):

    # find spaxels above threshold
    halpha_thresh_lim = 5
    halpha_threshmask_comp1 = np.ma.masked_where(halpha_comp1 < halpha_thresh_lim, halpha_comp1)
    halpha_threshmask_comp2 = np.ma.masked_where(halpha_comp2 < halpha_thresh_lim, halpha_comp2)
 
    threshmask_comp1 = halpha_threshmask_comp1.mask
    threshmask_comp2 = halpha_threshmask_comp2.mask
    
    threshmask = np.ma.mask_or(threshmask_comp1, threshmask_comp2)

    # also create mask on sig cut
    # find spaxels above sig cut
    halpha_sigcut = 3
    halpha_sigmask_comp1 = np.ma.masked_where((halpha_comp1 / halpha_err_comp1) < halpha_sigcut, halpha_comp1)
    halpha_sigmask_comp2 = np.ma.masked_where((halpha_comp2 / halpha_err_comp2) < halpha_sigcut, halpha_comp2)

    sigmask_comp1 = halpha_sigmask_comp1.mask
    sigmask_comp2 = halpha_sigmask_comp2.mask

    sigmask = np.ma.mask_or(sigmask_comp1, sigmask_comp2)

    # combine the two
    threshsig_mask = np.ma.mask_or(sigmask, threshmask)

    # also combine other masks
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    bridge_mask = np.ma.make_mask(bridge_mask)
    north_mask = np.ma.make_mask(north_mask)
    south_mask = np.ma.make_mask(south_mask)

    bridge_mask = np.logical_not(bridge_mask)
    north_mask = np.logical_not(north_mask)
    south_mask = np.logical_not(south_mask)

    ns_mask = np.ma.mask_or(north_mask, south_mask)
    nsb_mask = np.ma.mask_or(ns_mask, bridge_mask)
    nsb_mask = np.logical_not(nsb_mask)

    new_mask = np.ma.mask_or(threshsig_mask, nsb_mask)  # THIS MASK IS TO BE USED WITH TWO COMP FITS ONLY!

    #plt.imshow(new_mask, origin='lower', interpolation='nearest')
    #plt.colorbar()
    
    #mask_idx = np.where(new_mask)
    #fig, ax = overplot_mask(mask_idx)
    #plt.show()

    return new_mask

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # Read in stitched cube
    stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

    # read in masks for single and two comp fit
    all_cases = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = all_cases['COMP1_INV'].data
    comp2_inv_idx = all_cases['COMP2_INV'].data
    single_idx = all_cases['SINGLE_IDX'].data
    diffmean_idx = all_cases['DIFFMEAN_IDX'].data
    diffstd_idx = all_cases['DIFFSTD_IDX'].data
    diffboth_idx = all_cases['DIFFBOTH_IDX'].data

    # make contours for all derived quantities from the line fits
    # read in saved fit params
    amp1 = np.load(savedir + 'amp_halpha_comp1.npy')
    vel1 = np.load(savedir + 'vel_halpha_comp1.npy')
    std1 = np.load(savedir + 'std_halpha_comp1.npy')

    amp2 = np.load(savedir + 'amp_halpha_comp2.npy')
    vel2 = np.load(savedir + 'vel_halpha_comp2.npy')
    std2 = np.load(savedir + 'std_halpha_comp2.npy')

    # one comp fits
    amp0 = np.load(savedir + 'amp_halpha_onecomp.npy')
    vel0 = np.load(savedir + 'vel_halpha_onecomp.npy')
    std0 = np.load(savedir + 'std_halpha_onecomp.npy')

    # read in stitched maps for intg flux, vel, and vdisp
    intg_flux_comp1_hdu = fits.open(savedir + 'intg_flux_cube_comp1.fits')
    vel_comp1_hdu = fits.open(savedir + 'vel_cube_comp1.fits')
    vdisp_comp1_hdu = fits.open(savedir + 'vdisp_cube_comp1.fits')

    intg_flux_comp2_hdu = fits.open(savedir + 'intg_flux_cube_comp2.fits')
    vel_comp2_hdu = fits.open(savedir + 'vel_cube_comp2.fits')
    vdisp_comp2_hdu = fits.open(savedir + 'vdisp_cube_comp2.fits')

    intg_flux_comp1 = intg_flux_comp1_hdu[0].data
    vel_comp1 = vel_comp1_hdu[0].data
    vdisp_comp1 = vdisp_comp1_hdu[0].data

    intg_flux_comp2 = intg_flux_comp2_hdu[0].data
    vel_comp2 = vel_comp2_hdu[0].data
    vdisp_comp2 = vdisp_comp2_hdu[0].data

    # Make figure
    fig = plt.figure(figsize=(6,6))

    # color palette from colorbrewer2.org
    cm = vcm.get_colorbrewer_cm()

    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    h, wcs_lzifu = vcm.get_lzifu_products()

    # overlay the contours on sdss image
    # get the correct sized image and normalization and plot
    norm = ImageNormalize(sdss_i[0].data, interval=ZScaleInterval(), stretch=LogStretch())
    orig_cmap = mpl.cm.Greys
    shifted_cmap = vcm.shiftedColorMap(orig_cmap, midpoint=0.6, name='shifted')
    im = ax.imshow(sdss_i[0].data, origin='lower', cmap=shifted_cmap, vmin=-0.2, vmax=6, norm=norm)
    # FYI it looks okay even without the shifted cmap but the ability to shift it is awesome.

    ax.set_autoscale_on(False)  # to stop matplotlib from changing zoom level and able actually overplot the image and contours

    lon = ax.coords[0]
    lat = ax.coords[1]

    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')

    ax.coords.frame.set_color('None')

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # Levels taken interactively from ds9
    levels = np.array([1500, 2500, 6000, 12000, 20000, 35000, 50000, 75000])

    c = ax.contour(X, Y, intg_flux_comp1, transform=ax.get_transform(wcs_lzifu), levels=levels, cmap=cm)
    ax.clabel(c, inline=True, inline_spacing=0, fontsize=5, fmt='%1.1f', lw=3, ls='-')

    plt.show()

    # close all open fits files
    stitched_cube.close()
    all_cases.close()

    intg_flux_comp1.close()
    vel_comp1.close()
    vdisp_comp1.close()

    intg_flux_comp2.close()
    vel_comp2.close()
    vdisp_comp2.close()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)