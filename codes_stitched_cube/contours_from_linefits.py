from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from astropy.wcs import WCS
from astropy.visualization import ManualInterval, ZScaleInterval, LogStretch, ImageNormalize
from astropy.convolution import convolve, Gaussian2DKernel

import os
import sys
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import bpt_plots as bpt
import bpt_velo_comp as bptv
import vel_channel_map as vcm
import stitch_map as sm

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

    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

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

    # Plotting
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    # plot mask if desired; will also have to set which mask to plot
    plot_mask = True
    if plot_mask == True:
        # get x and y coords in mask
        mask_to_plot = single_idx
        mask_x = np.where(mask_to_plot)[1]
        mask_y = np.where(mask_to_plot)[0]

        # also get the 2-comp nan pixels which one show 1-comp
        nan_single_comp_arr = sm.get_nan_arr()

        nan_x = zip(*nan_single_comp_arr)[1]
        nan_y = zip(*nan_single_comp_arr)[0]

        mask_x = np.append(mask_x, nan_x)
        mask_y = np.append(mask_y, nan_y)

        # get rid of some repeated x,y pairs
        # from SO solution: https://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
        zipped_mask = zip(mask_x, mask_y)
        mask_unique = np.vstack({tuple(row) for row in zipped_mask})

        mask_unique_x = zip(*mask_unique)[0]
        mask_unique_y = zip(*mask_unique)[1]

        # highlight spaxels based on mask
        im = ax.scatter(mask_unique_x, mask_unique_y, s=34, c='r', marker='s',\
         alpha=0.3, edgecolors='none', transform=ax.get_transform(wcs_lzifu))
        # had to use scatter instead of using another imshow on the same axes
        # it was ignoring the transform on the second imshow

        # get mask for spaxels with high standard deviation
        mask_to_plot = diffstd_idx
        mask_x = np.where(mask_to_plot)[1]
        mask_y = np.where(mask_to_plot)[0]

        im = ax.scatter(mask_x, mask_y, s=34, c='g', marker='s',\
         alpha=0.3, edgecolors='none', transform=ax.get_transform(wcs_lzifu))

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    cm = vcm.get_colorbrewer_cm()

    # Levels taken interactively from ds9
    # uncomment as needed
    #levels = np.array([2500, 6000, 12000, 20000, 35000, 50000, 75000])  # intg flux comp1
    #levels = np.array([3000, 6000, 12000, 20000, 30000, 60000])  # intg flux comp2
    #levels = np.array([-250, -200, -100, 0, 100, 150, 200, 350])  # vel comp1
    levels = np.array([-250, -200, -110, 0, 100, 140, 180, 220, 350])  # vel comp2
    #levels = np.array([70, 85, 120, 180, 240, 400, 500])  # vdisp comp 1 
    #levels = np.array([70, 90, 110, 180, 240, 400, 500])  # vdisp comp 2

    # select contour map to plot and set the variables 
    # set the variables and limtis below and also the levels above
    con_map_type = 'vel'
    con_map_comp = 'comp2'
    con_map = vel_comp2

    # apply min and max limits
    minlim = -600
    maxlim = 600
    minidx = np.where(con_map < minlim)
    maxidx = np.where(con_map > maxlim)
    con_map[minidx] = np.nan
    con_map[maxidx] = np.nan

    # change all nan to None to get closed contours
    # this will go wrong for velocities because 0 is a perfectly valid velocity
    # this is only good for integrated fluxes and velocity dispersion maps
    #con_map = np.nan_to_num(con_map)

    # try smoothing the map to get smoother contours
    # define kernel
    kernel = Gaussian2DKernel(stddev=1.0)
    con_map = convolve(con_map, kernel, boundary='extend')

    c = ax.contour(X, Y, con_map, transform=ax.get_transform(wcs_lzifu),\
     levels=levels, cmap=cm, linewidths=1.5, interpolation='None')
    ax.clabel(c, inline=True, inline_spacing=0, fontsize=5, fmt='%1.1f', lw=4, ls='-')

    # add colorbar inside figure
    cbaxes = inset_axes(ax, width='30%', height='3%', loc=3)
    plt.colorbar(c, cax=cbaxes, ticks=[min(levels), max(levels)], orientation='horizontal')

    # save the figure
    fig.savefig(taffy_extdir + 'figures_stitched_cube/' \
        + con_map_type + '_' + con_map_comp + '_contour_smooth.png', dpi=150, bbox_inches='tight')
    #plt.show()

    # close all open fits files
    intg_flux_comp1_hdu.close()
    vel_comp1_hdu.close()
    vdisp_comp1_hdu.close()

    intg_flux_comp2_hdu.close()
    vel_comp2_hdu.close()
    vdisp_comp2_hdu.close()

    h.close()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)