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
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon

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

def make_diff_map():

    # read in indices file
    h = fits.open(savedir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    intg_flux_comp1_hdu = fits.open(savedir + 'intg_flux_cube_comp1.fits')
    intg_flux_comp2_hdu = fits.open(savedir + 'intg_flux_cube_comp2.fits')

    intg_flux_comp1 = intg_flux_comp1_hdu[0].data
    intg_flux_comp2 = intg_flux_comp2_hdu[0].data

    diffmap = intg_flux_comp2 - intg_flux_comp1

    # NaN out some spaxels that don't seem right
    # ds9 coords [x,y]
    nan_list = [[20,47], [20,48], [20,49]]

    for coord_to_nan in nan_list:
        i = coord_to_nan[1] - 1
        j = coord_to_nan[0] - 1
        print "DS9 coords:", coord_to_nan, "   Array coords:", i,j
        diffmap[i,j] = np.nan

    # Plotting
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # Save as fits file
    hdr = lzifu_hdulist['B_LINE'].header
    # Create file
    hdu = fits.PrimaryHDU(diffmap, header=hdr)
    # write
    hdu.writeto(taffy_extdir + 'intg_flux_diff_map.fits', overwrite=True)

    # ------------ also save indices file with proper header ------------ #
    indices_hdu = fits.PrimaryHDU()
    indices_hdul = fits.HDUList(indices_hdu)

    hdr['EXTNAME'] = 'COMP1_INV'
    indices_hdul.append(fits.ImageHDU(data=comp1_inv_idx, header=hdr))

    hdr['EXTNAME'] = 'COMP2_INV'
    indices_hdul.append(fits.ImageHDU(data=comp2_inv_idx, header=hdr))

    hdr['EXTNAME'] = 'SINGLE_IDX'
    indices_hdul.append(fits.ImageHDU(data=single_idx, header=hdr))

    hdr['EXTNAME'] = 'DIFFMEAN_IDX'
    indices_hdul.append(fits.ImageHDU(data=diffmean_idx, header=hdr))

    hdr['EXTNAME'] = 'DIFFSTD_IDX'
    indices_hdul.append(fits.ImageHDU(data=diffstd_idx, header=hdr))

    hdr['EXTNAME'] = 'DIFFBOTH_IDX'
    indices_hdul.append(fits.ImageHDU(data=diffboth_idx, header=hdr))
    indices_hdul.writeto(taffy_extdir + 'all_cases_indices_with_wcs.fits', overwrite=True)

    sys.exit(0)

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    colorbrewer_cm = vcm.get_colorbrewer_cm('coolwarm')

    kernel = Gaussian2DKernel(stddev=0.9)
    diffmap = convolve(diffmap, kernel, boundary='extend')

    c = ax.contour(X, Y, diffmap, transform=ax.get_transform(wcs_lzifu),\
     cmap=colorbrewer_cm, linewidths=2.0, interpolation='None')
    ax.clabel(c, inline=True, inline_spacing=2, fontsize=8, fmt='%1.1f', lw=4, ls='-')

    # add colorbar inside figure
    cbaxes = inset_axes(ax, width='30%', height='3%', loc=8, bbox_to_anchor=[0.02, 0.08, 1, 1], bbox_transform=ax.transAxes)
    #cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels), max(levels)], orientation='horizontal')
    #cb.ax.get_children()[0].set_linewidths(10.0)
    #cb.ax.set_xlabel(r'$\mathrm{Integrated\ flux [erg\, s^{-1}\, cm^{-2}\, \AA^{-1} * km\, s^{-1}]}$', fontsize=12)

    plt.show()

    return None

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    #make_diff_map()
    #sys.exit(0)

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

    # Change MPL RC params first to stop using TeX for all text
    # becasue I can't use bold text with TeX.
    mpl.rcParams['text.usetex'] = False
    # Bold text
    f = FontProperties()
    f.set_weight('bold')

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
         alpha=0.7, edgecolors='none', transform=ax.get_transform(wcs_lzifu))

    # Plot polygon to show area covered by double line profiles
    double_poly_points = [[0.4279817,23.4983800],[0.4203425,23.5050142],[0.4147988,23.5043322],[0.4134471,23.4987525],\
    [0.4116217,23.4960244],[0.4113375,23.4929686],[0.4146883,23.4922179],[0.4137292,23.4897778],[0.4109517,23.4892753],\
    [0.4088000,23.4925042],[0.4050446,23.4918222],[0.4035975,23.4887367],[0.4054367,23.4831758],[0.4100133,23.4840372],\
    [0.4114129,23.4877788],[0.4134233,23.4870962],[0.4122038,23.4826739],[0.4118607,23.4812233],[0.4167749,23.4793119],\
    [0.4162721,23.4821358],[0.4172892,23.4849342],[0.4198692,23.4850511],[0.4197337,23.4898247],[0.4246657,23.4935163],\
    [0.4266016,23.4932433],[0.4267506,23.4931750],[0.4281167,23.4949703]]
    p = Polygon(np.array(double_poly_points), edgecolor='r', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0, zorder=5)
    ax.add_patch(p)

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    # 'blues' for intg_flux
    # 'blues2yellow' for vel and vdisp
    colorbrewer_cm = vcm.get_colorbrewer_cm('blues')

    # select contour map to plot and set the variables 
    # set the variables, levels, and limits below
    con_map_type = 'intg_flux'
    con_map_comp = 'comp2'
    con_map = intg_flux_comp2

    # apply min and max limits
    # 0 to 35000 for intg_flux
    # 0 to 500 for vdisp
    # -400 to 400 for vel
    minlim = 0
    maxlim = 35000
    minidx = np.where(con_map < minlim)
    maxidx = np.where(con_map > maxlim)
    con_map[minidx] = np.nan
    con_map[maxidx] = np.nan

    # Levels taken interactively from ds9
    # uncomment as needed
    levels = np.array([2000, 3000, 6000, 12000, 15000, 20000, 25000, 30000])  # intg flux 
    #levels = np.array([-350, -250, -200, -150, -100, 0, 100, 150, 200, 250, 350])  # vel both comp
    # both velocity compoennts have the same levels to be consistent. 
    # The ds9 maps also have the same range i.e. -350 to +350 km/s
    #levels = np.array([50, 70, 90, 130, 160, 190, 230])  # vdisp 
    
    # change all nan to None to get closed contours
    # this will go wrong for velocities because 0 is a perfectly valid velocity
    # this is only good for integrated fluxes and velocity dispersion maps
    #con_map = np.nan_to_num(con_map)  # This is NOT USED ANYMORE. We are okay with open contours.

    # try smoothing the map to get smoother contours
    # define kernel
    kernel = Gaussian2DKernel(stddev=0.9)
    con_map = convolve(con_map, kernel, boundary='extend')

    c = ax.contour(X, Y, con_map, transform=ax.get_transform(wcs_lzifu),\
     levels=levels, cmap=colorbrewer_cm, linewidths=2.0, interpolation='None')
    #ax.clabel(c, inline=True, inline_spacing=2, fontsize=8, fmt='%1.1f', lw=4, ls='-')

    # add colorbar inside figure
    cbaxes = inset_axes(ax, width='30%', height='3%', loc=8, bbox_to_anchor=[0.02, 0.08, 1, 1], bbox_transform=ax.transAxes)
    cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels), max(levels)], orientation='horizontal')
    cb.ax.get_children()[0].set_linewidths(20.0)
    #cb.ax.set_xlabel(r'$\mathrm{Integrated\ flux\, [erg\, s^{-1}\, cm^{-2}\, \AA^{-1} * km\, s^{-1}]}$', fontsize=15)
    #cb.ax.set_xlabel(r'$\mathrm{Velocity\ dispersion\, [km\, s^{-1}]}$', fontsize=15)
    #cb.ax.set_xlabel(r'$\mathrm{Radial\ Velocity\, [km\, s^{-1}]}$', fontsize=15)
    # uncomment one of the above lines as required for the contour colorbar label
    # linewidths required
    # this depends on the number of levels you want
    # so if you change the number of levels then 
    # this will have to change too (by trial and error)
    # 20.0 for intg_flux
    # 23.0 for vdisp
    # 16.0 for vel

    # save the figure
    fig.savefig(taffy_extdir + 'figures_stitched_cube/' \
        + con_map_type + '_' + con_map_comp + '_contour_smooth.pdf', dpi=150, bbox_inches='tight')
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