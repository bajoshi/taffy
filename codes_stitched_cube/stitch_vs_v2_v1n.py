from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = home + '/Desktop/ipac/taffy_lzifu/products_work/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
taffydir = home + '/Desktop/ipac/taffy/'
gaussfits_dir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import stitch_map as sm

speed_of_light = 299792.458  # km/s
redshift = 0.0145

def genmap():

    # ---------------- Read in Inputs ---------------- #
    # read in indices file
    h = fits.open(gaussfits_dir + 'all_cases_indices.fits')

    comp1_inv_idx = h['COMP1_INV'].data
    comp2_inv_idx = h['COMP2_INV'].data
    single_idx = h['SINGLE_IDX'].data
    diffmean_idx = h['DIFFMEAN_IDX'].data
    diffstd_idx = h['DIFFSTD_IDX'].data
    diffboth_idx = h['DIFFBOTH_IDX'].data

    # read in velocity and vdisp maps for each component
    # and also read in lzifu result for single comp fit
    one_comp = fits.open(taffy_extdir + 'Taffy_1_comp_patched.fits')
    one_comp_vel = one_comp['V'].data[1]
    # put red line fit in array
    #r_line = one_comp['R_LINE_COMP1'].data

    # Read in velocities from fitting results
    map_comp1 = np.load(gaussfits_dir + 'vel_halpha_comp1.npy')
    map_comp2 = np.load(gaussfits_dir + 'vel_halpha_comp2.npy')
    map_onecomp = np.load(gaussfits_dir + 'vel_halpha_onecomp.npy')

    # ---------------- Prep for stitching ---------------- #
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # also define line wav
    halpha_air_wav = 6562.8  # Angstroms # in air

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # get nan spaxel array from stich_map code
    nan_single_comp_arr = sm.get_nan_arr()

    # Read in both velocity component maps
    # Comp2 map already has Vs + V2 + V2b stitched
    # You only need to replace the V2b spaxels wtih V1n
    # You can simply replace those spaxels from the comp1 map
    map_vel_comp2_hdu = fits.open(gaussfits_dir + 'vel_cube_comp2.fits')
    map_vel_comp1_hdu = fits.open(gaussfits_dir + 'vel_cube_comp1.fits')

    map_vel_comp1 = map_vel_comp1_hdu[0].data
    map_vel_comp2 = map_vel_comp2_hdu[0].data

    # Create empty new map
    map_vel_comp2_new = np.zeros(map_vel_comp2.shape)

    for i in range(58):
        for j in range(58):

            map_vel_comp2_new[i,j] = map_vel_comp2[i,j]

            if diffstd_idx[i,j]:
                map_vel_comp2_new[i,j] = map_vel_comp1[i,j]

    # write out map
    # get header from lzifu output
    hdr = one_comp['CHI2'].header
    hdr['EXTNAME'] = 'vel_comp_vs_v2_v1n'
    hdu = fits.PrimaryHDU(data=map_vel_comp2_new, header=hdr)
    hdu.writeto(gaussfits_dir + 'vel_cube_vs_v2_v1n.fits', overwrite=True)

    return None

def main():
    """
    The purpose of this code is to compare the Vs+V1+V1n map
    (this is currently the top middle panel in Fig 7) with the
    Vs+V2+V1n map. This is one of the plots the referee wanted 
    to see. I think they think that this will display the rotation
    better. I'm not entirely sure.

    This code is completely based on the codes --
    1. stitche_vel_vdisp_cube.py and 2. contours_from_linefits.py.
    """

    # Generate the Vs + V2 + V1n map
    #genmap()

    # Read in the new map and make contours
    map_vel_comp2_hdu = fits.open(gaussfits_dir + 'vel_cube_comp2.fits')
    map_vel_comp2_new_hdu = fits.open(gaussfits_dir + 'vel_cube_vs_v2_v1n.fits')

    map_vel_comp2 = map_vel_comp2_hdu[0].data
    map_vel_comp2_new = map_vel_comp2_new_hdu[0].data

    # Diagnostic figures
    # Do not delete this code block. Useful for checking.
    """
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.imshow(map_vel_comp2, origin='lower', vmin=-350, vmax=350)
    ax2.imshow(map_vel_comp2_new, origin='lower', vmin=-350, vmax=350)

    plt.show()
    """

    # ------- Contours --------- # 
    # Plotting
    # Change MPL RC params first to stop using TeX for all text
    # becasue I can't use bold text with TeX.
    mpl.rcParams['text.usetex'] = False
    # Bold text
    f = FontProperties()
    f.set_weight('bold')

    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    colorbrewer_cm = vcm.get_colorbrewer_cm('blue2yellow')

    # select contour map to plot and set the variables 
    # set the variables, levels, and limits below
    con_map = map_vel_comp2_new

    # apply min and max limits
    minlim = -500
    maxlim = 500
    minidx = np.where(con_map < minlim)
    maxidx = np.where(con_map > maxlim)
    con_map[minidx] = np.nan
    con_map[maxidx] = np.nan

    levels = np.array([-350, -250, -200, -150, -100, 0, 100, 150, 200, 250, 350])  # vel both comp

    # try smoothing the map to get smoother contours
    # define kernel
    kernel = Gaussian2DKernel(stddev=0.9)
    con_map = convolve(con_map, kernel, boundary='extend')

    c = ax.contour(X, Y, con_map, transform=ax.get_transform(wcs_lzifu),\
     levels=levels, cmap=colorbrewer_cm, linewidths=2.0, interpolation='None')

    # add colorbar inside figure
    cbaxes = inset_axes(ax, width='30%', height='3%', loc=8, bbox_to_anchor=[0.02, 0.08, 1, 1], bbox_transform=ax.transAxes)
    cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels), max(levels)], orientation='horizontal')
    cb.ax.get_children()[0].set_linewidths(16.0)
    cb.ax.set_xlabel(r'$\mathrm{Radial\ Velocity\ [km\, s^{-1}]}$', fontsize=15)

    plt.show()

    return None

if __name__ == '__main__':
    main()
    sys.exit(0)