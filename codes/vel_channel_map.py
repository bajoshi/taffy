from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText
from matplotlib.colors import LinearSegmentedColormap

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

import bpt_plots as bpt

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # constants:
    lightspeed = 299792.458  # km/s

    # read in i band SDSS image
    sdss_i = fits.open(taffy_dir + 'SDSS/taffyi_sdss.fits')
    wcs_sdss = WCS(sdss_i[0].header)

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort.fits')
    wcs_lzifu = WCS(h['B_LINE'].header)  # from the header it seems like the wcs is same for both red and blue channels; as it should be 

    # assign line arrays for the total and each component
    # -------------- total -------------- #
    b_line_total = h['B_LINE'].data
    r_line_total = h['R_LINE'].data

    # -------------- component 1 -------------- #
    b_line_comp1 = h['B_LINE_COMP1'].data
    r_line_comp1 = h['R_LINE_COMP1'].data

    # -------------- component 2 -------------- #
    b_line_comp2 = h['B_LINE_COMP2'].data
    r_line_comp2 = h['R_LINE_COMP2'].data

    # create wavelength array
    # I read these data from the corresponding headers
    delt_b = 0.3  # i.e. the wav axis is sampled at 0.3A
    blue_wav_start = 4662.0
    total_blue_res_elem = 2227

    blue_wav_arr = [blue_wav_start + delt_b*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # mask elements where LZIFU gave NaNs
    region_file = open(taffy_extdir + 'all_possibly_notnan_pixels.reg')
    region_list = np.array(region_file.readlines()[-1].split('(')[1].split(')')[0].split(','))
    region_list = region_list.astype(np.float64)
    region_file.close()

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    region_pn = pg.Polygon(pn_list)
    all_mask = bpt.getregionmask(region_pn, (58,58), "region.")

    # select a line to draw contours for
    # make sure that its rest wavelength is its wavelength in air!! the IFU data was taken on the ground
    # then make sure that you are looking for the line in the appropriate wav arr 
    # e.g. hbeta in blue and halpha in red 

    # find line index in wavelength array
    redshift = 0.0145  # average z 
    sys_vel = redshift * lightspeed  # avg systemic velocity is z*c = 4350 km/s
    hbeta_air_wav = 4861.363
    hbeta_wav = hbeta_air_wav*(1+redshift)
    hbeta_idx = np.argmin(abs(blue_wav_arr - hbeta_wav))

    # convert wavelengths to heliocentric velocities
    helio_vel_arr = ((blue_wav_arr - hbeta_air_wav) / hbeta_air_wav) * lightspeed

    # create figure and axis grid
    # this has to be done before looping over each vel range
    # also defining custom color list
    gs = gridspec.GridSpec(4,5)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure(figsize=(10,8))

    colors = ['#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594','#08306b']
    cm = LinearSegmentedColormap.from_list('blues', colors, N=8)

    # arbitrarily choosing -400 to +400 km/s as the velo range over which to show channel maps 
    # --> and the step I want for channel maps is 50 km/s to Nyquist sample the channels since the resolution 
    # at hbeta is 98 km/s
    low_vel_lim = sys_vel - 450.0
    high_vel_lim = sys_vel + 450.0
    vel_step = 50.0
    
    vel_range_arr = np.arange(low_vel_lim, high_vel_lim+vel_step, vel_step)

    for j in range(len(vel_range_arr)):

        vel_range_low = low_vel_lim + vel_step*j
        vel_range_high = low_vel_lim + vel_step*(j+1)
        vel_range_idx = np.where((helio_vel_arr >= vel_range_low) & (helio_vel_arr <= vel_range_high))

        #print vel_range_low - sys_vel, vel_range_high - sys_vel, vel_range_idx

        # avg along spectrum axis; has the same shape as only spatial dimensions # (58,58) for taffy
        vel_mean_arr = np.mean(b_line_total[vel_range_idx], axis=0)
        vel_mean_arr = ma.array(vel_mean_arr, mask=all_mask)

        # make contour plot
        row, col = divmod(j, 5)  # dividing by number of columns
        ax = fig.add_subplot(gs[row, col], projection=wcs_lzifu)

        x = np.arange(vel_mean_arr.shape[0])
        y = np.arange(vel_mean_arr.shape[1])
        X, Y = np.meshgrid(x, y)
        X, Y = X.T, Y.T

        # levels have to be defined by hand for each vel range; no other option there
        # taken after first verifying with ds9
        print np.mean(blue_wav_arr[vel_range_idx])
        levels = np.array([1,3,5,8,9,12,15])

        c = ax.contour(X, Y, vel_mean_arr, levels=levels, cmap=cm)
        ax.clabel(c, inline=True, inline_spacing=2, fontsize=5, fmt='%1.2f')

        # overlay the contours on sdss image
        # get the correct sized image and normalization and plot
        ax = WCSAxes(fig, [], wcs=wcs_sdss)
        fig.add_axes(ax)
        cutout = sdss_i[0].data[45:-80,100:-50]
        im = ax.imshow(cutout, origin='lower', cmap='Greys', vmin=-0.0545, vmax=0.5)
        # [45:-80,100:-50] is hte array slice that will be zoomed in on just the galaxies 
        # but the WCS coordinates are incorrect then
        # ok in this case, because I don't need to have tick marks for the coordinates

        # remove all ticks, ticklabels, and spines
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])

        ax.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="off", left="off", right="off", labelleft="off")

    plt.show()

    sys.exit(0)