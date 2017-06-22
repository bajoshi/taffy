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
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText
from matplotlib.colors import LinearSegmentedColormap

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

import bpt_plots as bpt

# this function to stretch the colormap by 
# shifting its midpoint is from SO user Paul H.
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # constants:
    lightspeed = 299792.458  # km/s

    # read in i band SDSS image
    sdss_i = fits.open(taffy_dir + 'SDSS/sdss_i_cutout.fits')
    wcs_sdss = WCS(sdss_i[0].header)

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort.fits')
    wcs_lzifu = WCS(h['B_LINE'].header)
    wcs_lzifu = wcs_lzifu.sub(['longitude', 'latitude'])
    # from the header it seems like the wcs is same for both red and blue channels; as it should be 

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
    gs = gridspec.GridSpec(3,3)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure(figsize=(8,8))

    # color palette from colorbrewer2.org
    colors = ['#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#084594','#08306b']
    cm = LinearSegmentedColormap.from_list('blues', colors, N=8)

    # arbitrarily choosing -400 to +400 km/s as the velo range over which to show channel maps 
    # --> and the step I want for channel maps is 50 km/s to Nyquist sample the channels since the resolution 
    # at hbeta is 98 km/s
    maptype = 'total'

    if maptype == 'total' or maptype == 'comp1':
        low_vel_lim = sys_vel - 400.0
        high_vel_lim = sys_vel + 400.0

    if maptype == 'comp2':
        low_vel_lim = sys_vel - 300.0
        high_vel_lim = sys_vel + 500.0

    vel_step = 100.0
    
    vel_range_arr = np.arange(low_vel_lim, high_vel_lim+vel_step, vel_step)

    for j in range(len(vel_range_arr)):

        # first plot the sdss image
        row, col = divmod(j, 3)  # dividing by number of columns
        ax = fig.add_subplot(gs[row, col], projection=wcs_sdss)

        # overlay the contours on sdss image
        # get the correct sized image and normalization and plot
        norm = ImageNormalize(sdss_i[0].data, interval=ZScaleInterval(), stretch=LogStretch())
        orig_cmap = mpl.cm.Greys
        shifted_cmap = shiftedColorMap(orig_cmap, midpoint=0.6, name='shifted')
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

        # make contour plot
        vel_range_low = low_vel_lim + vel_step*j
        vel_range_high = low_vel_lim + vel_step*(j+1)
        vel_range_idx = np.where((helio_vel_arr >= vel_range_low) & (helio_vel_arr <= vel_range_high))

        # avg along spectrum axis; has the same shape as only spatial dimensions # (58,58) for taffy
        if maptype == 'total':
            vel_mean_arr = np.mean(b_line_total[vel_range_idx], axis=0)
        elif maptype == 'comp1':
            vel_mean_arr = np.mean(b_line_comp1[vel_range_idx], axis=0)
        elif maptype == 'comp2':
            vel_mean_arr = np.mean(b_line_comp2[vel_range_idx], axis=0)

        vel_mean_arr = ma.array(vel_mean_arr, mask=all_mask)
        vel_mean_arr = ma.filled(vel_mean_arr, fill_value=0.0)
        vel_mean_arr = np.nan_to_num(vel_mean_arr)

        x = np.arange(vel_mean_arr.shape[1])
        y = np.arange(vel_mean_arr.shape[0])
        X, Y = np.meshgrid(x, y)

        # levels have to be defined by hand for each vel range; no other option there
        # taken after first verifying interactively with ds9
        print j, np.mean(blue_wav_arr[vel_range_idx])

        if maptype == 'total':
            if j==0: 
                levels = np.array([7,15,30,45])
            elif j==1: 
                levels = np.array([4,10,20,40,60,80])
            elif j==2: 
                levels = np.array([4,10,20,40,60,80])
            elif j==3: 
                levels = np.array([1,4,10,20,40,60,80])
            elif j==4: 
                levels = np.array([1,4,10,20,40,60,80])
            elif j==5: 
                levels = np.array([2,9,20,40,60,80])
            elif j==6: 
                levels = np.array([2,9,18,29])
            elif j==7: 
                levels = np.array([2,9,18,29])
            elif j==8: 
                levels = np.array([3,7,15,22])

        elif maptype == 'comp1' or maptype == 'comp2':
            if j==0: 
                levels = np.array([4,7,15,30])
            elif j==1: 
                levels = np.array([4,7,15,30,60])
            elif j==2: 
                levels = np.array([3,10,15,30,50])
            elif j==3: 
                levels = np.array([1,4,10,15,30])
            elif j==4: 
                levels = np.array([1,4,10,20,30])
            elif j==5: 
                levels = np.array([1,5,10,20])
            elif j==6: 
                levels = np.array([1,5,10,19])
            elif j==7: 
                levels = np.array([1,5,15])
            elif j==8: 
                levels = np.array([1,3])

        c = ax.contour(X, Y, vel_mean_arr, transform=ax.get_transform(wcs_lzifu), levels=levels, cmap=cm)
        ax.clabel(c, inline=True, inline_spacing=0, fontsize=5, fmt='%1.1f', lw=3, ls='-')

        # remove all ticks, ticklabels, and spines
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])

        ax.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="off", left="off", right="off", labelleft="off")

        # put in velocity range label
        vel_range_str = str(int(vel_range_low)) + " to " + str(int(vel_range_high))
        vel_rangebox = TextArea(vel_range_str, textprops=dict(color='k', size=9))
        anc_vel_rangebox = AnchoredOffsetbox(loc=2, child=vel_rangebox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.07, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_vel_rangebox)

    fig.savefig(taffy_dir + 'figures/vel_channel_hbeta_' + maptype + '_withclabel.eps', dpi=150, bbox_inches='tight')
    #plt.show()

    sys.exit(0)