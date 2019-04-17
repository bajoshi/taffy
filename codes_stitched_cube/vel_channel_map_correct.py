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
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

sys.path.append(taffy_dir + 'codes/')
import vel_channel_map as vcm

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # constants:
    lightspeed = 299792.458  # km/s

    # Plotting
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # assign array for the total line flux
    #r_line_total = lzifu_hdulist['R_LINE'].data

    # Using my own line fits instead of lzifu
    # calling it r_line_total so I don't have to change any of the rest of the code
    r_line_total = np.load(taffy_extdir + 'halpha_profile_mylinefits.npy')

    # create wavelength array
    # I read these data from the corresponding headers
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # select a line to draw contours for
    # make sure that its rest wavelength is its wavelength in air!! the IFU data was taken on the ground
    # then make sure that you are looking for the line in the appropriate wav arr 
    # e.g. hbeta in blue and halpha in red 
    # find line index in wavelength array
    redshift = 0.0145  # average z 
    sys_vel = redshift * lightspeed  # avg systemic velocity is z*c = 4350 km/s
    print "Systemic Velocity (km/s):", sys_vel
    halpha_air_wav = 6562.80
    halpha_wav = halpha_air_wav*(1+redshift)
    halpha_idx = np.argmin(abs(red_wav_arr - halpha_wav))

    # convert wavelengths to heliocentric velocities
    helio_vel_arr = ((red_wav_arr - halpha_air_wav) / halpha_air_wav) * lightspeed

    # create figure and axis grid
    # this has to be done before looping over each vel range
    # also defining custom color list
    gs = gridspec.GridSpec(4,4)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure(figsize=(12,12))  # small for testing; make large for final saved figure # figsize=(12,12))

    # color palette from colorbrewer2.org
    cm = vcm.get_colorbrewer_cm('blues')

    # arbitrarily choosing -450 to +450 km/s as the velo range over which to show channel maps 
    # --> and the step I want for channel maps is 35 km/s to Nyquist sample the channels since the resolution 
    # at halpha is 70 km/s. Also I need it to be on a grid (that I chose to be 5x5).
    # On new figure (as requested by referee) velocity range is -540 to +510 km/s
    vel_step = 70.0
    low_vel_lim = sys_vel - 540.0
    high_vel_lim = sys_vel + 470.0

    vel_range_arr = np.arange(low_vel_lim, high_vel_lim+vel_step, vel_step)

    print "Central velocities wrt heliocentric for channels being plotted", vel_range_arr
    print "Central velocities wrt systemic for channels being plotted:", vel_range_arr - sys_vel
    print "Total channels to plot:", len(vel_range_arr)

    # For bold velocity label
    # modify font related rc Params (mostly to get bold text)
    # Cuz you can't use bold with TeX for some reason
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
    mpl.rcParams["text.usetex"] = False
    mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

    f = FontProperties()
    f.set_weight('bold')

    for j in range(len(vel_range_arr)):

        # first plot the sdss image
        row, col = divmod(j, 4)  # dividing by number of columns
        ax = fig.add_subplot(gs[row, col], projection=wcs_sdss)

        # overlay the contours on sdss image
        # get the correct sized image and normalization and plot
        norm = ImageNormalize(sdss_i[0].data, interval=ZScaleInterval(), stretch=LogStretch())
        orig_cmap = mpl.cm.Greys
        shifted_cmap = vcm.shiftedColorMap(orig_cmap, midpoint=0.6, name='shifted')
        im = ax.imshow(sdss_i[0].data, origin='lower', cmap=shifted_cmap, vmin=0, vmax=3, norm=norm)
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
        vel_mean_arr = np.mean(r_line_total[vel_range_idx], axis=0)

        vel_mean_arr = ma.array(vel_mean_arr, mask=all_mask)
        vel_mean_arr = ma.filled(vel_mean_arr, fill_value=0.0)
        vel_mean_arr = np.nan_to_num(vel_mean_arr)

        x = np.arange(vel_mean_arr.shape[1])
        y = np.arange(vel_mean_arr.shape[0])
        X, Y = np.meshgrid(x, y)

        levels = np.array([1,4,10,20,40,80,105,135])

        c = ax.contour(X, Y, vel_mean_arr, transform=ax.get_transform(wcs_lzifu), levels=levels, cmap=cm)
        #ax.clabel(c, inline=True, inline_spacing=0, fontsize=5, fmt='%1.1f', lw=3, ls='-')

        # remove all ticks, ticklabels, and spines
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.get_xaxis().set_ticklabels([])
        ax.get_yaxis().set_ticklabels([])

        ax.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=False, left=False, right=False, labelleft=False)

        # put in velocity range label
        vel_range_str = str(int(vel_range_low)) + " to " + str(int(vel_range_high))
        ax.text(0.07, 0.12, vel_range_str, verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=15)

        # Save as a fits file for velocity range 
        if int(vel_range_low) == 4156 and int(vel_range_high) == 4226:
            print "Saving mean velocity array for range 4156 to 4226 [km/s]."
            # get a header from one of the data cubes for the WCS
            hdr = lzifu_hdulist['R_LINE'].header

            # Set all values that are exactly zero back to NaN again
            # This is to be able to see it better in ds9
            set_to_nan_idx = np.where(vel_mean_arr == 0.0)
            vel_mean_arr[set_to_nan_idx] = np.nan

            # Create file and save
            hdu = fits.PrimaryHDU(data=vel_mean_arr, header=hdr)
            hdu.writeto(taffy_extdir + 'halpha_vel_range_4156_to_4226.fits', overwrite=True)

    #plt.show()
    fig.savefig(taffy_extdir + 'figures_stitched_cube/vel_channel_halpha_revised.pdf', dpi=250, bbox_inches='tight')

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)