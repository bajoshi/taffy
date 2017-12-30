from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.visualization import ManualInterval, ZScaleInterval, LogStretch, ImageNormalize
from astropy.convolution import Gaussian1DKernel, convolve

import sys
import os
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # does not have a trailing slash
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

sys.path.append(taffy_dir + 'codes/')
import vel_channel_map as vcm

# this following rgb_to_hex function came from
# https://stackoverflow.com/a/214657
def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def plotspectrum(axesblue, axesred1, axesred2, bluewav, bluespec, redwav, redspec, \
    b_ylow, b_yhigh, r1_ylow, r1_yhigh, r2_ylow, r2_yhigh, \
    b_xlow, b_xhigh, regionname):

    # try smoothing the spectra first so that features are easily visible
    # Create kernel
    gauss1d = Gaussian1DKernel(stddev=4)
    
    # Convolve data
    blue_smooth = convolve(bluespec, gauss1d)
    red_smooth = convolve(redspec, gauss1d)

    # plot spectra
    axesblue.plot(bluewav, blue_smooth, color=myblue)

    # plot the red spectra to show [OI] lines on one axes and 
    # H-alpha + [NII] and [SII] line on the other axes
    axesred1.plot(xgal_hii_red['wav'], red_smooth, color=myred)
    axesred2.plot(xgal_hii_red['wav'], red_smooth, color=myred)

    # set limits
    # limits have to be set separately for both broken red axes
    axesblue.set_ylim(b_ylow, b_yhigh + b_yhigh * 0.1)
    axesred1.set_ylim(r1_ylow, r1_yhigh)
    axesred2.set_ylim(r2_ylow, r2_yhigh)

    axesblue.set_xlim(b_xlow, b_xhigh)
    axesred1.set_xlim(6375,6475)
    axesred2.set_xlim(6600,6855)

    # commands to make the plot pretty
    # do not show the spines except at the bottom
    axesblue.spines["top"].set_visible(False)
    #axesblue.spines["bottom"].set_visible(False)
    axesblue.spines["right"].set_visible(False)
    axesblue.spines["left"].set_visible(False)

    axesred1.spines["top"].set_visible(False)
    #axesred1.spines["bottom"].set_visible(False)
    axesred1.spines["right"].set_visible(False)
    axesred1.spines["left"].set_visible(False)
    axesred2.spines["top"].set_visible(False)
    #axesred2.spines["bottom"].set_visible(False)
    axesred2.spines["right"].set_visible(False)
    axesred2.spines["left"].set_visible(False)

    # do not show the ticks on the left and bottom
    axesblue.get_yaxis().set_ticklabels([])
    axesblue.get_yaxis().set_ticks([])

    axesred1.get_yaxis().set_ticklabels([])
    axesred1.get_yaxis().set_ticks([])

    axesred2.get_yaxis().set_ticklabels([])
    axesred2.get_yaxis().set_ticks([])

    # better xaxis ticklabels
    axesblue.set_xticklabels(axesblue.get_xticks().tolist(), size=6, rotation=30)
    axesblue.tick_params('both', which='major', pad=1)

    axesred1.set_xticklabels(axesred1.get_xticks().tolist(), size=6, rotation=30)
    axesred1.tick_params('both', which='major', pad=1)

    axesred2.set_xticklabels(axesred2.get_xticks().tolist(), size=6, rotation=30)
    axesred2.tick_params('both', which='major', pad=1)

    # show region label
    box = TextArea(regionname, textprops=dict(color='k', size=9))
    anc_box = AnchoredOffsetbox(loc=2, child=box, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.03, 0.9),\
                                         bbox_transform=axesblue.transAxes, borderpad=0.0)
    axesblue.add_artist(anc_box)

    return axesblue, axesred1, axesred2

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in g band SDSS image
    #sdss_g = fits.open(home + '/Desktop/ipac/taffy_lzifu/taffy_xliners_figs_misc_data/taffy/SDSS/taffyg_sdss.fits')
    #wcs = WCS(sdss_g[0].header)

    # read in all integrated spectra
    # north nucleus
    nnuc_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/ngal_nuc_blue.dat', dtype=None, names=['wav','flux'])
    nnuc_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/ngal_nuc_red.dat', dtype=None, names=['wav','flux'])

    # south nucleus
    snuc_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/sgal_nuc_blue.dat', dtype=None, names=['wav','flux'])
    snuc_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/sgal_nuc_red.dat', dtype=None, names=['wav','flux'])

    # extragalactic HII region
    xgal_hii_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/xhii_blue.dat', dtype=None, names=['wav','flux'])
    xgal_hii_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/xhii_red.dat', dtype=None, names=['wav','flux'])

    # north galaxy west region 
    nw_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/nw_blue.dat', dtype=None, names=['wav','flux'])
    nw_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/nw_red.dat', dtype=None, names=['wav','flux'])

    # north galaxy east region 
    ne_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/ne_blue.dat', dtype=None, names=['wav','flux'])
    ne_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/ne_red.dat', dtype=None, names=['wav','flux'])

    # south galaxy west region 
    sw_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/sw_blue.dat', dtype=None, names=['wav','flux'])
    sw_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/sw_red.dat', dtype=None, names=['wav','flux'])

    # bridge west region 
    bw_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bw_blue.dat', dtype=None, names=['wav','flux'])
    bw_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bw_red.dat', dtype=None, names=['wav','flux'])

    # bridge east region
    be_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/be_blue.dat', dtype=None, names=['wav','flux'])
    be_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/be_red.dat', dtype=None, names=['wav','flux'])

    # Make figure
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # create figure and plot
    # the sdss image is plotted according the the function in the 
    # vel_channel_map code it is slightly modified here to work with gridspec.
    fig = plt.figure(figsize=(9,9))  # figsize=(width, height)

    # modify rc Params
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
    #mpl.rcParams["text.usetex"] = True  
    # I've kept this above commented out line here as a reminder. 
    # This line will stop any bold text from appearing anywhere on your figure.
    # Didn't realize this for a couple hours lol...
    mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"

    # define colors
    myblue = rgb_to_hex(0, 100, 180)
    myred = rgb_to_hex(214, 39, 40)  # tableau 20 red

    # create grid to plot
    gs = gridspec.GridSpec(24,40)
    gs.update(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=30, hspace=40)

    # define axes using above grid
    ax = fig.add_subplot(gs[5:20,0:20], projection=wcs_sdss)

    # the axes are ordered here going colockwise on the actual plot
    # create another axes for each red spectrum axis to show 
    # broken axes effect
    ax_nnuc_b = fig.add_subplot(gs[0:4,0:8])
    ax_nnuc_r1 = fig.add_subplot(gs[0:4,8:14])
    ax_nnuc_r2 = fig.add_subplot(gs[0:4,14:20])

    ax_xgal_b = fig.add_subplot(gs[0:4,20:28])
    ax_xgal_r1 = fig.add_subplot(gs[0:4,28:34])
    ax_xgal_r2 = fig.add_subplot(gs[0:4,34:40])

    ax_nw_b = fig.add_subplot(gs[4:8,20:28])
    ax_nw_r1 = fig.add_subplot(gs[4:8,28:34])
    ax_nw_r2 = fig.add_subplot(gs[4:8,34:40])

    ax_ne_b = fig.add_subplot(gs[8:12,20:28])
    ax_ne_r1 = fig.add_subplot(gs[8:12,28:34])
    ax_ne_r2 = fig.add_subplot(gs[8:12,34:40])

    ax_sw_b = fig.add_subplot(gs[12:16,20:28])
    ax_sw_r1 = fig.add_subplot(gs[12:16,28:34])
    ax_sw_r2 = fig.add_subplot(gs[12:16,34:40])

    ax_bw_b = fig.add_subplot(gs[16:20,20:28])
    ax_bw_r1 = fig.add_subplot(gs[16:20,28:34])
    ax_bw_r2 = fig.add_subplot(gs[16:20,34:40])

    ax_snuc_b = fig.add_subplot(gs[20:24,0:8])
    ax_snuc_r1 = fig.add_subplot(gs[20:24,8:14])
    ax_snuc_r2 = fig.add_subplot(gs[20:24,14:20])

    ax_be_b = fig.add_subplot(gs[20:24,20:28])
    ax_be_r1 = fig.add_subplot(gs[20:24,28:34])
    ax_be_r2 = fig.add_subplot(gs[20:24,34:40])

    # ------------ Plot galaxy ------------ #
    # SDSS g band image with proper stretch
    # set axis labels
    ax.set_xlabel('Right Ascension', fontsize=8)
    ax.set_ylabel('Declination', fontsize=8)

    # get the correct sized image and normalization and plot
    norm = ImageNormalize(sdss_i[0].data, stretch=LogStretch())
    orig_cmap = mpl.cm.Greys
    shifted_cmap = vcm.shiftedColorMap(orig_cmap, midpoint=0.5, name='shifted')
    #im = ax.imshow(sdss_i[0].data, origin='lower', cmap=shifted_cmap, vmin=0, vmax=3, norm=norm)
    im = ax.imshow(sdss_i[0].data, origin='lower', cmap=shifted_cmap, vmin=0.05, vmax=8, norm=norm)

    ax.set_autoscale_on(False)

    lon = ax.coords[0]
    lat = ax.coords[1]

    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')

    ax.coords.frame.set_color('None')

    # ------------------------ PLOT THE REGIONS ------------------------ #
    # the format for the coordinates here is 
    # Rectangle((lower_left_x, lower_left_y), width, height, ....)
    # the -1 is to go from ds9 x,y to numpy row,col
    # i've translated the center coords that i got from ds9
    # to the coords of the lower left corner
    north_nuc = SphericalCircle((0.4245000045 * u.deg, 23.49619722 * u.deg), 0.0017047 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    south_nuc = SphericalCircle((0.4099666595 * u.deg, 23.48356417 * u.deg), 0.0017153 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(north_nuc)
    ax.add_patch(south_nuc)

    xgal_hii = SphericalCircle((0.4208500067 * u.deg, 23.49327528 * u.deg), 0.0020958 * u.degree, \
        edgecolor='green', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(xgal_hii)

    nw_patch = SphericalCircle((0.4192749977 * u.deg, 23.50056333 * u.deg), 0.0018061 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ne_patch = SphericalCircle((0.4283333302 * u.deg, 23.49243528 * u.deg), 0.0018061 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    sw_patch = SphericalCircle((0.4066708406 * u.deg, 23.49171278 * u.deg), 0.0018061 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(nw_patch)
    ax.add_patch(ne_patch)
    ax.add_patch(sw_patch)

    bw_patch = SphericalCircle((0.4149416606 * u.deg, 23.49478361 * u.deg), 0.0018061 * u.degree, \
        edgecolor='deeppink', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    be_patch = SphericalCircle((0.4228208383 * u.deg, 23.48801056 * u.deg), 0.0018061 * u.degree, \
        edgecolor='deeppink', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(bw_patch)
    ax.add_patch(be_patch)

    # add text to figure to indicate region name
    f = FontProperties()
    f.set_weight('bold')

    ax.text(0.36, 0.78, 'North galaxy nucleus', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.68, 0.41, 'South galaxy' + '\n' + 'nucleus', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.41, 0.61, 'Extragalactic' + '\n' + 'HII region', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.48, 0.87, 'North galaxy' + '\n' + 'west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.05, 0.59, 'North galaxy' + '\n' + 'east', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.76, 0.64, 'South galaxy' + '\n' + 'west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.57, 0.72, 'Bridge west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    ax.text(0.26, 0.49, 'Bridge east', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=11)

    # ------------------------ PLOT THE SPECTRA ------------------------ #
    # ------------ Plot north nuclear region spectrum ------------ #
    ax_nnuc_b, ax_nnuc_r1, ax_nnuc_r2 = \
    plotspectrum(ax_nnuc_b, ax_nnuc_r1, ax_nnuc_r2, nnuc_blue['wav'], nnuc_blue['flux'], nnuc_red['wav'], nnuc_red['flux'],\
        25,53, 58,78, 48,150, 4700,5320, 'North galaxy nucleus')

    # ------------ Plot south nuclear region spectrum ------------ #
    ax_snuc_b, ax_snuc_r1, ax_snuc_r2 = \
    plotspectrum(ax_snuc_b, ax_snuc_r1, ax_snuc_r2, snuc_blue['wav'], snuc_blue['flux'], snuc_red['wav'], snuc_red['flux'],\
        50,90, 70,105, 88,162, 4700,5320, 'South galaxy nucleus')

    # ------------ Plot X-Galactic HII region spectrum ------------ #
    ax_xgal_b, ax_xgal_r1, ax_xgal_r2 = \
    plotspectrum(ax_xgal_b, ax_xgal_r1, ax_xgal_r2, xgal_hii_blue['wav'], xgal_hii_blue['flux'], xgal_hii_red['wav'], xgal_hii_red['flux'],\
        0,30, 0,50, -10,100, 4850,5150, 'Extragalactic HII region')

    # ------------ Plot north galaxy west region spectrum ------------ #
    ax_nw_b, ax_nw_r1, ax_nw_r2 = \
    plotspectrum(ax_nw_b, ax_nw_r1, ax_nw_r2, nw_blue['wav'], nw_blue['flux'], nw_red['wav'], nw_red['flux'],\
        10,50, 15,40, 10,105, 4700,5320, 'North galaxy west')

    # ------------ Plot north galaxy east region spectrum ------------ #
    ax_ne_b, ax_ne_r1, ax_ne_r2 = \
    plotspectrum(ax_ne_b, ax_ne_r1, ax_ne_r2, ne_blue['wav'], ne_blue['flux'], ne_red['wav'], ne_red['flux'],\
        15,50, 15,45, 10,55, 4700,5320, 'North galaxy east')

    # ------------ Plot south galaxy west region spectrum ------------ #
    ax_sw_b, ax_sw_r1, ax_sw_r2 = \
    plotspectrum(ax_sw_b, ax_sw_r1, ax_sw_r2, sw_blue['wav'], sw_blue['flux'], sw_red['wav'], sw_red['flux'],\
        10,65, 15,30, 10,150, 4700,5320, 'South galaxy west')

    # ------------ Plot bridge west region spectrum ------------ #
    ax_bw_b, ax_bw_r1, ax_bw_r2 = \
    plotspectrum(ax_bw_b, ax_bw_r1, ax_bw_r2, bw_blue['wav'], bw_blue['flux'], bw_red['wav'], bw_red['flux'],\
        -5,5, -5,10, -5,20, 4850,5150, 'Bridge west')

    # ------------ Plot bridge east region spectrum ------------ #
    ax_be_b, ax_be_r1, ax_be_r2 = \
    plotspectrum(ax_be_b, ax_be_r1, ax_be_r2, be_blue['wav'], be_blue['flux'], be_red['wav'], be_red['flux'],\
        -5,5, -5,10, -8,10, 4850,5150, 'Bridge east')

    # ------------ Save figure ------------ #
    fig.savefig(taffy_extdir + 'figures_stitched_cube/taffy_intg_spec.eps', dpi=600, bbox_inches='tight')
    #plt.show()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)