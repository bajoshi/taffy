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
from matplotlib.patches import Rectangle, Polygon
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

    # try smoothing the spectra first so 
    # that features are easily visible
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
    axesblue.spines["right"].set_visible(False)
    axesblue.spines["left"].set_visible(False)

    axesred1.spines["top"].set_visible(False)
    axesred1.spines["right"].set_visible(False)
    axesred1.spines["left"].set_visible(False)

    axesred2.spines["top"].set_visible(False)
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
    axesblue.text(0.03, 0.9, regionname, verticalalignment='top', horizontalalignment='left', \
        transform=axesblue.transAxes, color='k', size=9)

    # plot line labels according to region name
    if 'HII' in regionname:
        axesblue.text(0.3, 0.63, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.34, 0.45, r'$\mathrm{[OIII]}$' + '\n' +  r'$\lambda 4959$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.81, 0.6, r'$\mathrm{[OIII]}$' + '\n' +  r'$\lambda 5007$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.15, 0.65, r'$\mathrm{[OI]\lambda 6300}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.27, 0.75, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(-0.05, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.35, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.6, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'Bridge east' in regionname:
        axesred1.text(0.15, 0.75, r'$\mathrm{[OI]\lambda}$' + '\n' + r'$6300$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.25, 0.84, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(-0.05, 0.7, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.33, 0.7, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.84, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'Bridge west' in regionname:
        axesblue.text(0.3, 0.63, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.1, 0.3, r'$\mathrm{[OI]\lambda}$' + '\n' + r'$6300$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.25, 0.75, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(-0.05, 0.55, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.35, 0.55, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.78, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'South galaxy east' in regionname:
        axesblue.text(0.4, 0.57, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.5, 0.43, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.25, 0.6, r'$\mathrm{[OI]\lambda}$' + '\n' + r'$6300$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.27, 0.85, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0, 0.55, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.35, 0.55, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.7, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'South galaxy west' in regionname:
        axesblue.text(0.24, 0.57, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.47, 0.52, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred2.text(0.23, 0.63, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(-0.05, 0.38, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.33, 0.38, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.5, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'North galaxy east' in regionname:
        axesblue.text(0.35, 0.15, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.5, 0.3, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred2.text(0.2, 0.72, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.32, 0.65, r'$\mathrm{[NII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.8, 0.65, r'$\mathrm{[SII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'North galaxy west' in regionname:
        axesblue.text(0.3, 0.15, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.5, 0.2, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.15, 0.7, r'$\mathrm{[OI]\lambda 6300}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred1.text(0.7, 0.6, r'$\mathrm{[OI]\lambda 6364}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.27, 0.85, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.38, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.58, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'nucleus' in regionname and 'South' in regionname:
        axesblue.text(0.3, 0.2, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.5, 0.2, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.03, 0.15, r'$\mathrm{[OI]\lambda 6300}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.2, 0.9, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0, 0.5, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.35, 0.5, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.8, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    elif 'nucleus' in regionname and 'North' in regionname:
        axesblue.text(0.3, 0.2, r'$\mathrm{H\beta}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesblue.text(0.5, 0.3, r'$\mathrm{[OIII]}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesblue.transAxes, color='k', size=7)
        axesred1.text(0.15, 0.65, r'$\mathrm{[OI]\lambda 6300}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred1.text(0.7, 0.5, r'$\mathrm{[OI]\lambda 6364}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred1.transAxes, color='k', size=7)
        axesred2.text(0.27, 0.85, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' + r'$6548$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.38, 0.4, r'$\mathrm{[NII]}\lambda$' + '\n' +  r'$6583$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)
        axesred2.text(0.75, 0.53, r'$\mathrm{[SII]\lambda\lambda}$' + '\n' + r'$6717,$' + '\n' + r'$6731$', \
            verticalalignment='top', horizontalalignment='left', \
            transform=axesred2.transAxes, color='k', size=7)

    return axesblue, axesred1, axesred2

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

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

    # south galaxy west region 
    se_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/se_blue.dat', dtype=None, names=['wav','flux'])
    se_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/se_red.dat', dtype=None, names=['wav','flux'])

    # bridge west region 
    bw_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bw_new_blue.dat', dtype=None, names=['wav','flux'])
    bw_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bw_new_red.dat', dtype=None, names=['wav','flux'])

    # bridge east region
    be_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/be_new_blue.dat', dtype=None, names=['wav','flux'])
    be_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/be_new_red.dat', dtype=None, names=['wav','flux'])

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
    # This line will stop any bold text from appearing anywhere 
    # on your figure. Didn't realize this for a couple hours lol...
    # Also, this single line increases the run time for this code
    # from ~6 seconds to ~20 seconds.
    mpl.rcParams["text.usetex"] = False
    # Explicitly setting this to False now
    # because my rc file has this set to True
    mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"

    # define colors
    myblue = rgb_to_hex(0, 100, 180)
    myred = rgb_to_hex(214, 39, 40)  # tableau 20 red

    # create grid to plot
    gs = gridspec.GridSpec(60,60)
    gs.update(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=10, hspace=10)

    # define axes using above grid
    ax = fig.add_subplot(gs[10:40,0:30], projection=wcs_sdss)

    # the axes are ordered here to make sure the
    # x axis tick labels show up (overlapped) on 
    # the axis space below.
    ax_nnuc_b = fig.add_subplot(gs[0:10,0:10])
    ax_nnuc_r1 = fig.add_subplot(gs[0:10,10:20])
    ax_nnuc_r2 = fig.add_subplot(gs[0:10,20:30])

    ax_bw_b = fig.add_subplot(gs[50:60,30:40])
    ax_bw_r1 = fig.add_subplot(gs[50:60,40:50])
    ax_bw_r2 = fig.add_subplot(gs[50:60,50:60])

    ax_se_b = fig.add_subplot(gs[40:50,30:40])
    ax_se_r1 = fig.add_subplot(gs[40:50,40:50])
    ax_se_r2 = fig.add_subplot(gs[40:50,50:60])

    ax_sw_b = fig.add_subplot(gs[30:40,30:40])
    ax_sw_r1 = fig.add_subplot(gs[30:40,40:50])
    ax_sw_r2 = fig.add_subplot(gs[30:40,50:60])

    ax_ne_b = fig.add_subplot(gs[20:30,30:40])
    ax_ne_r1 = fig.add_subplot(gs[20:30,40:50])
    ax_ne_r2 = fig.add_subplot(gs[20:30,50:60])

    ax_nw_b = fig.add_subplot(gs[10:20,30:40])
    ax_nw_r1 = fig.add_subplot(gs[10:20,40:50])
    ax_nw_r2 = fig.add_subplot(gs[10:20,50:60])

    ax_snuc_b = fig.add_subplot(gs[0:10,30:40])
    ax_snuc_r1 = fig.add_subplot(gs[0:10,40:50])
    ax_snuc_r2 = fig.add_subplot(gs[0:10,50:60])

    ax_be_b = fig.add_subplot(gs[50:60,0:10])
    ax_be_r1 = fig.add_subplot(gs[50:60,10:20])
    ax_be_r2 = fig.add_subplot(gs[50:60,20:30])

    ax_xgal_b = fig.add_subplot(gs[40:50,0:10])
    ax_xgal_r1 = fig.add_subplot(gs[40:50,10:20])
    ax_xgal_r2 = fig.add_subplot(gs[40:50,20:30])

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
    se_patch = SphericalCircle((0.4134277924 * u.deg, 23.47796085 * u.deg), 0.0018061 * u.degree, \
        edgecolor='dodgerblue', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(nw_patch)
    ax.add_patch(ne_patch)
    ax.add_patch(sw_patch)
    ax.add_patch(se_patch)

    bw_patch = SphericalCircle((0.4162353178 * u.deg, 23.49744626 * u.deg), 0.0018061 * u.degree, \
        edgecolor='deeppink', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    be_patch = SphericalCircle((0.4177249349 * u.deg, 23.49056863 * u.deg), 0.0018557 * u.degree, \
        edgecolor='deeppink', facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(bw_patch)
    ax.add_patch(be_patch)

    # Add polygon patches for galaxies and bridge
    north_poly_points = [[0.422200173,23.50396718],[0.4152304334,23.50387652],[0.4155502808,23.50226454],\
    [0.414591942,23.49908922],[0.4180543574,23.49850344],[0.4199720955,23.49723355],[0.4259616693,23.49127482],\
    [0.4286967795,23.48904169],[0.4291162163,23.49004197],[0.4313178761,23.49288664],[0.4291871345,23.49689221]]
    south_poly_points = [[0.4071288157,23.49766402],[0.4027494899,23.49619723],[0.4031346441,23.48922165],\
    [0.4060123993,23.48189481],[0.4078664731,23.47872971],[0.410943916,23.47614935],[0.4158554494,23.47630534],\
    [0.4202012426,23.47791063],[0.4217348365,23.48048999],[0.4213999114,23.48310946],[0.417944079,23.48275696],\
    [0.4133603491,23.48707351],[0.4114718625,23.48932533],[0.4090425678,23.49400985]]
    bridge_poly_points = [[0.4071288316,23.49758635],[0.4090459205,23.49391271],[0.4114969024,23.48922357],\
    [0.4133614547,23.48697676],[0.4179534303,23.48268663],[0.421404655,23.48303866],[0.4251755324,23.48233549],\
    [0.4271566586,23.48509069],[0.4286953339,23.48900495],[0.4258300899,23.49137329],[0.4199687248,23.49723481],\
    [0.4180519391,23.49850339],[0.412549028,23.49960668]]

    pn = Polygon(np.array(north_poly_points), edgecolor='k', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(pn)
    ps = Polygon(np.array(south_poly_points), edgecolor='k', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(ps)
    pb = Polygon(np.array(bridge_poly_points), edgecolor='k', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=1.5)
    ax.add_patch(pb)

    # add text to figure to indicate region name
    f = FontProperties()
    f.set_weight('bold')

    ax.text(0.1, 0.82, 'North galaxy' + '\n' + ' nucleus', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.68, 0.41, 'South galaxy' + '\n' + 'nucleus', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.45, 0.71, 'Extragalactic' + '\n' + 'HII region', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.48, 0.89, 'North galaxy' + '\n' + 'west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.07, 0.64, 'North' + '\n' + 'galaxy' + '\n' + 'east', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.76, 0.64, 'South galaxy' + '\n' + 'west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.61, 0.29, 'South galaxy' + '\n' + 'east', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.54, 0.79, 'Bridge' + '\n' + 'west', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    ax.text(0.34, 0.55, 'Bridge east', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=9)

    # Add text for galaxies and bridge
    ax.text(0.16, 0.94, 'NORTH' + '\n' + 'GALAXY', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=10)
    ax.text(0.27, 0.29, 'SOUTH' + '\n' + 'GALAXY', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=10)
    ax.text(0.1, 0.45, 'BRIDGE', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=10)

    # ------------------------ PLOT THE SPECTRA ------------------------ #
    # these spectra are plotted in the same order as the axes above
    # the axes oreder DOES matter for the x axis ticks labels. 
    # see above comment on axes ordering.
    # the plotting order here doesn't really matter. 
    # I just did it to be consistent with the axes order.
    # ------------ Plot north nuclear region spectrum ------------ #
    ax_nnuc_b, ax_nnuc_r1, ax_nnuc_r2 = \
    plotspectrum(ax_nnuc_b, ax_nnuc_r1, ax_nnuc_r2, nnuc_blue['wav'], nnuc_blue['flux'], nnuc_red['wav'], nnuc_red['flux'],\
        25,53, 58,78, 48,150, 4700,5320, 'North galaxy' + '\n' + 'nucleus')

    # ------------ Plot bridge west region spectrum ------------ #
    ax_bw_b, ax_bw_r1, ax_bw_r2 = \
    plotspectrum(ax_bw_b, ax_bw_r1, ax_bw_r2, bw_blue['wav'], bw_blue['flux'], bw_red['wav'], bw_red['flux'],\
        -5,10, -5,10, -5,20, 4850,5150, 'Bridge west')

    # ------------ Plot south galaxy east region spectrum ------------ #
    ax_se_b, ax_se_r1, ax_se_r2 = \
    plotspectrum(ax_se_b, ax_se_r1, ax_se_r2, se_blue['wav'], se_blue['flux'], se_red['wav'], se_red['flux'],\
        10,60, 28,43, 10,150, 4700,5320, 'South galaxy east')

    # ------------ Plot south galaxy west region spectrum ------------ #
    ax_sw_b, ax_sw_r1, ax_sw_r2 = \
    plotspectrum(ax_sw_b, ax_sw_r1, ax_sw_r2, sw_blue['wav'], sw_blue['flux'], sw_red['wav'], sw_red['flux'],\
        10,65, 20,27, 10,140, 4700,5320, 'South galaxy west')

    # ------------ Plot north galaxy east region spectrum ------------ #
    ax_ne_b, ax_ne_r1, ax_ne_r2 = \
    plotspectrum(ax_ne_b, ax_ne_r1, ax_ne_r2, ne_blue['wav'], ne_blue['flux'], ne_red['wav'], ne_red['flux'],\
        15,50, 25,38, 10,55, 4700,5320, 'North galaxy east')

    # ------------ Plot north galaxy west region spectrum ------------ #
    ax_nw_b, ax_nw_r1, ax_nw_r2 = \
    plotspectrum(ax_nw_b, ax_nw_r1, ax_nw_r2, nw_blue['wav'], nw_blue['flux'], nw_red['wav'], nw_red['flux'],\
        10,50, 15,40, 10,105, 4700,5320, 'North galaxy west')

    # ------------ Plot south nuclear region spectrum ------------ #
    ax_snuc_b, ax_snuc_r1, ax_snuc_r2 = \
    plotspectrum(ax_snuc_b, ax_snuc_r1, ax_snuc_r2, snuc_blue['wav'], snuc_blue['flux'], snuc_red['wav'], snuc_red['flux'],\
        50,90, 85,105, 88,162, 4700,5320, 'South galaxy' + '\n' + 'nucleus')

    # ------------ Plot bridge east region spectrum ------------ #
    ax_be_b, ax_be_r1, ax_be_r2 = \
    plotspectrum(ax_be_b, ax_be_r1, ax_be_r2, be_blue['wav'], be_blue['flux'], be_red['wav'], be_red['flux'],\
        -5,8, -5,10, -8,14, 4850,5150, 'Bridge east')

    # ------------ Plot X-Galactic HII region spectrum ------------ #
    ax_xgal_b, ax_xgal_r1, ax_xgal_r2 = \
    plotspectrum(ax_xgal_b, ax_xgal_r1, ax_xgal_r2, xgal_hii_blue['wav'], xgal_hii_blue['flux'], xgal_hii_red['wav'], xgal_hii_red['flux'],\
        0,30, 0,35, -10,100, 4850,5150, 'Extragalactic' + '\n' + 'HII region')

    # ------------ Save figure ------------ #
    fig.savefig(taffy_extdir + 'figures_stitched_cube/taffy_intg_spec.eps', dpi=300, bbox_inches='tight')
    #plt.show()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)