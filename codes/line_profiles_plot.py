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

def plot_line_profiles(axesblue, axesred, bluewav, bluespec, redwav, redspec, \
    b_ylow, b_yhigh, r_ylow, r_yhigh):

    # plot spectra
    axesblue.plot(bluewav, bluespec, color=myblue)
    axesred.plot(redwav, redspec, color=myred)

    # Set limits
    axesblue.set_xlim(5070, 5090)
    axesred.set_xlim(6645, 6675)

    axesblue.set_ylim(b_ylow, b_yhigh)
    axesred.set_ylim(r_ylow, r_yhigh)

    # commands to make the plot pretty
    # do not show the spines except at the bottom
    axesblue.spines["top"].set_visible(False)
    axesblue.spines["right"].set_visible(False)
    axesblue.spines["left"].set_visible(False)

    axesred.spines["top"].set_visible(False)
    axesred.spines["right"].set_visible(False)
    axesred.spines["left"].set_visible(False)

    # do not show the ticks
    axesblue.get_yaxis().set_ticklabels([])
    axesblue.get_yaxis().set_ticks([])

    axesred.get_yaxis().set_ticklabels([])
    axesred.get_yaxis().set_ticks([])

    # for now set xticklabels to empty # it will be edited by hand in __main__
    axesblue.get_xaxis().set_ticklabels([])
    axesred.get_xaxis().set_ticklabels([])

    # Draw a vertical line at the systemic recessional velocity
    axesblue.axvline(x=5079.6, ls='--', color='grey')
    axesred.axvline(x=6658.2, ls='--', color='grey')

    return axesblue, axesred

if __name__ == '__main__':

    # --------------------------------------------------- Read in all integrated spectra --------------------------------------------------- #
    # These are averaged spectra (NOT summed) over the regions shown.
    # I don't think it matters too much since we're only trying to see
    # the line profiles with this plot.
    north_majax_1_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_1_blue.dat', dtype=None, names=['wav','flux'])
    north_majax_2_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_2_blue.dat', dtype=None, names=['wav','flux'])
    north_majax_3_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_3_blue.dat', dtype=None, names=['wav','flux'])
    north_majax_4_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_4_blue.dat', dtype=None, names=['wav','flux'])
    north_majax_5_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_5_blue.dat', dtype=None, names=['wav','flux'])

    north_majax_1_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_1_red.dat', dtype=None, names=['wav','flux'])
    north_majax_2_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_2_red.dat', dtype=None, names=['wav','flux'])
    north_majax_3_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_3_red.dat', dtype=None, names=['wav','flux'])
    north_majax_4_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_4_red.dat', dtype=None, names=['wav','flux'])
    north_majax_5_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/north_majax_5_red.dat', dtype=None, names=['wav','flux'])

    # ------------------------------- #
    south_majax_1_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_1_blue.dat', dtype=None, names=['wav','flux'])
    south_majax_2_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_2_blue.dat', dtype=None, names=['wav','flux'])
    south_majax_3_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_3_blue.dat', dtype=None, names=['wav','flux'])
    south_majax_4_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_4_blue.dat', dtype=None, names=['wav','flux'])
    south_majax_5_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_5_blue.dat', dtype=None, names=['wav','flux'])

    south_majax_1_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_1_red.dat', dtype=None, names=['wav','flux'])
    south_majax_2_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_2_red.dat', dtype=None, names=['wav','flux'])
    south_majax_3_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_3_red.dat', dtype=None, names=['wav','flux'])
    south_majax_4_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_4_red.dat', dtype=None, names=['wav','flux'])
    south_majax_5_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/south_majax_5_red.dat', dtype=None, names=['wav','flux'])

    # ------------------------------- #
    bridge_1_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_1_blue.dat', dtype=None, names=['wav','flux'])
    bridge_2_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_2_blue.dat', dtype=None, names=['wav','flux'])
    bridge_3_blue = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_3_blue.dat', dtype=None, names=['wav','flux'])

    bridge_1_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_1_red.dat', dtype=None, names=['wav','flux'])
    bridge_2_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_2_red.dat', dtype=None, names=['wav','flux'])
    bridge_3_red = np.genfromtxt(taffy_extdir + 'rawspectra_for_paperplot/bridge_3_red.dat', dtype=None, names=['wav','flux'])

    # --------------------------------------------------- Plotting --------------------------------------------------- #
    # ----------- Plot SDSS image first ----------- #
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # create figure and plot
    # the sdss image is plotted according the the function in the 
    # vel_channel_map code it is slightly modified here to work with gridspec.
    fig = plt.figure(figsize=(16, 12))  # figsize=(width, height)

    # modify rc Params
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
    mpl.rcParams["text.usetex"] = False
    mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
    mpl.rcParams["xtick.top"] = False

    # define colors
    myblue = rgb_to_hex(0, 100, 180)
    myred = rgb_to_hex(214, 39, 40)  # tableau 20 red

    # create grid to plot
    gs = gridspec.GridSpec(60,84)
    gs.update(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=10, hspace=10)

    # define axes using above grid
    ax = fig.add_subplot(gs[:36,:36], projection=wcs_sdss)

    # ----------- Assign gridspec locations ----------- #
    # Bridge
    ax_br_1_blue = fig.add_subplot(gs[36:48,:12])
    ax_br_1_red  = fig.add_subplot(gs[48:,:12])

    ax_br_2_blue = fig.add_subplot(gs[36:48,12:24])
    ax_br_2_red  = fig.add_subplot(gs[48:,12:24])

    ax_br_3_blue = fig.add_subplot(gs[36:48,24:36])
    ax_br_3_red  = fig.add_subplot(gs[48:,24:36])

    # North galaxy
    ax_n_1_blue = fig.add_subplot(gs[:12,36:48])
    ax_n_2_blue = fig.add_subplot(gs[12:24,36:48])
    ax_n_3_blue = fig.add_subplot(gs[24:36,36:48])
    ax_n_4_blue = fig.add_subplot(gs[36:48,36:48])
    ax_n_5_blue = fig.add_subplot(gs[48:,36:48])

    ax_n_1_red  = fig.add_subplot(gs[:12,48:60])
    ax_n_2_red  = fig.add_subplot(gs[12:24,48:60])
    ax_n_3_red  = fig.add_subplot(gs[24:36,48:60])
    ax_n_4_red  = fig.add_subplot(gs[36:48,48:60])
    ax_n_5_red  = fig.add_subplot(gs[48:,48:60])

    # South galaxy
    ax_s_1_blue = fig.add_subplot(gs[:12,60:72])
    ax_s_2_blue = fig.add_subplot(gs[12:24,60:72])
    ax_s_3_blue = fig.add_subplot(gs[24:36,60:72])
    ax_s_4_blue = fig.add_subplot(gs[36:48,60:72])
    ax_s_5_blue = fig.add_subplot(gs[48:,60:72])

    ax_s_1_red  = fig.add_subplot(gs[:12,72:])
    ax_s_2_red  = fig.add_subplot(gs[12:24,72:])
    ax_s_3_red  = fig.add_subplot(gs[24:36,72:])
    ax_s_4_red  = fig.add_subplot(gs[36:48,72:])
    ax_s_5_red  = fig.add_subplot(gs[48:,72:])

    # ------------ Plot galaxy ------------ #
    # SDSS g band image with proper stretch
    # set axis labels
    ax.set_xlabel('Right Ascension', fontsize=8)
    ax.set_ylabel('Declination', fontsize=8)

    # get the correct sized image and normalization and plot
    norm = ImageNormalize(sdss_i[0].data, stretch=LogStretch())
    orig_cmap = mpl.cm.Greys
    shifted_cmap = vcm.shiftedColorMap(orig_cmap, midpoint=0.5, name='shifted')
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

    # Add text to figure to indicate region name
    f = FontProperties()
    f.set_weight('bold')

    ax.text(0.18, 0.7, 'N1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.23, 0.75, 'N2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.27, 0.8, 'N3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.32, 0.85, 'N4', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.38, 0.9, 'N5', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)

    ax.text(0.38, 0.63, 'B1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.52, 0.81, 'B2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.47, 0.52, 'B3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)

    ax.text(0.54, 0.3, 'S1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.65, 0.4, 'S2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.71, 0.49, 'S3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.75, 0.56, 'S4', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)
    ax.text(0.75, 0.64, 'S5', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=12)

    ax.text(0.05, 0.1, r'$\mathrm{[OIII]\lambda 5007}$', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color=myblue, size=16)
    ax.text(0.05, 0.06, r'$\mathrm{H\alpha}$', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color=myred, size=16)

    # ------------------------ PLOT THE REGIONS ------------------------ #
    reg_file = open(taffy_extdir + 'rawspectra_for_paperplot/regions_for_line_profiles.reg')
    circ_size = 3.776 / 3600 

    for ln_count, ln in enumerate(reg_file.readlines()):
        if ln_count < 3:
            continue
        dec = float(ln.split(',')[1])
        ra = float(ln.split(',')[0].split('circle(')[1])

        reg_color = 'forestgreen'
        if 'color=blue' in ln:
            reg_color = 'midnightblue'

        region = SphericalCircle((ra * u.deg, dec * u.deg), circ_size * u.degree, \
            edgecolor=reg_color, facecolor='none', transform=ax.get_transform('fk5'), lw=1.5)
        ax.add_patch(region)

    # ------------------------ PLOT THE SPECTRA ------------------------ #
    ax_n_1_blue, ax_n_1_red = plot_line_profiles(ax_n_1_blue, ax_n_1_red, north_majax_1_blue['wav'], north_majax_1_blue['flux'], north_majax_1_red['wav'], north_majax_1_red['flux'], 55,90, 40,175)
    ax_n_2_blue, ax_n_2_red = plot_line_profiles(ax_n_2_blue, ax_n_2_red, north_majax_2_blue['wav'], north_majax_2_blue['flux'], north_majax_2_red['wav'], north_majax_2_red['flux'], 40,80, 60,320)
    ax_n_3_blue, ax_n_3_red = plot_line_profiles(ax_n_3_blue, ax_n_3_red, north_majax_3_blue['wav'], north_majax_3_blue['flux'], north_majax_3_red['wav'], north_majax_3_red['flux'], 20,48, 20,235)
    ax_n_4_blue, ax_n_4_red = plot_line_profiles(ax_n_4_blue, ax_n_4_red, north_majax_4_blue['wav'], north_majax_4_blue['flux'], north_majax_4_red['wav'], north_majax_4_red['flux'], 15,50, 15,95)
    ax_n_5_blue, ax_n_5_red = plot_line_profiles(ax_n_5_blue, ax_n_5_red, north_majax_5_blue['wav'], north_majax_5_blue['flux'], north_majax_5_red['wav'], north_majax_5_red['flux'], 20,75, 0,175)

    ax_s_1_blue, ax_s_1_red = plot_line_profiles(ax_s_1_blue, ax_s_1_red, south_majax_1_blue['wav'], south_majax_1_blue['flux'], south_majax_1_red['wav'], south_majax_1_red['flux'], 24,58, 30,350)
    ax_s_2_blue, ax_s_2_red = plot_line_profiles(ax_s_2_blue, ax_s_2_red, south_majax_2_blue['wav'], south_majax_2_blue['flux'], south_majax_2_red['wav'], south_majax_2_red['flux'], 100,160, 140,260)
    ax_s_3_blue, ax_s_3_red = plot_line_profiles(ax_s_3_blue, ax_s_3_red, south_majax_3_blue['wav'], south_majax_3_blue['flux'], south_majax_3_red['wav'], south_majax_3_red['flux'], 25,63, 40,275)
    ax_s_4_blue, ax_s_4_red = plot_line_profiles(ax_s_4_blue, ax_s_4_red, south_majax_4_blue['wav'], south_majax_4_blue['flux'], south_majax_4_red['wav'], south_majax_4_red['flux'], 40,120, 40,340)
    ax_s_5_blue, ax_s_5_red = plot_line_profiles(ax_s_5_blue, ax_s_5_red, south_majax_5_blue['wav'], south_majax_5_blue['flux'], south_majax_5_red['wav'], south_majax_5_red['flux'], 25,75, 20,130)

    ax_br_1_blue, ax_br_1_red = plot_line_profiles(ax_br_1_blue, ax_br_1_red, bridge_1_blue['wav'], bridge_1_blue['flux'], bridge_1_red['wav'], bridge_1_red['flux'], -7,75, -15,200)
    ax_br_2_blue, ax_br_2_red = plot_line_profiles(ax_br_2_blue, ax_br_2_red, bridge_2_blue['wav'], bridge_2_blue['flux'], bridge_2_red['wav'], bridge_2_red['flux'], -10,10, -10,45)
    ax_br_3_blue, ax_br_3_red = plot_line_profiles(ax_br_3_blue, ax_br_3_red, bridge_3_blue['wav'], bridge_3_blue['flux'], bridge_3_red['wav'], bridge_3_red['flux'], -10,25, -10,50)

    # Wavelength axis labels
    ax_n_5_blue.set_xticklabels(ax_n_5_blue.get_xticks().tolist(), size=10, rotation=20)
    ax_n_5_blue.tick_params('both', which='major', pad=-1)
    ax_n_5_red.set_xticklabels(ax_n_5_red.get_xticks().tolist(), size=10, rotation=20)
    ax_n_5_red.tick_params('both', which='major', pad=1)

    ax_s_5_blue.set_xticklabels(ax_s_5_blue.get_xticks().tolist(), size=10, rotation=20)
    ax_s_5_blue.tick_params('both', which='major', pad=1)
    ax_s_5_red.set_xticklabels(ax_s_5_red.get_xticks().tolist(), size=10, rotation=20)
    ax_s_5_red.tick_params('both', which='major', pad=-1)

    #ax_br_1_blue.set_xticklabels(ax_br_1_blue.get_xticks().tolist(), size=10, rotation=20)
    ax_br_1_red.set_xticklabels(ax_br_1_red.get_xticks().tolist(), size=10, rotation=20)

    #ax_br_2_blue.set_xticklabels(ax_br_2_blue.get_xticks().tolist(), size=10, rotation=20)
    ax_br_2_red.set_xticklabels(ax_br_2_red.get_xticks().tolist(), size=10, rotation=20)

    #ax_br_3_blue.set_xticklabels(ax_br_3_blue.get_xticks().tolist(), size=10, rotation=20)
    ax_br_3_red.set_xticklabels(ax_br_3_red.get_xticks().tolist(), size=10, rotation=20)

    # ------------------ Text for regionname ------------------- #
    ax_n_1_red.text(0.0, 0.98, 'N1', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_1_red.transAxes, color='k', fontproperties=f, size=12)
    ax_n_2_red.text(0.0, 0.98, 'N2', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_2_red.transAxes, color='k', fontproperties=f, size=12)
    ax_n_3_red.text(0.0, 0.98, 'N3', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_3_red.transAxes, color='k', fontproperties=f, size=12)
    ax_n_4_red.text(0.0, 0.98, 'N4', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_4_red.transAxes, color='k', fontproperties=f, size=12)
    ax_n_5_red.text(0.0, 0.98, 'N5', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_5_red.transAxes, color='k', fontproperties=f, size=12)

    ax_s_1_red.text(0.0, 0.98, 'S1', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_1_red.transAxes, color='k', fontproperties=f, size=12)
    ax_s_2_red.text(0.0, 0.98, 'S2', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_2_red.transAxes, color='k', fontproperties=f, size=12)
    ax_s_3_red.text(0.0, 0.98, 'S3', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_3_red.transAxes, color='k', fontproperties=f, size=12)
    ax_s_4_red.text(0.0, 0.98, 'S4', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_4_red.transAxes, color='k', fontproperties=f, size=12)
    ax_s_5_red.text(0.0, 0.98, 'S5', verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_5_red.transAxes, color='k', fontproperties=f, size=12)

    ax_br_1_blue.text(0.05, 0.98, 'B1', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_1_blue.transAxes, color='k', fontproperties=f, size=12)
    ax_br_2_blue.text(0.05, 0.98, 'B2', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_2_blue.transAxes, color='k', fontproperties=f, size=12)
    ax_br_3_blue.text(0.05, 0.98, 'B3', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_3_blue.transAxes, color='k', fontproperties=f, size=12)

    # Save figure
    fig.savefig(taffy_extdir + 'figures_stitched_cube/line_profile_plot.eps', dpi=300, bbox_inches='tight')

    sys.exit(0)