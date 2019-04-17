from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel

import os
import sys
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
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

    # read in halpha map from the stitched cube
    halpha = fits.open(taffy_extdir + 'stitched_cube_HALPHA.fits')
    halpha_total = halpha[0].data[0]

    # Plotting
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # modify rc Params
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.sans-serif"] = ["Computer Modern Sans"]
    mpl.rcParams["text.usetex"] = False
    mpl.rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss)

    # draw contours
    x = np.arange(58)
    y = np.arange(58)
    X, Y = np.meshgrid(x,y)

    # get colormap
    colorbrewer_cm = vcm.get_colorbrewer_cm('blues')

    # plot contours
    levels = np.array([30, 50, 100, 200, 300, 400, 600, 800, 1000, 1200])
    # try smoothing the map to get smoother contours
    # define kernel
    kernel = Gaussian2DKernel(stddev=0.5)
    halpha_total = convolve(halpha_total, kernel, boundary='extend')

    c = ax.contour(X, Y, halpha_total, transform=ax.get_transform(wcs_lzifu),\
     levels=levels, cmap=colorbrewer_cm, linewidths=2.0, interpolation='None')
    #ax.clabel(c, inline=True, inline_spacing=0, fontsize=12, fmt='%1.1f', lw=4, ls='-')

    # add colorbar inside figure
    cbaxes = inset_axes(ax, width='30%', height='3%', loc=8, bbox_to_anchor=[0.02, 0.08, 1, 1], bbox_transform=ax.transAxes)
    cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels), max(levels)], orientation='horizontal')
    cb.ax.tick_params(labelsize=13)
    cb.ax.get_children()[0].set_linewidths(15.0)
    cb.ax.set_xlabel(r'$\mathbf{Total\ H\alpha \ flux\ [10^{-18}\, erg\, s^{-1}\, cm^{-2}]}$', fontsize=18)

    # add rectangle to show IFU coverage
    ifu_cover = Rectangle((181.46044-1 - 292.26768/2, 206.09294-1 - 292.62374/2), 292.26768, 292.62374,\
     edgecolor='red', facecolor='none', lw=2, ls='--')
    ax.add_patch(ifu_cover)

    # OTher formatting stuff
    lon = ax.coords[0]
    lat = ax.coords[1]

    lon.set_ticks_visible(True)
    lon.set_ticklabel_visible(True)
    lat.set_ticks_visible(True)
    lat.set_ticklabel_visible(True)
    
    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    
    lon.set_axislabel('Right Ascension', fontsize=16)
    lat.set_axislabel('Declination', fontsize=16)

    ax.coords.frame.set_color('k')
    ax.grid(color='gray', ls='dashed', lw=0.7)

    # Add text for figure panel
    f = FontProperties()
    f.set_weight('bold')
    ax.text(0.03, 0.97, '(a)', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=20)

    # save the figure
    fig.savefig(taffy_extdir + 'figures_stitched_cube/halpha_contour_smooth.pdf', dpi=300, bbox_inches='tight')

    sys.exit(0)