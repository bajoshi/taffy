from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle

import sys
import os
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

home = os.getenv('HOME')  # does not have a trailing slash
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_data_dir = '/Volumes/Bhavins_backup/ipac/Taffy_2/'

# this following rgb_to_hex function came from
# https://stackoverflow.com/a/214657
def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in g band SDSS image
    sdss_g = fits.open(home + '/Desktop/ipac/taffy_lzifu/taffy_xliners_figs_misc_data/taffy/SDSS/taffyg_sdss.fits')
    wcs = WCS(sdss_g[0].header)
    print wcs

    # read in all integrated spectra
    """
    nnuc_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/north_gal_nuc_blue_ifu.dat', dtype=None, names=['wav','flux'])
    nnuc_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/north_gal_nuc_red_ifu.dat', dtype=None, names=['wav','flux'])

    snuc_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/south_gal_nuc_blue_ifu.dat', dtype=None, names=['wav','flux'])
    snuc_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/south_gal_nuc_red_ifu.dat', dtype=None, names=['wav','flux'])

    b1_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b1_blue_ifu.dat', dtype=None, names=['wav','flux'])
    b1_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b1_red_ifu.dat', dtype=None, names=['wav','flux'])

    b2_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b2_blue_ifu.dat', dtype=None, names=['wav','flux'])
    b2_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b2_red_ifu.dat', dtype=None, names=['wav','flux'])

    b3_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b3_blue_ifu.dat', dtype=None, names=['wav','flux'])
    b3_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b3_red_ifu.dat', dtype=None, names=['wav','flux'])

    b4_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b4_blue_ifu.dat', dtype=None, names=['wav','flux'])
    b4_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/b4_red_ifu.dat', dtype=None, names=['wav','flux'])

    xgal_hii_blue = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/xgal_hii_blue_ifu.dat', dtype=None, names=['wav','flux'])
    xgal_hii_red = np.genfromtxt(taffy_dir + 'raw_spectra_for_ifu_paper_regions/xgal_hii_red_ifu.dat', dtype=None, names=['wav','flux'])
    """

    # Make figure
    fig = plt.figure(figsize=(12,12))

    # modify rc Params
    mpl.rcParams['font.family'] = 'serif'

    # define colors
    myblue = rgb_to_hex(0, 100, 180)
    myred = rgb_to_hex(214, 39, 40)  # tableau 20 red

    # create grid to plot
    gs = gridspec.GridSpec(25,42)
    gs.update(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace=30, hspace=40)

    # define axes using above grid
    ax = fig.add_subplot(gs[5:20,0:20], projection=wcs)

    ax_nnuc_b = fig.add_subplot(gs[0:5,0:10])
    ax_nnuc_r = fig.add_subplot(gs[0:5,10:20])

    ax_xgal_b = fig.add_subplot(gs[0:5,22:32])
    ax_xgal_r = fig.add_subplot(gs[0:5,32:42])

    ax_b1_b = fig.add_subplot(gs[5:10,22:32])
    ax_b1_r = fig.add_subplot(gs[5:10,32:42])

    ax_b2_b = fig.add_subplot(gs[10:15,22:32])
    ax_b2_r = fig.add_subplot(gs[10:15,32:42])

    ax_b3_b = fig.add_subplot(gs[15:20,22:32])
    ax_b3_r = fig.add_subplot(gs[15:20,32:42])

    ax_snuc_b = fig.add_subplot(gs[20:25,22:32])
    ax_snuc_r = fig.add_subplot(gs[20:25,32:42])

    ax_b4_b = fig.add_subplot(gs[20:25,0:10])
    ax_b4_r = fig.add_subplot(gs[20:25,10:20])

    # ------------ Plot galaxy ------------ #
    # SDSS g band image with proper stretch
    # set axis labels
    ax.set_xlabel('Right Ascension')
    ax.set_ylabel('Declination')

    # get the correct sized image and normalization and plot
    cutout = sdss_g[0].data
    im = ax.imshow(cutout, origin='lower', cmap='Greys_r', vmin=-0.0545, vmax=0.5)
    # [45:-80,100:-50] is hte array slice that will be zoomed in on just the galaxies 
    # but the WCS coordinates are incorrect then

    ax.set_autoscale_on(False)
    plt.show()
    sys.exit(0)

    # the format for the coordinates here is 
    # Rectangle((lower_left_x, lower_left_y), width, height, ....)
    # the -1 is to go from ds9 x,y to numpy row,col
    # i've translated the center coords that i got from ds9
    # to the coords of the lower left corner

    b1 = Rectangle((362.33361-1 - 70.707071/2, 365.47836-1 - 42.424242/2), 70.707071, 42.424242, edgecolor='red', facecolor='none', lw=2)
    b2 = Rectangle((330.75788-1 - 49.025505/2, 322.80038-1 - 42.926768/2), 49.025505, 42.926768, edgecolor='red', facecolor='none', lw=2)
    b3 = Rectangle((284.55414-1 - 84.848485/2, 280.65364-1 - 42.424242/2), 84.848485, 42.424242, edgecolor='red', facecolor='none', lw=2)
    b4 = Rectangle((256.27391-1 - 56.565657/2, 238.23011-1 - 42.424242/2), 56.565657, 42.424242, edgecolor='red', facecolor='none', lw=2)
    ax.add_patch(b1)
    ax.add_patch(b2)
    ax.add_patch(b3)
    ax.add_patch(b4)

    xgal_hii = Rectangle((283.09567-1 - 46.296717/2, 323.63346-1 - 41.088131/2), 46.296717, 41.088131,\
     edgecolor='green', facecolor='none', lw=2)
    ax.add_patch(xgal_hii)

    north_nuc = SphericalCircle((0.42523865 * u.deg, 23.496016 * u.deg), 0.0015 * u.degree, edgecolor='blue', facecolor='none',\
        transform=ax.get_transform('fk5'), lw=2)
    south_nuc = SphericalCircle((0.40989012 * u.deg, 23.483233 * u.deg), 0.0015 * u.degree, edgecolor='blue', facecolor='none',\
        transform=ax.get_transform('fk5'), lw=2)
    ax.add_patch(north_nuc)
    ax.add_patch(south_nuc)

    # ------------ Plot north nuclear region spectrum ------------ #

    ax_nnuc_b.plot(nnuc_blue['wav'], nnuc_blue['flux'], color=myblue)
    ax_nnuc_r.plot(nnuc_red['wav'], nnuc_red['flux'], color=myred)

    ax_nnuc_b.spines["top"].set_visible(False)
    ax_nnuc_b.spines["right"].set_visible(False)
    ax_nnuc_r.spines["top"].set_visible(False)
    ax_nnuc_r.spines["right"].set_visible(False)

    ax_nnuc_b.get_xaxis().tick_bottom()
    ax_nnuc_b.get_yaxis().tick_left()
    ax_nnuc_r.get_xaxis().tick_bottom()
    ax_nnuc_r.get_yaxis().tick_left()

    ax_nnuc_b.set_ylim(800,1700)
    ax_nnuc_r.set_ylim(1200,4700)
    ax_nnuc_b.set_xlim(4800,5300)
    ax_nnuc_r.set_xlim(6370,6750)

    ax_nnuc_b.set_xticklabels(ax_nnuc_b.get_xticks().tolist(), size=6, rotation=30)
    ax_nnuc_r.set_xticklabels(ax_nnuc_r.get_xticks().tolist(), size=6, rotation=30)
    ax_nnuc_b.set_yticklabels(ax_nnuc_b.get_yticks().tolist(), size=6, rotation=30)
    ax_nnuc_r.set_yticklabels(ax_nnuc_r.get_yticks().tolist(), size=6, rotation=30)

    # I couldn't set the pad in the above lines, so I had to do it separately
    ax_nnuc_b.tick_params('both', which='major', pad=1)
    ax_nnuc_r.tick_params('both', which='major', pad=1)

    #ax.minorticks_on()
    #ax.tick_params('both', width=1, length=3, which='minor')
    #ax.tick_params('both', width=1, length=4.7, which='major')

    # ------------ Plot south nuclear region spectrum ------------ #

    ax_snuc_b.plot(snuc_blue['wav'], snuc_blue['flux'], color=myblue)
    ax_snuc_r.plot(snuc_red['wav'], snuc_red['flux'], color=myred)

    ax_snuc_b.set_ylim(1700,3000)
    ax_snuc_r.set_ylim(2000,5000)

    ax_snuc_b.set_xlim(4800,5300)
    ax_snuc_r.set_xlim(6370,6750)

    ax_snuc_b.spines["top"].set_visible(False)
    ax_snuc_b.spines["bottom"].set_visible(False)
    ax_snuc_b.spines["right"].set_visible(False)
    ax_snuc_b.spines["left"].set_visible(False)

    ax_snuc_r.spines["top"].set_visible(False)
    ax_snuc_r.spines["bottom"].set_visible(False)
    ax_snuc_r.spines["right"].set_visible(False)
    ax_snuc_r.spines["left"].set_visible(False)

    ax_snuc_b.get_xaxis().tick_bottom()
    ax_snuc_b.get_yaxis().tick_left()

    ax_snuc_r.get_xaxis().tick_bottom()
    ax_snuc_r.get_yaxis().tick_left()

    # ------------ Plot b1 region spectrum ------------ #

    ax_b1_b.plot(b1_blue['wav'], b1_blue['flux'], color=myblue)
    ax_b1_r.plot(b1_red['wav'], b1_red['flux'], color=myred)

    ax_b1_b.set_ylim(0,1200)
    ax_b1_r.set_ylim(-300,3200)

    ax_b1_b.set_xlim(4800,5300)
    ax_b1_r.set_xlim(6370,6750)

    ax_b1_b.spines["top"].set_visible(False)
    ax_b1_b.spines["bottom"].set_visible(False)
    ax_b1_b.spines["right"].set_visible(False)
    ax_b1_b.spines["left"].set_visible(False)

    ax_b1_r.spines["top"].set_visible(False)
    ax_b1_r.spines["bottom"].set_visible(False)
    ax_b1_r.spines["right"].set_visible(False)
    ax_b1_r.spines["left"].set_visible(False)

    ax_b1_b.get_xaxis().tick_bottom()
    ax_b1_b.get_yaxis().tick_left()

    ax_b1_r.get_xaxis().tick_bottom()
    ax_b1_r.get_yaxis().tick_left()

    # ------------ Plot b2 region spectrum ------------ #

    ax_b2_b.plot(b2_blue['wav'], b2_blue['flux'], color=myblue)
    ax_b2_r.plot(b2_red['wav'], b2_red['flux'], color=myred)

    ax_b2_b.set_ylim(-500,800)
    ax_b2_r.set_ylim(-500,3000)

    ax_b2_b.set_xlim(4800,5300)
    ax_b2_r.set_xlim(6370,6750)

    ax_b2_b.spines["top"].set_visible(False)
    ax_b2_b.spines["bottom"].set_visible(False)
    ax_b2_b.spines["right"].set_visible(False)
    ax_b2_b.spines["left"].set_visible(False)

    ax_b2_r.spines["top"].set_visible(False)
    ax_b2_r.spines["bottom"].set_visible(False)
    ax_b2_r.spines["right"].set_visible(False)
    ax_b2_r.spines["left"].set_visible(False)

    ax_b2_b.get_xaxis().tick_bottom()
    ax_b2_b.get_yaxis().tick_left()

    ax_b2_r.get_xaxis().tick_bottom()
    ax_b2_r.get_yaxis().tick_left()

    # ------------ Plot b3 region spectrum ------------ #

    ax_b3_b.plot(b3_blue['wav'], b3_blue['flux'], color=myblue)
    ax_b3_r.plot(b3_red['wav'], b3_red['flux'], color=myred)

    ax_b3_b.set_ylim(-400,1000)
    ax_b3_r.set_ylim(-500,2200)

    ax_b3_b.set_xlim(4800,5300)
    ax_b3_r.set_xlim(6370,6750)

    ax_b3_b.spines["top"].set_visible(False)
    ax_b3_b.spines["bottom"].set_visible(False)
    ax_b3_b.spines["right"].set_visible(False)
    ax_b3_b.spines["left"].set_visible(False)

    ax_b3_r.spines["top"].set_visible(False)
    ax_b3_r.spines["bottom"].set_visible(False)
    ax_b3_r.spines["right"].set_visible(False)
    ax_b3_r.spines["left"].set_visible(False)

    ax_b3_b.get_xaxis().tick_bottom()
    ax_b3_b.get_yaxis().tick_left()

    ax_b3_r.get_xaxis().tick_bottom()
    ax_b3_r.get_yaxis().tick_left()

    # ------------ Plot b4 region spectrum ------------ #

    ax_b4_b.plot(b4_blue['wav'], b4_blue['flux'], color=myblue)
    ax_b4_r.plot(b4_red['wav'], b4_red['flux'], color=myred)

    ax_b4_b.set_ylim(-200,1200)
    ax_b4_r.set_ylim(-500,1200)

    ax_b4_b.set_xlim(4800,5300)
    ax_b4_r.set_xlim(6370,6750)

    ax_b4_b.spines["top"].set_visible(False)
    ax_b4_b.spines["bottom"].set_visible(False)
    ax_b4_b.spines["right"].set_visible(False)
    ax_b4_b.spines["left"].set_visible(False)

    ax_b4_r.spines["top"].set_visible(False)
    ax_b4_r.spines["bottom"].set_visible(False)
    ax_b4_r.spines["right"].set_visible(False)
    ax_b4_r.spines["left"].set_visible(False)

    ax_b4_b.get_xaxis().tick_bottom()
    ax_b4_b.get_yaxis().tick_left()

    ax_b4_r.get_xaxis().tick_bottom()
    ax_b4_r.get_yaxis().tick_left()

    # ------------ Plot X-Galactic HII region spectrum ------------ #

    ax_xgal_b.plot(xgal_hii_blue['wav'], xgal_hii_blue['flux'], color=myblue)
    ax_xgal_r.plot(xgal_hii_red['wav'], xgal_hii_red['flux'], color=myred)

    ax_xgal_b.set_ylim(-125,1400)
    ax_xgal_r.set_ylim(-700,5000)

    ax_xgal_b.set_xlim(4800,5300)
    ax_xgal_r.set_xlim(6370,6750)

    ax_xgal_b.spines["top"].set_visible(False)
    ax_xgal_b.spines["bottom"].set_visible(False)
    ax_xgal_b.spines["right"].set_visible(False)
    ax_xgal_b.spines["left"].set_visible(False)

    ax_xgal_r.spines["top"].set_visible(False)
    ax_xgal_r.spines["bottom"].set_visible(False)
    ax_xgal_r.spines["right"].set_visible(False)
    ax_xgal_r.spines["left"].set_visible(False)

    ax_xgal_b.get_xaxis().tick_bottom()
    ax_xgal_b.get_yaxis().tick_left()

    ax_xgal_r.get_xaxis().tick_bottom()
    ax_xgal_r.get_yaxis().tick_left()

    # ------------ Save figure ------------ #
    fig.savefig(taffy_dir + 'figures/taffy_intg_spec.eps', dpi=600, bbox_inches='tight')
    #plt.show()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)