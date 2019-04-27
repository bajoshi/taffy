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
from matplotlib.patches import Rectangle, Polygon, Arrow
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # does not have a trailing slash
taffy_dir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

sys.path.append(taffy_dir + 'codes/')
import vel_channel_map as vcm

# Define constants
lightspeed = 3e5 # 299792.458  # km/s
redshift = 0.0145
sys_vel = redshift * lightspeed
halpha_air_wav = 6562.80
oiii_air_wav = 5006.84

# this following rgb_to_hex function came from
# https://stackoverflow.com/a/214657
def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def plot_line_profiles(axesblue, axesred, bluewav, bluespec, redwav, redspec, \
    b_ylow, b_yhigh, r_ylow, r_yhigh):

    # Convert air wav to heliocentric vel
    helio_vel_arr_oiii = ((bluewav - oiii_air_wav) / oiii_air_wav) * lightspeed
    helio_vel_arr_halpha = ((redwav - halpha_air_wav) / halpha_air_wav) * lightspeed

    # convert to velocities wrt to systemic
    vel_arr_oiii = helio_vel_arr_oiii - sys_vel
    vel_arr_halpha = helio_vel_arr_halpha - sys_vel

    # plot spectra
    axesblue.plot(vel_arr_oiii, bluespec, color=myblue)
    axesred.plot(vel_arr_halpha, redspec, color=myred)

    # Set limits
    b_low_idx = np.argmin(abs(bluewav - 5070))
    b_high_idx = np.argmin(abs(bluewav - 5090))
    axesblue.set_xlim(vel_arr_oiii[b_low_idx], vel_arr_oiii[b_high_idx])

    r_low_idx = np.argmin(abs(redwav - 6645))
    r_high_idx = np.argmin(abs(redwav - 6675))
    axesred.set_xlim(vel_arr_halpha[r_low_idx], vel_arr_halpha[r_high_idx])

    axesblue.set_ylim(b_ylow, b_yhigh)
    axesred.set_ylim(r_ylow, r_yhigh)

    # --- commands to make the plot pretty --- #
    # do not show the spines except at the bottom and left
    axesblue.spines["top"].set_visible(False)
    axesblue.spines["right"].set_visible(False)
    axesblue.spines["bottom"].set_visible(True)
    axesblue.spines["left"].set_visible(True)

    axesred.spines["top"].set_visible(False)
    axesred.spines["right"].set_visible(False)
    axesred.spines["bottom"].set_visible(True)
    axesred.spines["left"].set_visible(True)

    # Increase y axis ticklabel size
    # First convert to a numpy array of ints because it is float by default
    blue_ax_ytickslist = np.array(axesblue.get_yticks().tolist(), dtype=np.int)
    red_ax_ytickslist = np.array(axesred.get_yticks().tolist(), dtype=np.int)
    axesblue.set_yticklabels(blue_ax_ytickslist, size=12, rotation=20)
    axesred.set_yticklabels(red_ax_ytickslist, size=12, rotation=20)

    # do not show the velocity axis ticklabels 
    # instead we will just show the velocity ticks
    # in the bottom panels which will be done in __main__
    axesblue.get_xaxis().set_ticklabels([])
    axesred.get_xaxis().set_ticklabels([])

    # Turn on minor ticks
    axesblue.minorticks_on()
    axesred.minorticks_on()

    # turn off ticks on the right
    axesblue.tick_params(axis='y', which='both', right=False)
    axesred.tick_params(axis='y', which='both', right=False)

    # Draw a vertical line at the systemic recessional velocity
    axesblue.axvline(x=0.0, ls='--', color='grey')
    axesred.axvline(x=0.0, ls='--', color='grey')

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

    # ----------- Assign gridspec locations ----------- #
    # define axes using above grid
    ax = fig.add_subplot(gs[:36,24:60], projection=wcs_sdss)

    """
    These rightward increasing looking zorder values for hte subplots are to make sure that
    the axes of a graph do not get cut off by the plot to its left.
    Try making the plot without the zorder in the subplots and it'll be clear (check B1 in particular).
    I'm not entirely sure why the problem arises in the first place but this fixes it.
    """

    # Bridge
    ax_br_1_blue = fig.add_subplot(gs[36:48,24:36], zorder=5)
    ax_br_1_red  = fig.add_subplot(gs[48:,24:36], zorder=5)

    ax_br_2_blue = fig.add_subplot(gs[36:48,36:48], zorder=6)
    ax_br_2_red  = fig.add_subplot(gs[48:,36:48], zorder=6)

    ax_br_3_blue = fig.add_subplot(gs[36:48,48:60], zorder=7)
    ax_br_3_red  = fig.add_subplot(gs[48:,48:60], zorder=7)

    # North galaxy
    ax_n_1_blue = fig.add_subplot(gs[:12,:12], zorder=-2)
    ax_n_2_blue = fig.add_subplot(gs[12:24,:12], zorder=-2)
    ax_n_3_blue = fig.add_subplot(gs[24:36,:12], zorder=-2)
    ax_n_4_blue = fig.add_subplot(gs[36:48,:12])
    ax_n_5_blue = fig.add_subplot(gs[48:,:12])

    ax_n_1_red  = fig.add_subplot(gs[:12,12:24], zorder=-1)
    ax_n_2_red  = fig.add_subplot(gs[12:24,12:24], zorder=-1)
    ax_n_3_red  = fig.add_subplot(gs[24:36,12:24], zorder=-1)
    ax_n_4_red  = fig.add_subplot(gs[36:48,12:24])
    ax_n_5_red  = fig.add_subplot(gs[48:,12:24])

    """
    # Orginial Locations
    # define axes using above grid
    ax = fig.add_subplot(gs[:36,:36], projection=wcs_sdss)

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
    """

    # South galaxy
    ax_s_1_blue = fig.add_subplot(gs[:12,60:72])
    ax_s_2_blue = fig.add_subplot(gs[12:24,60:72])
    ax_s_3_blue = fig.add_subplot(gs[24:36,60:72])
    ax_s_4_blue = fig.add_subplot(gs[36:48,60:72], zorder=8)
    ax_s_5_blue = fig.add_subplot(gs[48:,60:72], zorder=8)

    ax_s_1_red  = fig.add_subplot(gs[:12,72:])
    ax_s_2_red  = fig.add_subplot(gs[12:24,72:])
    ax_s_3_red  = fig.add_subplot(gs[24:36,72:])
    ax_s_4_red  = fig.add_subplot(gs[36:48,72:], zorder=9)
    ax_s_5_red  = fig.add_subplot(gs[48:,72:], zorder=9)

    # ------------ Plot galaxy ------------ #
    # SDSS g band image with proper stretch
    # set axis labels
    ax.set_xlabel('Right Ascension (J2000)', fontsize=8)
    ax.set_ylabel('Declination (J2000)', fontsize=8)

    # get the correct sized image and normalization and plot
    norm = ImageNormalize(sdss_i[0].data, stretch=LogStretch())
    orig_cmap = mpl.cm.Greys
    shifted_cmap = vcm.shiftedColorMap(orig_cmap, midpoint=0.5, name='shifted')
    im = ax.imshow(sdss_i[0].data, origin='lower', cmap=shifted_cmap, vmin=0.05, vmax=8, norm=norm)

    ax.set_autoscale_on(False)

    lon = ax.coords[0]
    lat = ax.coords[1]

    #lon.set_ticks_visible(False)
    #lon.set_ticklabel_visible(False)
    #lat.set_ticks_visible(False)
    #lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')

    ax.coords.frame.set_color('k')
    ax.grid(color='gray', ls='dashed', lw=0.7)

    # Add scale and compass
    # Bold text
    f = FontProperties()
    f.set_weight('bold')

    scalebar_x = [0.435, 0.435]
    scalebar_y = [23.469, 23.469 + 60/3600]
    ax.plot(scalebar_x, scalebar_y, color='k', lw=4.0, transform=ax.get_transform('fk5'), zorder=5)
    ax.text(0.1, 0.1, '1 arcmin = 18 kpc', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)

    # Add polygon patches for galaxies and bridge.
    # I copied these directly from the ds9 regions file
    # There is a file called all_regions.reg in $HOME/Desktop/ipac/taffy_lzifu/
    # which has all of the regions in the required format.
    # These have to be in degrees as [dec, ra], i.e. a list of coordinate lists .
    north_poly_points = [[0.4221999963,23.50396722],[0.4134583314,23.50394611],[0.412683328,23.502875],\
    	[0.4127041658,23.50115528],[0.4168249925,23.49923944],[0.4199708303,23.49723361],[0.4259624958,23.49127472],\
    	[0.4286958377,23.48904167],[0.4291166623,23.49004194],[0.4313166618,23.49288667],[0.4291875045,23.49689222]]
    south_poly_points = [[0.4071291606,23.49766389],[0.4027499994,23.49619722],[0.4031333288,23.48922167],\
        [0.4060125033,23.48189472],[0.4078666687,23.47872972],[0.4109458288,23.47614944],[0.4158541679,23.47630528],\
        [0.4201999982,23.47791056],[0.4217333317,23.48049],[0.4214000066,23.48310944],[0.41794583,23.48275694],\
        [0.4133583387,23.48707361],[0.4114708265,23.48932528],[0.409041659,23.49400972]]
    bridge_poly_points = [[0.4071291606,23.49758639],[0.4114374955,23.4893825],[0.4133458296,23.48707306],\
        [0.4179583391,23.48273472],[0.4213874976,23.48308472],[0.4228458405,23.48838778],[0.4258291721,23.49137333],\
        [0.4199666659,23.49723472],[0.4167291641,23.49936111],[0.4127541701,23.50118361]]
    bridge_north_poly_points = [[0.4127250036,23.50113444],[0.4071333408,23.49758472],[0.4080041726,23.49578778],\
        [0.4114833355,23.48936833],[0.4143750032,23.49176806],[0.4193958282,23.49761556],[0.4167374929,23.49928694]]
    north_west_poly_points = [[0.4221999963,23.50395139],[0.4134458383,23.50394639],[0.4127000014,23.50286139],\
        [0.4127333323,23.50115639],[0.4167583307,23.49929667],[0.4193874995,23.49760556],[0.424620835,23.50152333]]
    south_nuc_poly_points = [[0.4109875043,23.48449194],[0.4090958277,23.48535917],[0.4080208302,23.48467111],\
        [0.4084750017,23.48272667],[0.4094875018,23.48191917],[0.4110541662,23.48191917],[0.4116416613,23.48263722]]

    pn = Polygon(np.array(north_poly_points), edgecolor='forestgreen', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(pn)
    ps = Polygon(np.array(south_poly_points), edgecolor='midnightblue', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(ps)
    pb = Polygon(np.array(bridge_poly_points), edgecolor='maroon', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(pb)
    pbn = Polygon(np.array(bridge_north_poly_points), edgecolor='darkorange', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(pbn)
    pnw = Polygon(np.array(north_west_poly_points), edgecolor='darkorchid', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(pnw)
    psnuc = Polygon(np.array(south_nuc_poly_points), edgecolor='midnightblue', facecolor='None', \
        closed=True, transform=ax.get_transform('fk5'), lw=3.0)
    ax.add_patch(psnuc)

    # Add text to figure to indicate region name
    ax.text(0.18, 0.7, 'N1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.23, 0.75, 'N2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.27, 0.8, 'N3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.32, 0.85, 'N4', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.38, 0.9, 'N5', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)

    ax.text(0.38, 0.63, 'B1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)
    ax.text(0.54, 0.79, 'B2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)
    ax.text(0.47, 0.52, 'B3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)

    ax.text(0.54, 0.3, 'S1', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.65, 0.4, 'S2', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.71, 0.49, 'S3', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.74, 0.6, 'S4', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)
    ax.text(0.73, 0.68, 'S5', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='w', fontproperties=f, size=13)

    #ax.text(0.02, 0.1, r'$\mathrm{[OIII]\lambda 5007}$' + '  ' + r'$\mathrm{Velocity\ scale\, (km/s):[-575,+623]}$', \
    #    verticalalignment='top', horizontalalignment='left', \
    #    transform=ax.transAxes, color=myblue, size=16)
    #ax.text(0.02, 0.06, r'$\mathrm{H\alpha}$' + '               ' + r'$\mathrm{Velocity\ scale\, (km/s):[-602,+770]}$', \
    #    verticalalignment='top', horizontalalignment='left', \
    #    transform=ax.transAxes, color=myred, size=16)

    # Add text for galaxies and bridge
    #r'$\mathrm{Taffy}$' + '-' + r'$\mathrm{N}$' + '\n' + r'$\mathrm{UGC\ 12915}$', \
    ax.text(0.02, 0.96, 'Taffy-N' + '\n' + 'UGC 12915', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=20)
    ax.text(0.67, 0.3, 'Taffy-S' + '\n' + 'UGC 12914', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=20)
    ax.text(0.59, 0.85, 'Taffy Bridge', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=20)

    # ------------------------ PLOT THE REGIONS ------------------------ #
    reg_file = open(taffy_extdir + 'rawspectra_for_paperplot/regions_for_line_profiles.reg')
    circ_size = 3.776 / 3600 

    for ln_count, ln in enumerate(reg_file.readlines()):
        if ln_count < 3:
            continue
        dec = float(ln.split(',')[1])
        ra = float(ln.split(',')[0].split('circle(')[1])

        reg_color = 'white'
        if 'color=blue' in ln:
            reg_color = 'k'  # Now coloring the bridge circles maroon to be consistent with teh BPT plots

        region = SphericalCircle((ra * u.deg, dec * u.deg), circ_size * u.degree, \
            edgecolor=reg_color, facecolor='none', transform=ax.get_transform('fk5'), lw=2.0)
        ax.add_patch(region)

    # ------------------------ PLOT THE SPECTRA ------------------------ #
    ax_n_1_blue, ax_n_1_red = plot_line_profiles(ax_n_1_blue, ax_n_1_red, north_majax_1_blue['wav'], north_majax_1_blue['flux'], north_majax_1_red['wav'], north_majax_1_red['flux'], 55,90, 40,175)
    ax_n_2_blue, ax_n_2_red = plot_line_profiles(ax_n_2_blue, ax_n_2_red, north_majax_2_blue['wav'], north_majax_2_blue['flux'], north_majax_2_red['wav'], north_majax_2_red['flux'], 40,80, 60,320)
    ax_n_3_blue, ax_n_3_red = plot_line_profiles(ax_n_3_blue, ax_n_3_red, north_majax_3_blue['wav'], north_majax_3_blue['flux'], north_majax_3_red['wav'], north_majax_3_red['flux'], 20,48, 20,235)
    ax_n_4_blue, ax_n_4_red = plot_line_profiles(ax_n_4_blue, ax_n_4_red, north_majax_4_blue['wav'], north_majax_4_blue['flux'], north_majax_4_red['wav'], north_majax_4_red['flux'], 15,50, 15,95)
    ax_n_5_blue, ax_n_5_red = plot_line_profiles(ax_n_5_blue, ax_n_5_red, north_majax_5_blue['wav'], north_majax_5_blue['flux'], north_majax_5_red['wav'], north_majax_5_red['flux'], 20,75, 0,175)

    ax_s_1_blue, ax_s_1_red = plot_line_profiles(ax_s_1_blue, ax_s_1_red, south_majax_1_blue['wav'], south_majax_1_blue['flux'], south_majax_1_red['wav'], south_majax_1_red['flux'], 24,56, 30,340)
    ax_s_2_blue, ax_s_2_red = plot_line_profiles(ax_s_2_blue, ax_s_2_red, south_majax_2_blue['wav'], south_majax_2_blue['flux'], south_majax_2_red['wav'], south_majax_2_red['flux'], 105,160, 140,260)
    ax_s_3_blue, ax_s_3_red = plot_line_profiles(ax_s_3_blue, ax_s_3_red, south_majax_3_blue['wav'], south_majax_3_blue['flux'], south_majax_3_red['wav'], south_majax_3_red['flux'], 25,63, 40,275)
    ax_s_4_blue, ax_s_4_red = plot_line_profiles(ax_s_4_blue, ax_s_4_red, south_majax_4_blue['wav'], south_majax_4_blue['flux'], south_majax_4_red['wav'], south_majax_4_red['flux'], 40,120, 40,340)
    ax_s_5_blue, ax_s_5_red = plot_line_profiles(ax_s_5_blue, ax_s_5_red, south_majax_5_blue['wav'], south_majax_5_blue['flux'], south_majax_5_red['wav'], south_majax_5_red['flux'], 25,75, 20,130)

    ax_br_1_blue, ax_br_1_red = plot_line_profiles(ax_br_1_blue, ax_br_1_red, bridge_1_blue['wav'], bridge_1_blue['flux'], bridge_1_red['wav'], bridge_1_red['flux'], -7,75, -15,200)
    ax_br_2_blue, ax_br_2_red = plot_line_profiles(ax_br_2_blue, ax_br_2_red, bridge_2_blue['wav'], bridge_2_blue['flux'], bridge_2_red['wav'], bridge_2_red['flux'], -10,10, -10,45)
    ax_br_3_blue, ax_br_3_red = plot_line_profiles(ax_br_3_blue, ax_br_3_red, bridge_3_blue['wav'], bridge_3_blue['flux'], bridge_3_red['wav'], bridge_3_red['flux'], -10,25, -10,50)

    # Velocity axis labels
    blue_ax_xtickslist_n5 = np.array(ax_n_5_blue.get_xticks().tolist(), dtype=np.int)
    red_ax_xtickslist_n5 = np.array(ax_n_5_red.get_xticks().tolist(), dtype=np.int)
    ax_n_5_blue.set_xticklabels(blue_ax_xtickslist_n5, size=12)
    ax_n_5_red.set_xticklabels(red_ax_xtickslist_n5, size=12)

    blue_ax_xtickslist_s5 = np.array(ax_s_5_blue.get_xticks().tolist(), dtype=np.int)
    red_ax_xtickslist_s5 = np.array(ax_s_5_red.get_xticks().tolist(), dtype=np.int)
    ax_s_5_blue.set_xticklabels(blue_ax_xtickslist_s5, size=12)
    ax_s_5_red.set_xticklabels(red_ax_xtickslist_s5, size=12)

    red_ax_xtickslist_b1 = np.array(ax_br_1_red.get_xticks().tolist(), dtype=np.int)
    red_ax_xtickslist_b2 = np.array(ax_br_2_red.get_xticks().tolist(), dtype=np.int)
    red_ax_xtickslist_b3 = np.array(ax_br_3_red.get_xticks().tolist(), dtype=np.int)
    ax_br_1_red.set_xticklabels(red_ax_xtickslist_b1, size=12)
    ax_br_2_red.set_xticklabels(red_ax_xtickslist_b2, size=12)
    ax_br_3_red.set_xticklabels(red_ax_xtickslist_b3, size=12)

    # ------------------ Text for regionname ------------------- #
    ax_n_1_blue.text(0.05, 0.75, 'N1', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_1_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_n_2_blue.text(0.05, 0.9, 'N2', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_2_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_n_3_blue.text(0.05, 0.9, 'N3', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_3_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_n_4_blue.text(0.05, 0.9, 'N4', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_4_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_n_5_blue.text(0.05, 0.9, 'N5', verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_5_blue.transAxes, color='k', fontproperties=f, size=16)

    ax_s_1_blue.text(0.05, 0.75, 'S1', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_1_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_s_2_blue.text(0.05, 0.9, 'S2', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_2_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_s_3_blue.text(0.05, 0.9, 'S3', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_3_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_s_4_blue.text(0.05, 0.9, 'S4', verticalalignment='top', horizontalalignment='left', \
       transform=ax_s_4_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_s_5_blue.text(0.05, 0.9, 'S5', verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_5_blue.transAxes, color='k', fontproperties=f, size=16)

    ax_br_1_blue.text(0.05, 0.9, 'B1', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_1_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_br_2_blue.text(0.05, 0.9, 'B2', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_2_blue.transAxes, color='k', fontproperties=f, size=16)
    ax_br_3_blue.text(0.05, 0.9, 'B3', verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_3_blue.transAxes, color='k', fontproperties=f, size=16)

    # Text for flux units
    """
    ax_n_1_blue.text(0.03, 1.06, r'$\times 10^{-18}$' + '\n' + r'$\mathrm{erg\, s^{-1}\, cm^{-2}\, \AA^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_1_blue.transAxes, color='k', size=14)
    ax_n_1_red.text(0.03, 1.06, r'$\times 10^{-18}$' + '\n' + r'$\mathrm{erg\, s^{-1}\, cm^{-2}\, \AA^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_1_red.transAxes, color='k', size=14)

    ax_s_1_blue.text(0.03, 1.06, r'$\times 10^{-18}$' + '\n' + r'$\mathrm{erg\, s^{-1}\, cm^{-2}\, \AA^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_1_blue.transAxes, color='k', size=14)
    ax_s_1_red.text(0.03, 1.06, r'$\times 10^{-18}$' + '\n' + r'$\mathrm{erg\, s^{-1}\, cm^{-2}\, \AA^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_1_red.transAxes, color='k', size=14)
    """

    # Text for [NII] lines seen in some cases
    ax_n_3_red.text(0.08, 0.44, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6548$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_3_red.transAxes, color='k', size=13)
    ax_n_4_red.text(0.07, 0.55, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6548$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_n_4_red.transAxes, color='k', size=13)

    ax_s_1_red.text(0.05, 0.43, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6548$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_1_red.transAxes, color='k', size=13)
    ax_s_3_red.text(0.75, 0.6, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6583$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_3_red.transAxes, color='k', size=13)
    ax_s_4_red.text(0.75, 0.85, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6583$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_4_red.transAxes, color='k', size=13)
    ax_s_5_red.text(0.75, 0.83, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6583$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_s_5_red.transAxes, color='k', size=13)

    ax_br_1_red.text(0.03, 0.43, r'$\mathrm{[NII]}$' + '\n' + r'$\lambda6548$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax_br_1_red.transAxes, color='k', size=13)

    # Save figure
    fig.savefig(taffy_extdir + 'figures_stitched_cube/line_profile_plot_newgrid.pdf', dpi=300, bbox_inches='tight')

    sys.exit(0)