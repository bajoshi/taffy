from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

import bpt_plots as bpt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'
ipac_taffy_figdir = home + "/Desktop/ipac/taffy/figures/"

def get_comp_fluxes_vdispsort():

    """
    These are values taken from ds9 for the taffy products cube that was sorted by velocity dispersion.
    """

    oiii_north_flux_comp1 = 13117.59
    oiii_north_flux_comp2 = 8725.64
    oiii_south_flux_comp1 = 11186.97
    oiii_south_flux_comp2 = 9897.07
    oiii_snuc_flux_comp1 = 1795.84
    oiii_snuc_flux_comp2 = 1941.88
    oiii_nw_flux_comp1 = 4348.22
    oiii_nw_flux_comp2 = 2224.27
    oiii_hii_flux_comp1 = 3055.85
    oiii_hii_flux_comp2 = 3653.12
    oiii_b1_flux_comp1 = 163.92
    oiii_b1_flux_comp2 = 100.39
    oiii_b2_flux_comp1 = 630.91
    oiii_b2_flux_comp2 = 436.42
    oiii_b3_flux_comp1 = 696.85
    oiii_b3_flux_comp2 = 912.21
    oiii_b4_flux_comp1 = 157.95
    oiii_b4_flux_comp2 = 152.32

    hbeta_north_flux_comp1 = 17556.37
    hbeta_north_flux_comp2 = 11665.15
    hbeta_south_flux_comp1 = 30024.34
    hbeta_south_flux_comp2 = 14455.98
    hbeta_snuc_flux_comp1 = 5452.40
    hbeta_snuc_flux_comp2 = 3144.23
    hbeta_nw_flux_comp1 = 4369.27
    hbeta_nw_flux_comp2 = 2641.42
    hbeta_hii_flux_comp1 = 2045.0
    hbeta_hii_flux_comp2 = 2801.58
    hbeta_b1_flux_comp1 = 390.61
    hbeta_b1_flux_comp2 = 413.88
    hbeta_b2_flux_comp1 = 613.06
    hbeta_b2_flux_comp2 = 626.42
    hbeta_b3_flux_comp1 = 1230.62
    hbeta_b3_flux_comp2 = 745.55
    hbeta_b4_flux_comp1 = 457.34
    hbeta_b4_flux_comp2 = 231.37

    halpha_north_flux_comp1 = 74997.6
    halpha_north_flux_comp2 = 88537.91
    halpha_south_flux_comp1 = 106122.27
    halpha_south_flux_comp2 = 59624.53
    halpha_snuc_flux_comp1 = 15957.03
    halpha_snuc_flux_comp2 = 10202.88
    halpha_nw_flux_comp1 = 18480.32
    halpha_nw_flux_comp2 = 15097.33
    halpha_hii_flux_comp1 = 11739.45
    halpha_hii_flux_comp2 = 17380.49
    halpha_b1_flux_comp1 = 1631.40
    halpha_b1_flux_comp2 = 2364.26
    halpha_b2_flux_comp1 = 5077.37
    halpha_b2_flux_comp2 = 3493.85
    halpha_b3_flux_comp1 = 4630.74
    halpha_b3_flux_comp2 = 4509.51
    halpha_b4_flux_comp1 = 1503.37
    halpha_b4_flux_comp2 = 184.79

    nii_north_flux_comp1 = 33994.54
    nii_north_flux_comp2 = 41761.1
    nii_south_flux_comp1 = 59953.46
    nii_south_flux_comp2 = 36410.42
    nii_snuc_flux_comp1 = 9289.65
    nii_snuc_flux_comp2 = 8673.91
    nii_nw_flux_comp1 = 7049.91
    nii_nw_flux_comp2 = 5652.74
    nii_hii_flux_comp1 = 4567.76
    nii_hii_flux_comp2 = 8168.01
    nii_b1_flux_comp1 = 1108.52
    nii_b1_flux_comp2 = 1378.92
    nii_b2_flux_comp1 = 3655.37
    nii_b2_flux_comp2 = 2858.14
    nii_b3_flux_comp1 = 2404.76
    nii_b3_flux_comp2 = 2578.40
    nii_b4_flux_comp1 = 810.75
    nii_b4_flux_comp2 = 184.74

    oi_north_flux_comp1 = 4680.97
    oi_north_flux_comp2 = 9266.8
    oi_south_flux_comp1 = 6538.53
    oi_south_flux_comp2 = 6413.16
    oi_snuc_flux_comp1 = 1262.49
    oi_snuc_flux_comp2 = 2227.15
    oi_nw_flux_comp1 = 837.02
    oi_nw_flux_comp2 = 2945.78
    oi_hii_flux_comp1 = 887.21
    oi_hii_flux_comp2 = 2181.4
    oi_b1_flux_comp1 = 144.20
    oi_b1_flux_comp2 = 317.19
    oi_b2_flux_comp1 = 679.55
    oi_b2_flux_comp2 = 631.39
    oi_b3_flux_comp1 = 404.89
    oi_b3_flux_comp2 = 858.04
    oi_b4_flux_comp1 = 336.7
    oi_b4_flux_comp2 = 180.92

    sii6716_north_flux_comp1 = 20638.0
    sii6716_north_flux_comp2 = 25331.97
    sii6716_south_flux_comp1 = 27721.44
    sii6716_south_flux_comp2 = 25614.13
    sii6716_snuc_flux_comp1 = 3662.66
    sii6716_snuc_flux_comp2 = 8037.55
    sii6716_nw_flux_comp1 = 4957.45
    sii6716_nw_flux_comp2 = 4194.68
    sii6716_hii_flux_comp1 = 3339.78
    sii6716_hii_flux_comp2 = 5757.96
    sii6716_b1_flux_comp1 = 584.1
    sii6716_b1_flux_comp2 = 949.91
    sii6716_b2_flux_comp1 = 1988.14
    sii6716_b2_flux_comp2 = 1463.57
    sii6716_b3_flux_comp1 = 1574.78
    sii6716_b3_flux_comp2 = 2272.11
    sii6716_b4_flux_comp1 = 612.57
    sii6716_b4_flux_comp2 = 382.42

    sii6731_north_flux_comp1 = 13510.36
    sii6731_north_flux_comp2 = 16319.83
    sii6731_south_flux_comp1 = 22720.01
    sii6731_south_flux_comp2 = 29702.82
    sii6731_snuc_flux_comp1 = 2373.99
    sii6731_snuc_flux_comp2 = 8448.03
    sii6731_nw_flux_comp1 = 1924.52
    sii6731_nw_flux_comp2 = 4585.32
    sii6731_hii_flux_comp1 = 2128.57
    sii6731_hii_flux_comp2 = 3201.94
    sii6731_b1_flux_comp1 = 1807.68
    sii6731_b1_flux_comp2 = 1150.10
    sii6731_b2_flux_comp1 = 1373.99
    sii6731_b2_flux_comp2 = 1941.77
    sii6731_b3_flux_comp1 = 2138.45
    sii6731_b3_flux_comp2 = 1419.12
    sii6731_b4_flux_comp1 = 464.86
    sii6731_b4_flux_comp2 = 353.07

    # ------------------------------------------------------ # 

    return oiii_north_flux_comp1, oiii_north_flux_comp2, oiii_south_flux_comp1, oiii_south_flux_comp2, \
    oiii_snuc_flux_comp1, oiii_snuc_flux_comp2, oiii_nw_flux_comp1, oiii_nw_flux_comp2, oiii_hii_flux_comp1, \
    oiii_hii_flux_comp2, oiii_b1_flux_comp1, oiii_b1_flux_comp2, oiii_b2_flux_comp1, oiii_b2_flux_comp2, \
    oiii_b3_flux_comp1, oiii_b3_flux_comp2, oiii_b4_flux_comp1, oiii_b4_flux_comp2, \
    hbeta_north_flux_comp1, hbeta_north_flux_comp2, hbeta_south_flux_comp1, hbeta_south_flux_comp2, \
    hbeta_snuc_flux_comp1, hbeta_snuc_flux_comp2, hbeta_nw_flux_comp1, hbeta_nw_flux_comp2, hbeta_hii_flux_comp1, \
    hbeta_hii_flux_comp2, hbeta_b1_flux_comp1, hbeta_b1_flux_comp2, hbeta_b2_flux_comp1, hbeta_b2_flux_comp2, \
    hbeta_b3_flux_comp1, hbeta_b3_flux_comp2, hbeta_b4_flux_comp1, hbeta_b4_flux_comp2, \
    halpha_north_flux_comp1, halpha_north_flux_comp2, halpha_south_flux_comp1, halpha_south_flux_comp2, \
    halpha_snuc_flux_comp1, halpha_snuc_flux_comp2, halpha_nw_flux_comp1, halpha_nw_flux_comp2, \
    halpha_hii_flux_comp1, halpha_hii_flux_comp2, halpha_b1_flux_comp1, halpha_b1_flux_comp2, \
    halpha_b2_flux_comp1, halpha_b2_flux_comp2, halpha_b3_flux_comp1, halpha_b3_flux_comp2, \
    halpha_b4_flux_comp1, halpha_b4_flux_comp2, \
    nii_north_flux_comp1, nii_north_flux_comp2, nii_south_flux_comp1, nii_south_flux_comp2, \
    nii_snuc_flux_comp1, nii_snuc_flux_comp2, nii_nw_flux_comp1, nii_nw_flux_comp2, nii_hii_flux_comp1, \
    nii_hii_flux_comp2, nii_b1_flux_comp1, nii_b1_flux_comp2, nii_b2_flux_comp1, nii_b2_flux_comp2, \
    nii_b3_flux_comp1, nii_b3_flux_comp2, nii_b4_flux_comp1, nii_b4_flux_comp2,
    oi_north_flux_comp1, oi_north_flux_comp2, oi_south_flux_comp1, oi_south_flux_comp2, oi_snuc_flux_comp1, \
    oi_snuc_flux_comp2, oi_nw_flux_comp1, oi_nw_flux_comp2, oi_hii_flux_comp1, oi_hii_flux_comp2, \
    oi_b1_flux_comp1, oi_b1_flux_comp2, oi_b2_flux_comp1, oi_b2_flux_comp2, oi_b3_flux_comp1, \
    oi_b3_flux_comp2, oi_b4_flux_comp1, oi_b4_flux_comp2, \
    sii6716_north_flux_comp1, sii6716_north_flux_comp2, sii6716_south_flux_comp1, sii6716_south_flux_comp2, \
    sii6716_snuc_flux_comp1, sii6716_snuc_flux_comp2, sii6716_nw_flux_comp1, sii6716_nw_flux_comp2, \
    sii6716_hii_flux_comp1, sii6716_hii_flux_comp2, sii6716_b1_flux_comp1, sii6716_b1_flux_comp2, \
    sii6716_b2_flux_comp1, sii6716_b2_flux_comp2, sii6716_b3_flux_comp1, sii6716_b3_flux_comp2, \
    sii6716_b4_flux_comp1, sii6716_b4_flux_comp2, \
    sii6731_north_flux_comp1, sii6731_north_flux_comp2, sii6731_south_flux_comp1, sii6731_south_flux_comp2, \
    sii6731_snuc_flux_comp1, sii6731_snuc_flux_comp2, sii6731_nw_flux_comp1, sii6731_nw_flux_comp2, \
    sii6731_hii_flux_comp1, sii6731_hii_flux_comp2, sii6731_b1_flux_comp1, sii6731_b1_flux_comp2, \
    sii6731_b2_flux_comp1, sii6731_b2_flux_comp2, sii6731_b3_flux_comp1, sii6731_b3_flux_comp2, \
    sii6731_b4_flux_comp1, sii6731_b4_flux_comp2

def get_pn(region_list):

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    return pg.Polygon(pn_list)

def plotbpt(plottype, vel_comp, xarr_n, yarr_n, xarr_s, yarr_s, xarr_snuc, yarr_snuc, xarr_nw, yarr_nw, xarr_hii, yarr_hii,
 xarr_b1, yarr_b1, xarr_b2, yarr_b2, xarr_b3, yarr_b3, xarr_b4, yarr_b4):
    """
    All of the BPT classifications are taken from Kewley et al 2006, MNRAS, 372, 961
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    ax.plot(xarr_n, yarr_n, 'o', markersize=6, markeredgecolor='None', label='north')
    ax.plot(xarr_s, yarr_s, 'o', markersize=6, markeredgecolor='None', label='south')
    ax.plot(xarr_snuc, yarr_snuc, 'o', markersize=6, markeredgecolor='None', label='snuc')
    ax.plot(xarr_nw, yarr_nw, 'o', markersize=6, markeredgecolor='None', label='nw')
    ax.plot(xarr_hii, yarr_hii, 'o', markersize=6, markeredgecolor='None', label='hii')
    ax.plot(xarr_b1, yarr_b1, 'd', markersize=9, markeredgecolor='None', label='b1')
    ax.plot(xarr_b2, yarr_b2, 'd', markersize=9, markeredgecolor='None', label='b2')
    ax.plot(xarr_b3, yarr_b3, 'd', markersize=9, markeredgecolor='None', label='b3')
    ax.plot(xarr_b4, yarr_b4, 'd', markersize=9, markeredgecolor='None', label='b4')

    if plottype == 'nii':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
        y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

        ax.plot(np.arange(-1, 0, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.set_xlim(-1,0.3)
        ax.set_ylim(-1,1)

        # labels
        agnbox = TextArea('AGN', textprops=dict(color='k', size=16))
        anc_agnbox = AnchoredOffsetbox(loc=2, child=agnbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.57, 0.93),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_agnbox) 

        compbox = TextArea('HII-AGN Composite', textprops=dict(color='k', size=16))
        anc_compbox = AnchoredOffsetbox(loc=2, child=compbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.55, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_compbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.22, 0.3),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    elif plottype == 'oi':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[OI]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.33 + 0.73 / (np.arange(-2.5, -0.8, 0.01) + 0.59)
        y_liner_seyfert_line = 1.30 + 1.18 * np.arange(-1.1, 0, 0.01)

        ax.plot(np.arange(-2.5, -0.8, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-1.1, 0, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.set_xlim(-2.0,0)
        ax.set_ylim(-1,1)

        # labels
        seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
        anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.35, 0.93),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_seyfertbox) 

        linerbox = TextArea('LINER', textprops=dict(color='k', size=16))
        anc_linerbox = AnchoredOffsetbox(loc=2, child=linerbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.8, 0.45),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_linerbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.2, 0.2),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    elif plottype == 'sii':
        ax.set_xlabel(r'$\mathrm{log\left( \frac{[SII]}{H\alpha} \right)}$', fontsize=15)

        y_agn_hii_line = 1.3 + 0.72 / (np.arange(-1, 0.1, 0.01) - 0.32)
        y_liner_seyfert_line = 0.76 + 1.89 * np.arange(-0.3, 1, 0.01)

        ax.plot(np.arange(-1, 0.1, 0.01), y_agn_hii_line, '-', color='k')
        ax.plot(np.arange(-0.3, 1, 0.01), y_liner_seyfert_line, '--', color='k')

        ax.set_xlim(-1,0.5)
        ax.set_ylim(-1,1)

        # labels
        seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
        anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.35, 0.93),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_seyfertbox) 

        linerbox = TextArea('LINER', textprops=dict(color='k', size=16))
        anc_linerbox = AnchoredOffsetbox(loc=2, child=linerbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.75, 0.45),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_linerbox) 

        hiibox = TextArea('HII', textprops=dict(color='k', size=16))
        anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.22, 0.3),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_hiibox)

    ax.legend(loc=0, numpoints=1, prop={'size':10})

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(ipac_taffy_figdir + 'bpt_regions_' + plottype + '_comp' + vel_comp + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()
    
    return None

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'old_cube_2_comp_velsort.fits')

    # assign line arrays for each component and their errors
    # -------------- component 1 -------------- #
    halpha_comp1 = h['HALPHA'].data[1]
    hbeta_comp1 = h['HBETA'].data[1]
    nii6583_comp1 = h['NII6583'].data[1]
    oiii5007_comp1 = h['OIII5007'].data[1]
    oi6300_comp1 = h['OI6300'].data[1]
    oi6364_comp1 = h['OI6364'].data[1]
    sii6716_comp1 = h['SII6716'].data[1]
    sii6731_comp1 = h['SII6731'].data[1]

    halpha_err_comp1 = h['HALPHA_ERR'].data[1]
    hbeta_err_comp1 = h['HBETA_ERR'].data[1]
    nii6583_err_comp1 = h['NII6583_ERR'].data[1]
    oiii5007_err_comp1 = h['OIII5007_ERR'].data[1]
    oi6300_err_comp1 = h['OI6300_ERR'].data[1]
    oi6364_err_comp1 = h['OI6364_ERR'].data[1]
    sii6716_err_comp1 = h['SII6716_ERR'].data[1]
    sii6731_err_comp1 = h['SII6731_ERR'].data[1]

    # -------------- component 2 -------------- #
    halpha_comp2 = h['HALPHA'].data[2]
    hbeta_comp2 = h['HBETA'].data[2]
    nii6583_comp2 = h['NII6583'].data[2]
    oiii5007_comp2 = h['OIII5007'].data[2]
    oi6300_comp2 = h['OI6300'].data[2]
    oi6364_comp2 = h['OI6364'].data[2]
    sii6716_comp2 = h['SII6716'].data[2]
    sii6731_comp2 = h['SII6731'].data[2]

    halpha_err_comp2 = h['HALPHA_ERR'].data[2]
    hbeta_err_comp2 = h['HBETA_ERR'].data[2]
    nii6583_err_comp2 = h['NII6583_ERR'].data[2]
    oiii5007_err_comp2 = h['OIII5007_ERR'].data[2]
    oi6300_err_comp2 = h['OI6300_ERR'].data[2]
    oi6364_err_comp2 = h['OI6364_ERR'].data[2]
    sii6716_err_comp2 = h['SII6716_ERR'].data[2]
    sii6731_err_comp2 = h['SII6731_ERR'].data[2]
    
    # add lines which are doublets
    sii_comp1 = sii6716_comp1 + sii6731_comp1
    sii_err_comp1 = np.sqrt((sii6716_err_comp1)**2 + (sii6731_err_comp1)**2)

    sii_comp2 = sii6716_comp2 + sii6731_comp2
    sii_err_comp2 = np.sqrt((sii6716_err_comp2)**2 + (sii6731_err_comp2)**2)

    # apply sig and baseline cuts
    # sig cut for comp 1
    nii_halpha_withcut_comp1, oiii_hbeta_withcut_comp1, oi_halpha_withcut_comp1, sii_halpha_withcut_comp1,\
    halpha_withcut_comp1, hbeta_withcut_comp1, oiii5007_withcut_comp1, oi6300_withcut_comp1, nii6583_withcut_comp1, sii_withcut_comp1 = \
    bpt.get_arr_withsigcut(3, halpha_comp1, halpha_err_comp1, hbeta_comp1, hbeta_err_comp1, oiii5007_comp1, oiii5007_err_comp1,\
    nii6583_comp1, nii6583_err_comp1, oi6300_comp1, oi6300_err_comp1, sii_comp1, sii_err_comp1, (58,58))
    
    # sig cut for comp 2
    nii_halpha_withcut_comp2, oiii_hbeta_withcut_comp2, oi_halpha_withcut_comp2, sii_halpha_withcut_comp2,\
    halpha_withcut_comp2, hbeta_withcut_comp2, oiii5007_withcut_comp2, oi6300_withcut_comp2, nii6583_withcut_comp2, sii_withcut_comp2 = \
    bpt.get_arr_withsigcut(3, halpha_comp2, halpha_err_comp2, hbeta_comp2, hbeta_err_comp2, oiii5007_comp2, oiii5007_err_comp2,\
    nii6583_comp2, nii6583_err_comp2, oi6300_comp2, oi6300_err_comp2, sii_comp2, sii_err_comp2, (58,58))

    # get region mask for region defined first in ds9
    # see process to do this detailed in the comments in the bpt_plots.py code.
    region_file = open(taffydir + '2017new_xy.reg')
    
    north_list = []
    south_list = []
    snuc_list = []
    hii_list = []
    nw_list = []
    b1_list = []
    b2_list = []
    b3_list = []
    b4_list = []

    for line in region_file.readlines()[3:]:

        line_reg = line.split('(')[1].split(')')[1].split(' ')[-1].rstrip()
        region_list = np.array(line.split('(')[1].split(')')[0].split(','))

        if line_reg == 'north':
            north_pn = get_pn(region_list.astype(np.float64))
            north_mask = bpt.getregionmask(north_pn, (58,58), "North galaxy region.")

        elif line_reg == 'south':
            south_pn = get_pn(region_list.astype(np.float64))
            south_mask = bpt.getregionmask(south_pn, (58,58), "South galaxy region.")

        elif line_reg == 'SNucleus':
            snuc_pn = get_pn(region_list.astype(np.float64))
            snuc_mask = bpt.getregionmask(snuc_pn, (58,58), "South nucleus region.")

        elif line_reg == 'HII':
            hii_pn = get_pn(region_list.astype(np.float64))
            hii_mask = bpt.getregionmask(hii_pn, (58,58), "HII region.")

        elif line_reg == 'NW':
            nw_pn = get_pn(region_list.astype(np.float64))
            nw_mask = bpt.getregionmask(nw_pn, (58,58), "NW region.")

        elif line_reg == 'B1':
            b1_pn = get_pn(region_list.astype(np.float64))
            b1_mask = bpt.getregionmask(b1_pn, (58,58), "B1 region.")

        elif line_reg == 'B2':
            b2_pn = get_pn(region_list.astype(np.float64))
            b2_mask = bpt.getregionmask(b2_pn, (58,58), "B2 region.")

        elif line_reg == 'B3':
            b2_pn = get_pn(region_list.astype(np.float64))
            b2_mask = bpt.getregionmask(b2_pn, (58,58), "B2 region.")

        elif line_reg == 'B4':
            b2_pn = get_pn(region_list.astype(np.float64))
            b2_mask = bpt.getregionmask(b2_pn, (58,58), "B2 region.")

    region_file.close()

    # ---------------------------------------------------------------------------------- #
    # I realized here that I didn't need to apply masks and do the integration myself because 
    # lzifu had alraedy done it.
    # I was trying to see if I could get the same result by applying hte mask and summing 
    # but I don't know why its not the same. I'm just going to use the lzifu results for now.
    #halpha_north_flux_comp1 = np.nansum(ma.array(halpha_withcut_comp1, mask=north_mask))
    #halpha_north_flux_comp2 = np.nansum(ma.array(halpha_withcut_comp2, mask=north_mask))
    #print halpha_north_flux_comp1 + halpha_north_flux_comp2
    #sys.exit(0)

    # these are values taken from the maps that lzifu generates
    # ----------------- total ----------------- #
    oiii_north_flux = 21843.23
    oiii_south_flux = 21084.03
    oiii_snuc_flux = 3737.7152
    oiii_nw_flux = 6572.4901
    oiii_hii_flux = 6708.9706
    oiii_b1_flux = 264.3138
    oiii_b2_flux = 1067.3249
    oiii_b3_flux = 1609.0673
    oiii_b4_flux = 310.27

    hbeta_north_flux = 29221.52
    hbeta_south_flux = 44523.56
    hbeta_snuc_flux = 8596.63
    hbeta_nw_flux = 7010.69
    hbeta_hii_flux = 4846.58
    hbeta_b1_flux = 804.48
    hbeta_b2_flux = 1239.48
    hbeta_b3_flux = 1976.17
    hbeta_b4_flux = 688.71

    halpha_north_flux = 163535.51
    halpha_south_flux = 165746.79
    halpha_snuc_flux = 26159.91
    halpha_nw_flux = 33577.65
    halpha_hii_flux = 29119.93
    halpha_b1_flux = 3995.66
    halpha_b2_flux = 8571.21
    halpha_b3_flux = 9140.25
    halpha_b4_flux = 1688.16

    nii_north_flux = 75755.64
    nii_south_flux = 96363.88
    nii_snuc_flux = 17963.56
    nii_nw_flux = 12702.65
    nii_hii_flux = 12735.77
    nii_b1_flux = 2487.44
    nii_b2_flux = 5436.90
    nii_b3_flux = 4983.16
    nii_b4_flux = 995.5

    oi_north_flux = 13947.77
    oi_south_flux = 12951.69
    oi_snuc_flux = 3489.65
    oi_nw_flux = 3782.80
    oi_hii_flux = 3068.61
    oi_b1_flux = 461.39
    oi_b2_flux = 1310.94
    oi_b3_flux = 1262.94
    oi_b4_flux = 517.61

    sii6716_north_flux = 45969.97
    sii6716_south_flux = 53335.57
    sii6716_snuc_flux = 11700.21
    sii6716_nw_flux = 9152.13
    sii6716_hii_flux = 9097.75
    sii6716_b1_flux = 1534.0
    sii6716_b2_flux = 3451.71
    sii6716_b3_flux = 3846.88
    sii6716_b4_flux = 994.99

    sii6731_north_flux = 29830.19
    sii6731_south_flux = 52422.83
    sii6731_snuc_flux = 10822.02
    sii6731_nw_flux = 6509.84
    sii6731_hii_flux = 5330.52
    sii6731_b1_flux = 2957.79
    sii6731_b2_flux = 3315.76
    sii6731_b3_flux = 3557.57
    sii6731_b4_flux = 817.93

    sii_north_flux = sii6716_north_flux + sii6731_north_flux
    sii_south_flux = sii6716_south_flux + sii6731_south_flux
    sii_snuc_flux = sii6716_snuc_flux + sii6731_snuc_flux
    sii_nw_flux = sii6716_nw_flux + sii6731_nw_flux
    sii_hii_flux = sii6716_hii_flux + sii6731_hii_flux
    sii_b1_flux = sii6716_b1_flux + sii6731_b1_flux
    sii_b2_flux = sii6716_b2_flux + sii6731_b2_flux
    sii_b3_flux = sii6716_b3_flux + sii6731_b3_flux
    sii_b4_flux = sii6716_b4_flux + sii6731_b4_flux

    # ----------------- components ----------------- #
    oiii_north_flux_comp1 = 8006.58
    oiii_north_flux_comp2 = 11692.3
    oiii_south_flux_comp1 = 9801.9
    oiii_south_flux_comp2 = 9924.80
    oiii_snuc_flux_comp1 = 1519.94
    oiii_snuc_flux_comp2 = 2231.78
    oiii_nw_flux_comp1 = 1539.23
    oiii_nw_flux_comp2 = 4384.11
    oiii_hii_flux_comp1 = 2439.05
    oiii_hii_flux_comp2 = 3428.53
    oiii_b1_flux_comp1 = 292.48
    oiii_b1_flux_comp2 = 32.09
    oiii_b2_flux_comp1 = 781.52
    oiii_b2_flux_comp2 = 219.29
    oiii_b3_flux_comp1 = 590.68
    oiii_b3_flux_comp2 = 706.58
    oiii_b4_flux_comp1 = 64.0
    oiii_b4_flux_comp2 = 74.71

    hbeta_north_flux_comp1 = 13254.99
    hbeta_north_flux_comp2 = 13489.56
    hbeta_south_flux_comp1 = 27616.5
    hbeta_south_flux_comp2 = 14154.16
    hbeta_snuc_flux_comp1 = 4366.42
    hbeta_snuc_flux_comp2 = 4223.38
    hbeta_nw_flux_comp1 = 1971.03
    hbeta_nw_flux_comp2 = 4313.1
    hbeta_hii_flux_comp1 = 1287.81
    hbeta_hii_flux_comp2 = 3032.56
    hbeta_b1_flux_comp1 = 620.46
    hbeta_b1_flux_comp2 = 154.70
    hbeta_b2_flux_comp1 = 449.53
    hbeta_b2_flux_comp2 = 596.98
    hbeta_b3_flux_comp1 = 971.38
    hbeta_b3_flux_comp2 = 645.29
    hbeta_b4_flux_comp1 = 209.81
    hbeta_b4_flux_comp2 = 207.29

    halpha_north_flux_comp1 = 60450.61
    halpha_north_flux_comp2 = 89235.34
    halpha_south_flux_comp1 = 102365.44
    halpha_south_flux_comp2 = 53491.50
    halpha_snuc_flux_comp1 = 14747.68
    halpha_snuc_flux_comp2 = 11392.43
    halpha_nw_flux_comp1 = 11302.96
    halpha_nw_flux_comp2 = 18922.63
    halpha_hii_flux_comp1 = 8318.10
    halpha_hii_flux_comp2 = 16875.69
    halpha_b1_flux_comp1 = 2999.96
    halpha_b1_flux_comp2 = 652.94
    halpha_b2_flux_comp1 = 4449.25
    halpha_b2_flux_comp2 = 2671.71
    halpha_b3_flux_comp1 = 4013.86
    halpha_b3_flux_comp2 = 3783.28
    halpha_b4_flux_comp1 = 437.60
    halpha_b4_flux_comp2 = 617.03

    nii_north_flux_comp1 = 27539.57
    nii_north_flux_comp2 = 42102.21
    nii_south_flux_comp1 = 58989.62
    nii_south_flux_comp2 = 32498.1
    nii_snuc_flux_comp1 = 10201.75
    nii_snuc_flux_comp2 = 7738.58
    nii_nw_flux_comp1 = 4366.56
    nii_nw_flux_comp2 = 7032.29
    nii_hii_flux_comp1 = 2615.32
    nii_hii_flux_comp2 = 8401.42
    nii_b1_flux_comp1 = 1391.08
    nii_b1_flux_comp2 = 887.48
    nii_b2_flux_comp1 = 3305.13
    nii_b2_flux_comp2 = 2175.46
    nii_b3_flux_comp1 = 2289.91
    nii_b3_flux_comp2 = 2113.09
    nii_b4_flux_comp1 = 171.47
    nii_b4_flux_comp2 = 440.56

    oi_north_flux_comp1 = 7104.49
    oi_north_flux_comp2 = 5665.37
    oi_south_flux_comp1 = 6674.42
    oi_south_flux_comp2 = 5781.27
    oi_snuc_flux_comp1 = 1647.43
    oi_snuc_flux_comp2 = 1857.39
    oi_nw_flux_comp1 = 2644.68
    oi_nw_flux_comp2 = 980.03
    oi_hii_flux_comp1 = 378.98
    oi_hii_flux_comp2 = 2082.59
    oi_b1_flux_comp1 = 204.02
    oi_b1_flux_comp2 = 161.65
    oi_b2_flux_comp1 = 619.97
    oi_b2_flux_comp2 = 531.53
    oi_b3_flux_comp1 = 307.4
    oi_b3_flux_comp2 = 709.12
    oi_b4_flux_comp1 = 65.95
    oi_b4_flux_comp2 = 363.14

    sii6716_north_flux_comp1 = 18096.38
    sii6716_north_flux_comp2 = 23708.17
    sii6716_south_flux_comp1 = 27330.28
    sii6716_south_flux_comp2 = 23355.61
    sii6716_snuc_flux_comp1 = 5536.94
    sii6716_snuc_flux_comp2 = 6147.53
    sii6716_nw_flux_comp1 = 3337.32
    sii6716_nw_flux_comp2 = 4849.76
    sii6716_hii_flux_comp1 = 2293.01
    sii6716_hii_flux_comp2 = 5379.43
    sii6716_b1_flux_comp1 = 843.71
    sii6716_b1_flux_comp2 = 589.41
    sii6716_b2_flux_comp1 = 1931.12
    sii6716_b2_flux_comp2 = 1157.39
    sii6716_b3_flux_comp1 = 1510.16
    sii6716_b3_flux_comp2 = 2008.59
    sii6716_b4_flux_comp1 = 227.17
    sii6716_b4_flux_comp2 = 489.64

    sii6731_north_flux_comp1 = 14453.93
    sii6731_north_flux_comp2 = 13554.63
    sii6731_south_flux_comp1 = 23469.77
    sii6731_south_flux_comp2 = 26016.76
    sii6731_snuc_flux_comp1 = 3811.16
    sii6731_snuc_flux_comp2 = 7080.46
    sii6731_nw_flux_comp1 = 3887.60
    sii6731_nw_flux_comp2 = 2257.74
    sii6731_hii_flux_comp1 = 1076.64
    sii6731_hii_flux_comp2 = 3509.79
    sii6731_b1_flux_comp1 = 1107.85
    sii6731_b1_flux_comp2 = 2325.75
    sii6731_b2_flux_comp1 = 1109.61
    sii6731_b2_flux_comp2 = 1558.74
    sii6731_b3_flux_comp1 = 1862.15
    sii6731_b3_flux_comp2 = 1176.41
    sii6731_b4_flux_comp1 = 282.76
    sii6731_b4_flux_comp2 = 283.4

    sii_north_flux_comp1 = sii6716_north_flux_comp1 + sii6731_north_flux_comp1
    sii_south_flux_comp1 = sii6716_south_flux_comp1 + sii6731_south_flux_comp1
    sii_snuc_flux_comp1 = sii6716_snuc_flux_comp1 + sii6731_snuc_flux_comp1
    sii_nw_flux_comp1 = sii6716_nw_flux_comp1 + sii6731_nw_flux_comp1
    sii_hii_flux_comp1 = sii6716_hii_flux_comp1 + sii6731_hii_flux_comp1
    sii_b1_flux_comp1 = sii6716_b1_flux_comp1 + sii6731_b1_flux_comp1
    sii_b2_flux_comp1 = sii6716_b2_flux_comp1 + sii6731_b2_flux_comp1
    sii_b3_flux_comp1 = sii6716_b3_flux_comp1 + sii6731_b3_flux_comp1
    sii_b4_flux_comp1 = sii6716_b4_flux_comp1 + sii6731_b4_flux_comp1

    sii_north_flux_comp2 = sii6716_north_flux_comp2 + sii6731_north_flux_comp2
    sii_south_flux_comp2 = sii6716_south_flux_comp2 + sii6731_south_flux_comp2
    sii_snuc_flux_comp2 = sii6716_snuc_flux_comp2 + sii6731_snuc_flux_comp2
    sii_nw_flux_comp2 = sii6716_nw_flux_comp2 + sii6731_nw_flux_comp2
    sii_hii_flux_comp2 = sii6716_hii_flux_comp2 + sii6731_hii_flux_comp2
    sii_b1_flux_comp2 = sii6716_b1_flux_comp2 + sii6731_b1_flux_comp2
    sii_b2_flux_comp2 = sii6716_b2_flux_comp2 + sii6731_b2_flux_comp2
    sii_b3_flux_comp2 = sii6716_b3_flux_comp2 + sii6731_b3_flux_comp2
    sii_b4_flux_comp2 = sii6716_b4_flux_comp2 + sii6731_b4_flux_comp2

    # make figure
    # first define ratios
    # north
    oiii_hbeta_north = np.log10(oiii_north_flux / hbeta_north_flux)
    nii_halpha_north = np.log10(nii_north_flux / halpha_north_flux)
    oi_halpha_north = np.log10(oi_north_flux / halpha_north_flux)
    sii_halpha_north = np.log10(sii_north_flux / halpha_north_flux)

    oiii_hbeta_north_comp1 = np.log10(oiii_north_flux_comp1 / hbeta_north_flux_comp1)
    nii_halpha_north_comp1 = np.log10(nii_north_flux_comp1 / halpha_north_flux_comp1)
    oi_halpha_north_comp1 = np.log10(oi_north_flux_comp1 / halpha_north_flux_comp1)
    sii_halpha_north_comp1 = np.log10(sii_north_flux_comp1 / halpha_north_flux_comp1)

    oiii_hbeta_north_comp2 = np.log10(oiii_north_flux_comp2 / hbeta_north_flux_comp2)
    nii_halpha_north_comp2 = np.log10(nii_north_flux_comp2 / halpha_north_flux_comp2)
    oi_halpha_north_comp2 = np.log10(oi_north_flux_comp2 / halpha_north_flux_comp2)
    sii_halpha_north_comp2 = np.log10(sii_north_flux_comp2 / halpha_north_flux_comp2)

    # south
    oiii_hbeta_south = np.log10(oiii_south_flux / hbeta_south_flux)
    nii_halpha_south = np.log10(nii_south_flux / halpha_south_flux)
    oi_halpha_south = np.log10(oi_south_flux / halpha_south_flux)
    sii_halpha_south = np.log10(sii_south_flux / halpha_south_flux)

    oiii_hbeta_south_comp1 = np.log10(oiii_south_flux_comp1 / hbeta_south_flux_comp1)
    nii_halpha_south_comp1 = np.log10(nii_south_flux_comp1 / halpha_south_flux_comp1)
    oi_halpha_south_comp1 = np.log10(oi_south_flux_comp1 / halpha_south_flux_comp1)
    sii_halpha_south_comp1 = np.log10(sii_south_flux_comp1 / halpha_south_flux_comp1)

    oiii_hbeta_south_comp2 = np.log10(oiii_south_flux_comp2 / hbeta_south_flux_comp2)
    nii_halpha_south_comp2 = np.log10(nii_south_flux_comp2 / halpha_south_flux_comp2)
    oi_halpha_south_comp2 = np.log10(oi_south_flux_comp2 / halpha_south_flux_comp2)
    sii_halpha_south_comp2 = np.log10(sii_south_flux_comp2 / halpha_south_flux_comp2)

    # snuc
    oiii_hbeta_snuc = np.log10(oiii_snuc_flux / hbeta_snuc_flux)
    nii_halpha_snuc = np.log10(nii_snuc_flux / halpha_snuc_flux)
    oi_halpha_snuc = np.log10(oi_snuc_flux / halpha_snuc_flux)
    sii_halpha_snuc = np.log10(sii_snuc_flux / halpha_snuc_flux)

    oiii_hbeta_snuc_comp1 = np.log10(oiii_snuc_flux_comp1 / hbeta_snuc_flux_comp1)
    nii_halpha_snuc_comp1 = np.log10(nii_snuc_flux_comp1 / halpha_snuc_flux_comp1)
    oi_halpha_snuc_comp1 = np.log10(oi_snuc_flux_comp1 / halpha_snuc_flux_comp1)
    sii_halpha_snuc_comp1 = np.log10(sii_snuc_flux_comp1 / halpha_snuc_flux_comp1)

    oiii_hbeta_snuc_comp2 = np.log10(oiii_snuc_flux_comp2 / hbeta_snuc_flux_comp2)
    nii_halpha_snuc_comp2 = np.log10(nii_snuc_flux_comp2 / halpha_snuc_flux_comp2)
    oi_halpha_snuc_comp2 = np.log10(oi_snuc_flux_comp2 / halpha_snuc_flux_comp2)
    sii_halpha_snuc_comp2 = np.log10(sii_snuc_flux_comp2 / halpha_snuc_flux_comp2)

    # nw
    oiii_hbeta_nw = np.log10(oiii_nw_flux / hbeta_nw_flux)
    nii_halpha_nw = np.log10(nii_nw_flux / halpha_nw_flux)
    oi_halpha_nw = np.log10(oi_nw_flux / halpha_nw_flux)
    sii_halpha_nw = np.log10(sii_nw_flux / halpha_nw_flux)

    oiii_hbeta_nw_comp1 = np.log10(oiii_nw_flux_comp1 / hbeta_nw_flux_comp1)
    nii_halpha_nw_comp1 = np.log10(nii_nw_flux_comp1 / halpha_nw_flux_comp1)
    oi_halpha_nw_comp1 = np.log10(oi_nw_flux_comp1 / halpha_nw_flux_comp1)
    sii_halpha_nw_comp1 = np.log10(sii_nw_flux_comp1 / halpha_nw_flux_comp1)

    oiii_hbeta_nw_comp2 = np.log10(oiii_nw_flux_comp2 / hbeta_nw_flux_comp2)
    nii_halpha_nw_comp2 = np.log10(nii_nw_flux_comp2 / halpha_nw_flux_comp2)
    oi_halpha_nw_comp2 = np.log10(oi_nw_flux_comp2 / halpha_nw_flux_comp2)
    sii_halpha_nw_comp2 = np.log10(sii_nw_flux_comp2 / halpha_nw_flux_comp2)

    # hii
    oiii_hbeta_hii = np.log10(oiii_hii_flux / hbeta_hii_flux)
    nii_halpha_hii = np.log10(nii_hii_flux / halpha_hii_flux)
    oi_halpha_hii = np.log10(oi_hii_flux / halpha_hii_flux)
    sii_halpha_hii = np.log10(sii_hii_flux / halpha_hii_flux)

    oiii_hbeta_hii_comp1 = np.log10(oiii_hii_flux_comp1 / hbeta_hii_flux_comp1)
    nii_halpha_hii_comp1 = np.log10(nii_hii_flux_comp1 / halpha_hii_flux_comp1)
    oi_halpha_hii_comp1 = np.log10(oi_hii_flux_comp1 / halpha_hii_flux_comp1)
    sii_halpha_hii_comp1 = np.log10(sii_hii_flux_comp1 / halpha_hii_flux_comp1)

    oiii_hbeta_hii_comp2 = np.log10(oiii_hii_flux_comp2 / hbeta_hii_flux_comp2)
    nii_halpha_hii_comp2 = np.log10(nii_hii_flux_comp2 / halpha_hii_flux_comp2)
    oi_halpha_hii_comp2 = np.log10(oi_hii_flux_comp2 / halpha_hii_flux_comp2)
    sii_halpha_hii_comp2 = np.log10(sii_hii_flux_comp2 / halpha_hii_flux_comp2)

    # b1
    oiii_hbeta_b1 = np.log10(oiii_b1_flux / hbeta_b1_flux)
    nii_halpha_b1 = np.log10(nii_b1_flux / halpha_b1_flux)
    oi_halpha_b1 = np.log10(oi_b1_flux / halpha_b1_flux)
    sii_halpha_b1 = np.log10(sii_b1_flux / halpha_b1_flux)

    oiii_hbeta_b1_comp1 = np.log10(oiii_b1_flux_comp1 / hbeta_b1_flux_comp1)
    nii_halpha_b1_comp1 = np.log10(nii_b1_flux_comp1 / halpha_b1_flux_comp1)
    oi_halpha_b1_comp1 = np.log10(oi_b1_flux_comp1 / halpha_b1_flux_comp1)
    sii_halpha_b1_comp1 = np.log10(sii_b1_flux_comp1 / halpha_b1_flux_comp1)

    oiii_hbeta_b1_comp2 = np.log10(oiii_b1_flux_comp2 / hbeta_b1_flux_comp2)
    nii_halpha_b1_comp2 = np.log10(nii_b1_flux_comp2 / halpha_b1_flux_comp2)
    oi_halpha_b1_comp2 = np.log10(oi_b1_flux_comp2 / halpha_b1_flux_comp2)
    sii_halpha_b1_comp2 = np.log10(sii_b1_flux_comp2 / halpha_b1_flux_comp2)

    # b2
    oiii_hbeta_b2 = np.log10(oiii_b2_flux / hbeta_b2_flux)
    nii_halpha_b2 = np.log10(nii_b2_flux / halpha_b2_flux)
    oi_halpha_b2 = np.log10(oi_b2_flux / halpha_b2_flux)
    sii_halpha_b2 = np.log10(sii_b2_flux / halpha_b2_flux)

    oiii_hbeta_b2_comp1 = np.log10(oiii_b2_flux_comp1 / hbeta_b2_flux_comp1)
    nii_halpha_b2_comp1 = np.log10(nii_b2_flux_comp1 / halpha_b2_flux_comp1)
    oi_halpha_b2_comp1 = np.log10(oi_b2_flux_comp1 / halpha_b2_flux_comp1)
    sii_halpha_b2_comp1 = np.log10(sii_b2_flux_comp1 / halpha_b2_flux_comp1)

    oiii_hbeta_b2_comp2 = np.log10(oiii_b2_flux_comp2 / hbeta_b2_flux_comp2)
    nii_halpha_b2_comp2 = np.log10(nii_b2_flux_comp2 / halpha_b2_flux_comp2)
    oi_halpha_b2_comp2 = np.log10(oi_b2_flux_comp2 / halpha_b2_flux_comp2)
    sii_halpha_b2_comp2 = np.log10(sii_b2_flux_comp2 / halpha_b2_flux_comp2)

    # b3
    oiii_hbeta_b3 = np.log10(oiii_b3_flux / hbeta_b3_flux)
    nii_halpha_b3 = np.log10(nii_b3_flux / halpha_b3_flux)
    oi_halpha_b3 = np.log10(oi_b3_flux / halpha_b3_flux)
    sii_halpha_b3 = np.log10(sii_b3_flux / halpha_b3_flux)

    oiii_hbeta_b3_comp1 = np.log10(oiii_b3_flux_comp1 / hbeta_b3_flux_comp1)
    nii_halpha_b3_comp1 = np.log10(nii_b3_flux_comp1 / halpha_b3_flux_comp1)
    oi_halpha_b3_comp1 = np.log10(oi_b3_flux_comp1 / halpha_b3_flux_comp1)
    sii_halpha_b3_comp1 = np.log10(sii_b3_flux_comp1 / halpha_b3_flux_comp1)

    oiii_hbeta_b3_comp2 = np.log10(oiii_b3_flux_comp2 / hbeta_b3_flux_comp2)
    nii_halpha_b3_comp2 = np.log10(nii_b3_flux_comp2 / halpha_b3_flux_comp2)
    oi_halpha_b3_comp2 = np.log10(oi_b3_flux_comp2 / halpha_b3_flux_comp2)
    sii_halpha_b3_comp2 = np.log10(sii_b3_flux_comp2 / halpha_b3_flux_comp2)

    # b1
    oiii_hbeta_b4 = np.log10(oiii_b4_flux / hbeta_b4_flux)
    nii_halpha_b4 = np.log10(nii_b4_flux / halpha_b4_flux)
    oi_halpha_b4 = np.log10(oi_b4_flux / halpha_b4_flux)
    sii_halpha_b4 = np.log10(sii_b4_flux / halpha_b4_flux)

    oiii_hbeta_b4_comp1 = np.log10(oiii_b4_flux_comp1 / hbeta_b4_flux_comp1)
    nii_halpha_b4_comp1 = np.log10(nii_b4_flux_comp1 / halpha_b4_flux_comp1)
    oi_halpha_b4_comp1 = np.log10(oi_b4_flux_comp1 / halpha_b4_flux_comp1)
    sii_halpha_b4_comp1 = np.log10(sii_b4_flux_comp1 / halpha_b4_flux_comp1)

    oiii_hbeta_b4_comp2 = np.log10(oiii_b4_flux_comp2 / hbeta_b4_flux_comp2)
    nii_halpha_b4_comp2 = np.log10(nii_b4_flux_comp2 / halpha_b4_flux_comp2)
    oi_halpha_b4_comp2 = np.log10(oi_b4_flux_comp2 / halpha_b4_flux_comp2)
    sii_halpha_b4_comp2 = np.log10(sii_b4_flux_comp2 / halpha_b4_flux_comp2)

    # total
    plotbpt('nii', 'total', nii_halpha_north, oiii_hbeta_north, nii_halpha_south, oiii_hbeta_south, nii_halpha_snuc, oiii_hbeta_snuc,\
    nii_halpha_nw, oiii_hbeta_nw, nii_halpha_hii, oiii_hbeta_hii, nii_halpha_b1, oiii_hbeta_b1, nii_halpha_b2, oiii_hbeta_b2,\
    nii_halpha_b3, oiii_hbeta_b3, nii_halpha_b4, oiii_hbeta_b4)

    plotbpt('oi', 'total', oi_halpha_north, oiii_hbeta_north, oi_halpha_south, oiii_hbeta_south, oi_halpha_snuc, oiii_hbeta_snuc,\
    oi_halpha_nw, oiii_hbeta_nw, oi_halpha_hii, oiii_hbeta_hii, oi_halpha_b1, oiii_hbeta_b1, oi_halpha_b2, oiii_hbeta_b2,\
    oi_halpha_b3, oiii_hbeta_b3, oi_halpha_b4, oiii_hbeta_b4)

    plotbpt('sii', 'total', sii_halpha_north, oiii_hbeta_north, sii_halpha_south, oiii_hbeta_south, sii_halpha_snuc, oiii_hbeta_snuc,\
    sii_halpha_nw, oiii_hbeta_nw, sii_halpha_hii, oiii_hbeta_hii, sii_halpha_b1, oiii_hbeta_b1, sii_halpha_b2, oiii_hbeta_b2,\
    sii_halpha_b3, oiii_hbeta_b3, sii_halpha_b4, oiii_hbeta_b4)

    # comp 1
    plotbpt('nii', '1', nii_halpha_north_comp1, oiii_hbeta_north_comp1, nii_halpha_south_comp1, oiii_hbeta_south_comp1, nii_halpha_snuc_comp1, oiii_hbeta_snuc_comp1,\
    nii_halpha_nw_comp1, oiii_hbeta_nw_comp1, nii_halpha_hii_comp1, oiii_hbeta_hii_comp1, nii_halpha_b1_comp1, oiii_hbeta_b1_comp1, nii_halpha_b2_comp1, oiii_hbeta_b2_comp1,\
    nii_halpha_b3_comp1, oiii_hbeta_b3_comp1, nii_halpha_b4_comp1, oiii_hbeta_b4_comp1)

    plotbpt('oi', '1', oi_halpha_north_comp1, oiii_hbeta_north_comp1, oi_halpha_south_comp1, oiii_hbeta_south_comp1, oi_halpha_snuc_comp1, oiii_hbeta_snuc_comp1,\
    oi_halpha_nw_comp1, oiii_hbeta_nw_comp1, oi_halpha_hii_comp1, oiii_hbeta_hii_comp1, oi_halpha_b1_comp1, oiii_hbeta_b1_comp1, oi_halpha_b2_comp1, oiii_hbeta_b2_comp1,\
    oi_halpha_b3_comp1, oiii_hbeta_b3_comp1, oi_halpha_b4_comp1, oiii_hbeta_b4_comp1)

    plotbpt('sii', '1', sii_halpha_north_comp1, oiii_hbeta_north_comp1, sii_halpha_south_comp1, oiii_hbeta_south_comp1, sii_halpha_snuc_comp1, oiii_hbeta_snuc_comp1,\
    sii_halpha_nw_comp1, oiii_hbeta_nw_comp1, sii_halpha_hii_comp1, oiii_hbeta_hii_comp1, sii_halpha_b1_comp1, oiii_hbeta_b1_comp1, sii_halpha_b2_comp1, oiii_hbeta_b2_comp1,\
    sii_halpha_b3_comp1, oiii_hbeta_b3_comp1, sii_halpha_b4_comp1, oiii_hbeta_b4_comp1)

    # comp2
    plotbpt('nii', '2', nii_halpha_north_comp2, oiii_hbeta_north_comp2, nii_halpha_south_comp2, oiii_hbeta_south_comp2, nii_halpha_snuc_comp2, oiii_hbeta_snuc_comp2,\
    nii_halpha_nw_comp2, oiii_hbeta_nw_comp2, nii_halpha_hii_comp2, oiii_hbeta_hii_comp2, nii_halpha_b1_comp2, oiii_hbeta_b1_comp2, nii_halpha_b2_comp2, oiii_hbeta_b2_comp2,\
    nii_halpha_b3_comp2, oiii_hbeta_b3_comp2, nii_halpha_b4_comp2, oiii_hbeta_b4_comp2)

    plotbpt('oi', '2', oi_halpha_north_comp2, oiii_hbeta_north_comp2, oi_halpha_south_comp2, oiii_hbeta_south_comp1, oi_halpha_snuc_comp2, oiii_hbeta_snuc_comp2,\
    oi_halpha_nw_comp2, oiii_hbeta_nw_comp2, oi_halpha_hii_comp2, oiii_hbeta_hii_comp2, oi_halpha_b1_comp2, oiii_hbeta_b1_comp1, oi_halpha_b2_comp2, oiii_hbeta_b2_comp2,\
    oi_halpha_b3_comp2, oiii_hbeta_b3_comp2, oi_halpha_b4_comp2, oiii_hbeta_b4_comp2)

    plotbpt('sii', '2', sii_halpha_north_comp2, oiii_hbeta_north_comp2, sii_halpha_south_comp2, oiii_hbeta_south_comp2, sii_halpha_snuc_comp2, oiii_hbeta_snuc_comp2,\
    sii_halpha_nw_comp2, oiii_hbeta_nw_comp2, sii_halpha_hii_comp2, oiii_hbeta_hii_comp2, sii_halpha_b1_comp2, oiii_hbeta_b1_comp2, sii_halpha_b2_comp2, oiii_hbeta_b2_comp2,\
    sii_halpha_b3_comp2, oiii_hbeta_b3_comp2, sii_halpha_b4_comp2, oiii_hbeta_b4_comp2)

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)