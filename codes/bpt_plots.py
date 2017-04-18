from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg

import os
import sys
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')  # Does not have a trailing slash at the end
stacking_analysis_dir = home + "/Desktop/FIGS/stacking-analysis-pears/"
ipac_taffy_figdir = home + "/Desktop/ipac/taffy/figures/"
taffy_products = '/Volumes/Bhavins_backup/phil/TAFFY/products/'

sys.path.append(stacking_analysis_dir + 'codes/')
import fast_chi2_jackknife as fcj

def check_3darr_equality(arr, arrname):
    """
    This function checks if the first and second dimensions of a 3D array
    have the same elements.
    I wrote this function because I noticed that the line maps have shapes
    (2,58,58) but I had expected only (1, 58, 58) i.e. just (58,58). So this
    function checks if the zeroth and the first row in the first dimension
    of the array have all elements equal.

    I also noticed that many elements are NaN. So this function first checks
    if the elements are finite and then checks for equality. If the elements 
    are equal then it proceeds to the next pair otherwise if they are not
    equal then it prints the absolute value of the difference between them.
    """

    print "Checking", arrname, "now..."

    eq = np.allclose(arr[0], arr[1], equal_nan=True)
    print "Using Numpy's allclose which checks elementwise equality," +\
    " the statement that the" + " two rows in the zeorth dimension are equal is --", eq

    print "Now checking elementwise in a for loop."
    print "If it doesn't print anything then the answer is True."

    for i in range(len(arr[0][0])):
        for j in range(len(arr[0][1])):
            if np.isfinite(arr[0][i,j]) and np.isfinite(arr[1][i,j]):
                if arr[0][i,j] == arr[1][i,j]:
                    continue
                else:
                    diff = abs(arr[0][i,j] - arr[1][i,j])
                    if diff > 1e-4:
                        print diff

    return None

def getregionmask(pn, final_arr_shape, regionname):
    """
    This function will return a Boolean mask defining the region of interest.

    Arguments:
    pn: Object of Polygon class. 
    This is the polygon defining the region of interest.

    final_arr_shape: Numpy array shape.
    This is two integers separated by a comma and enclosed in parentheses.
    """

    regionmask = np.zeros(final_arr_shape, dtype=np.int)

    bbox = pn.boundingBox()
    bbox_xmin = int(bbox[0])
    bbox_xmax = int(bbox[1])
    bbox_ymin = int(bbox[2])
    bbox_ymax = int(bbox[3])

    if bbox_xmin < 0:
        bbox_xmin = 0
    if bbox_xmax > 58:
        bbox_xmax = 58
    if bbox_ymin < 0:
        bbox_ymin = 0
    if bbox_ymax > 58:
        bbox_ymax = 58

    # loop over all possible points inside the bounding box
    for i in range(bbox_xmin, bbox_xmax):
        for j in range(bbox_ymin, bbox_ymax):

            # convert from ds9 coords to array coords
            # this just needs to switch between x and y without the -1???
            arr_x = int(j)
            arr_y = int(i)

            regionmask[arr_x,arr_y] = pn.isInside(i,j)

    # loop over all vertices because otherwise Polygon.isInside() 
    # may or may not return true when the argument is a vertex point
    for k in range(pn.nPoints()):

        arr_x = int(pn[0][k][1])
        arr_y = int(pn[0][k][0])

        if arr_x < 0:
            arr_x = 0
        if arr_x > 57:
            arr_x = 57
        if arr_y < 0:
            arr_y = 0
        if arr_y > 57:
            arr_y = 57

        regionmask[arr_x, arr_y] = True

    reg_idx_zero = np.where(regionmask == 0.0)
    reg_idx_one = np.where(regionmask == 1.0)
    regionmask[reg_idx_zero] = 1.0
    regionmask[reg_idx_one] = 0.0

    print len(reg_idx_one[0].flatten()), "points in the", regionname

    return regionmask

def zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, orig_shape):

    zipped_ha_ind = zip(ha_valid_ind[0], ha_valid_ind[1])
    zipped_hb_ind = zip(hb_valid_ind[0], hb_valid_ind[1])
    zipped_oiii_ind = zip(oiii_valid_ind[0], oiii_valid_ind[1])
    zipped_nii_ind = zip(nii_valid_ind[0], nii_valid_ind[1])
    zipped_oi_ind = zip(oi_valid_ind[0], oi_valid_ind[1])
    zipped_sii_ind = zip(sii_valid_ind[0], sii_valid_ind[1])

    nii_halpha_withcut = np.zeros(orig_shape)
    oiii_hbeta_withcut = np.zeros(orig_shape)
    oi_halpha_withcut = np.zeros(orig_shape)
    sii_halpha_withcut = np.zeros(orig_shape)

    halpha_withcut = np.zeros(orig_shape)
    hbeta_withcut = np.zeros(orig_shape)
    oiii5007_withcut = np.zeros(orig_shape)
    oi6300_withcut = np.zeros(orig_shape)
    nii6583_withcut = np.zeros(orig_shape)
    sii_withcut = np.zeros(orig_shape)

    for i in range(len(zipped_ha_ind)):

        if zipped_ha_ind[i] in zipped_oi_ind:
            if zipped_ha_ind[i] in zipped_nii_ind:
                if zipped_ha_ind[i] in zipped_hb_ind:
                    if zipped_ha_ind[i] in zipped_oiii_ind:
                        if zipped_ha_ind[i] in zipped_sii_ind:
                            nii_halpha_withcut[zipped_ha_ind[i]] =\
                             np.log10(nii6583[0][zipped_ha_ind[i]] / halpha[0][zipped_ha_ind[i]])
                            oi_halpha_withcut[zipped_ha_ind[i]] =\
                             np.log10(oi6300[0][zipped_ha_ind[i]] / halpha[0][zipped_ha_ind[i]])
                            oiii_hbeta_withcut[zipped_ha_ind[i]] =\
                             np.log10(oiii5007[0][zipped_ha_ind[i]] / hbeta[0][zipped_ha_ind[i]])
                            sii_halpha_withcut[zipped_ha_ind[i]] =\
                             np.log10(sii[zipped_ha_ind[i]] / halpha[0][zipped_ha_ind[i]])

                            halpha_withcut[zipped_ha_ind[i]] = halpha[0][zipped_ha_ind[i]]
                            hbeta_withcut[zipped_ha_ind[i]] = hbeta[0][zipped_ha_ind[i]]
                            oiii5007_withcut[zipped_ha_ind[i]] = oiii5007[0][zipped_ha_ind[i]]
                            oi6300_withcut[zipped_ha_ind[i]] = oi6300[0][zipped_ha_ind[i]]
                            nii6583_withcut[zipped_ha_ind[i]] = nii6583[0][zipped_ha_ind[i]]
                            sii_withcut[zipped_ha_ind[i]] = sii[zipped_ha_ind[i]]

    return nii_halpha_withcut, oiii_hbeta_withcut, oi_halpha_withcut, sii_halpha_withcut,\
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut

def get_arr_withsigcut(sig_cut, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, orig_shape):

    # apply sigma cuts and create arrays for bpt diagram
    ha_valid_ind = np.where((halpha[0]/halpha_err[0]) > sig_cut)
    hb_valid_ind = np.where((hbeta[0]/hbeta_err[0]) > sig_cut)
    oiii_valid_ind = np.where((oiii5007[0]/oiii5007_err[0]) > sig_cut)
    nii_valid_ind = np.where((nii6583[0]/nii6583_err[0]) > sig_cut)
    oi_valid_ind = np.where((oi6300[0]/oi6300_err[0]) > sig_cut)
    sii_valid_ind = np.where((sii/sii_err) > sig_cut)

    nii_halpha_withsigcut, oiii_hbeta_withsigcut, oi_halpha_withsigcut, sii_halpha_withsigcut,\
    halpha_withsigcut, hbeta_withsigcut, oiii5007_withsigcut, oi6300_withsigcut, nii6583_withsigcut, sii_withsigcut =\
    zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, orig_shape)    

    return nii_halpha_withsigcut, oiii_hbeta_withsigcut, oi_halpha_withsigcut, sii_halpha_withsigcut,\
    halpha_withsigcut, hbeta_withsigcut, oiii5007_withsigcut, oi6300_withsigcut, nii6583_withsigcut, sii_withsigcut

def get_arr_basecut(halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, orig_shape):

    ha_valid_ind = np.where(halpha[0] >= 10)
    hb_valid_ind = np.where(hbeta[0] >= 4)
    oiii_valid_ind = np.where(oiii5007[0] >= 1)
    nii_valid_ind = np.where(nii6583[0] >= 8)
    oi_valid_ind = np.where(oi6300[0] >= 0.5)
    sii_valid_ind = np.where(sii >= 6)

    nii_halpha_withbasecut, oiii_hbeta_withbasecut, oi_halpha_withbasecut, sii_halpha_withbasecut,\
    halpha_withbasecut, hbeta_withbasecut, oiii5007_withbasecut, oi6300_withbasecut, nii6583_withbasecut, sii_withbasecut =\
    zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, orig_shape) 

    return nii_halpha_withbasecut, oiii_hbeta_withbasecut, oi_halpha_withbasecut, sii_halpha_withbasecut,\
    halpha_withbasecut, hbeta_withbasecut, oiii5007_withbasecut, oi6300_withbasecut, nii6583_withbasecut, sii_withbasecut

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in lzifu output file
    h = fits.open(taffy_products + 'TAFFY_2_comp.fits')
    hdu_vdisp = fits.open(taffy_products + 'TAFFY_2_comp_VDISP.fits')

    # loop over extensions 
    total_ext = fcj.get_total_extensions(h)

    # assign line arrays and check equality; read func definition
    halpha = h['HALPHA'].data
    hbeta = h['HBETA'].data
    nii6583 = h['NII6583'].data
    oiii5007 = h['OIII5007'].data
    oi6300 = h['OI6300'].data
    oi6364 = h['OI6364'].data
    sii6716 = h['SII6716'].data
    sii6731 = h['SII6731'].data
    vdisp_line1 = hdu_vdisp[0].data[1]
    vdisp_line2 = hdu_vdisp[0].data[2]

    halpha_err = h['HALPHA_ERR'].data
    hbeta_err = h['HBETA_ERR'].data
    nii6583_err = h['NII6583_ERR'].data
    oiii5007_err = h['OIII5007_ERR'].data
    oi6300_err = h['OI6300_ERR'].data
    oi6364_err = h['OI6364_ERR'].data
    sii6716_err = h['SII6716_ERR'].data
    sii6731_err = h['SII6731_ERR'].data

    # add lines which are doublets
    sii = sii6716[0] + sii6731[0]
    sii_err = np.sqrt((sii6716_err[0])**2 + (sii6731_err[0])**2)

    # get arrays with some baseline level cut off
    nii_halpha_withcut, oiii_hbeta_withcut, oi_halpha_withcut, sii_halpha_withcut,\
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut  =\
     get_arr_basecut(halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha[0].shape)

    # get arrays with significance cut applied
    nii_halpha_withcut, oiii_hbeta_withcut, oi_halpha_withcut, sii_halpha_withcut,\
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut =\
     get_arr_withsigcut(3, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, halpha[0].shape)

    # spatial mapping for interesting bpt points
    oi_bpt_array = np.zeros(oiii_hbeta_withcut.shape)
    oi_bpt_points = []

    oi_agn_hii = lambda x: 1.33 + 0.73 / (x + 0.59)
    oi_liner_seyfert = lambda x: 1.30 + 1.18 * x

    for u in range(oiii_hbeta_withcut.shape[0]):

        for v in range(oiii_hbeta_withcut.shape[1]):

            current_y = oiii_hbeta_withcut[u,v]
            current_x = oi_halpha_withcut[u,v]

            # for LINER region
            if (current_y > oi_agn_hii(current_x)) and (current_y < oi_liner_seyfert(current_x)):
                oi_bpt_points.append([current_x, current_y])
                #print current_x, current_y, u, v
                oi_bpt_array[u,v] = 1.0
                # add 1 to both u and v to convert to ds9 y and x respectively

            # for HII region
            # its currently giving all x values that are to the right of the HII line
            # which is not what I want
            # i think this is because the function has an asymptote somewhere to the right of the line drawn
            if (current_y <= oi_agn_hii(current_x)) and (current_y != 0.0):
                print current_x, current_y
                oi_bpt_array[u,v] = 2.0

            # for SEYFERT region
            if (current_y > oi_agn_hii(current_x)) and (current_y > oi_liner_seyfert(current_x)):
                oi_bpt_array[u,v] = 3.0

    plt.imshow(oi_bpt_array, origin='lower')
    plt.colorbar()
    plt.show()

    sys.exit(0)

    # arrays without any cut
    nii_halpha = np.log10(np.divide(nii6583[0], halpha[0]))
    oiii_hbeta = np.log10(np.divide(oiii5007[0], hbeta[0]))
    oi_halpha = np.log10(np.divide(oi6300[0], halpha[0]))

    # define regions and create corresponding masks
    # done by hand first in ds9 and copied over here
    bridge_region = [24.498158,42.619034,27.071083,47.087243,32.801407,51.460385,\
    39.498693,45.62965,42.281712,36.584839,42.142561,32.688612,28.923221,15.573047,\
    18.023204,19.340414,6.5625575,26.277121,18.475598,34.269415]

    bridge_list = []
    for i in range(0,len(bridge_region),2):
        bridge_list.append([int(round(bridge_region[i])),int(round(bridge_region[i+1]))])

    bridge_pn = pg.Polygon(bridge_list)
    bridge_mask = getregionmask(bridge_pn, halpha[0].shape, "bridge region.")

    # -------------------------- #
    north_galaxy_small_region = [9.9895105,46.297203,18.241259,52.348485,23.811189,56.061772,\
    30.068765,54.548951,30.687646,49.11655,23.467366,46.022145,17.553613,37.495338,\
    9.3706294,34.05711,7.6515152,36.395105,7.995338,40.38345]

    north_galaxy_small_list = []
    for i in range(0,len(north_galaxy_small_region),2):
        north_galaxy_small_list.append([int(round(north_galaxy_small_region[i])),int(round(north_galaxy_small_region[i+1]))])

    north_galaxy_small_pn = pg.Polygon(north_galaxy_small_list)
    north_mask = getregionmask(north_galaxy_small_pn, halpha[0].shape, "north galaxy small region.")

    # -------------------------- #
    #north_galaxy_region = [19.203963,58.812354,39.420746,58.881119,38.458042,49.047786,\
    #29.449883,49.254079,25.874126,45.334499,18.378788,39.214452,11.983683,32.200466,\
    #5.2447552,31.168998,6.5512821,41.277389,11.777389,46.572261,16.934732,53.517483]

    #north_galaxy_list = []
    #for i in range(0,len(north_galaxy_region),2):
    #    north_galaxy_list.append([int(round(north_galaxy_region[i])),int(round(north_galaxy_region[i+1]))])

    #north_galaxy_pn = pg.Polygon(north_galaxy_list)
    #north_mask = getregionmask(north_galaxy_pn, halpha[0].shape, "north galaxy region.")

    # -------------------------- #
    south_galaxy_region = [38.320513,44.715618,49.597902,43.271562,52.554779,35.019814,\
    51.660839,27.455711,47.878788,14.80303,37.564103,4.969697,30.618881,7.8578089,\
    33.438228,22.77972,42.0338,34.607226,41.758741,35.157343]

    south_galaxy_list = []
    for i in range(0,len(south_galaxy_region),2):
        south_galaxy_list.append([int(round(south_galaxy_region[i])),int(round(south_galaxy_region[i+1]))])

    south_galaxy_pn = pg.Polygon(south_galaxy_list)
    south_mask = getregionmask(south_galaxy_pn, halpha[0].shape, "south galaxy region.")

    # apply bridge mask
    nii_halpha_withcut_bridge = ma.array(nii_halpha_withcut, mask=bridge_mask)
    oiii_hbeta_withcut_bridge = ma.array(oiii_hbeta_withcut, mask=bridge_mask)
    oi_halpha_withcut_bridge = ma.array(oi_halpha_withcut, mask=bridge_mask)
    sii_halpha_withcut_bridge = ma.array(sii_halpha_withcut, mask=bridge_mask)
    halpha_withcut_bridge = ma.array(halpha_withcut, mask=bridge_mask)
    hbeta_withcut_bridge = ma.array(hbeta_withcut, mask=bridge_mask)
    oiii5007_withcut_bridge = ma.array(oiii5007_withcut, mask=bridge_mask)
    oi6300_withcut_bridge = ma.array(oi6300_withcut, mask=bridge_mask)
    nii6583_withcut_bridge = ma.array(nii6583_withcut, mask=bridge_mask)
    sii_withcut_bridge = ma.array(sii_withcut, mask=bridge_mask)

    # apply north mask
    nii_halpha_withcut_north = ma.array(nii_halpha_withcut, mask=north_mask)
    oiii_hbeta_withcut_north = ma.array(oiii_hbeta_withcut, mask=north_mask)
    oi_halpha_withcut_north = ma.array(oi_halpha_withcut, mask=north_mask)
    sii_halpha_withcut_north = ma.array(sii_halpha_withcut, mask=north_mask)
    halpha_withcut_north = ma.array(halpha_withcut, mask=north_mask)
    hbeta_withcut_north = ma.array(hbeta_withcut, mask=north_mask)
    oiii5007_withcut_north = ma.array(oiii5007_withcut, mask=north_mask)
    oi6300_withcut_north = ma.array(oi6300_withcut, mask=north_mask)
    nii6583_withcut_north = ma.array(nii6583_withcut, mask=north_mask)
    sii_withcut_north = ma.array(sii_withcut, mask=north_mask)

    # apply south mask
    nii_halpha_withcut_south = ma.array(nii_halpha_withcut, mask=south_mask)
    oiii_hbeta_withcut_south = ma.array(oiii_hbeta_withcut, mask=south_mask)
    oi_halpha_withcut_south = ma.array(oi_halpha_withcut, mask=south_mask)
    sii_halpha_withcut_south = ma.array(sii_halpha_withcut, mask=south_mask)
    halpha_withcut_south = ma.array(halpha_withcut, mask=south_mask)
    hbeta_withcut_south = ma.array(hbeta_withcut, mask=south_mask)
    oiii5007_withcut_south = ma.array(oiii5007_withcut, mask=south_mask)
    oi6300_withcut_south = ma.array(oi6300_withcut, mask=south_mask)
    nii6583_withcut_south = ma.array(nii6583_withcut, mask=south_mask)
    sii_withcut_south = ma.array(sii_withcut, mask=south_mask)

    #print len(np.where(np.isfinite(nii_halpha_withcut_bridge.flatten()))[0])

    #print len(np.where((nii_halpha_withcut_bridge.flatten() >= -1) & (nii_halpha_withcut_bridge.flatten() <= 0.5))[0])
    #print len(np.where((oiii_hbeta_withcut_bridge.flatten() >= -1) & (oiii_hbeta_withcut_bridge.flatten() <= 1.0))[0])
    print len(np.unique(oiii_hbeta_withcut_bridge.flatten()))
    print len(np.unique(nii_halpha_withcut_bridge.flatten()))

    print len(np.unique(oiii_hbeta_withcut_north.flatten()))
    print len(np.unique(nii_halpha_withcut_north.flatten()))

    print len(np.unique(oiii_hbeta_withcut_south.flatten()))
    print len(np.unique(nii_halpha_withcut_south.flatten()))

    #print len(np.where(abs(nii_halpha) <= 1)[0])
    #print len(np.where(abs(oiii_hbeta) <= 1)[0])
    #print len(np.where(np.isfinite(nii_halpha.flatten()))[0]), len(np.where(np.isfinite(oiii_hbeta.flatten()))[0])
    #print len(np.where(oi_halpha <= 0)[0])
    #print len(np.where(np.isfinite(oi_halpha.flatten()))[0])

    # plot bpt diagrams
    # BPT with [NII]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{log\left( \frac{[NII]}{H\alpha} \right)}$', fontsize=15)
    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    # use these lines if you want to plot points which are circled if they pass the sig cut
    #ax.plot(nii_halpha, oiii_hbeta, 'o', color='k', markersize=2, markeredgecolor='k')
    #ax.scatter(nii_halpha_withsigcut, oiii_hbeta_withsigcut, s=50, facecolors='None', edgecolors='r')

    # these classifications and the ones on the next plots are taken from Kewley et al 2006, MNRAS, 372, 961
    y_agn_hii_line = 1.3 + 0.61 / (np.arange(-1, 0, 0.01) - 0.05)
    y_liner_seyfert_line = 1.19 + 0.61 / (np.arange(-1, 0.4, 0.01) - 0.47)

    #ax.plot(nii_halpha_withcut, oiii_hbeta_withcut, 'o', color='k', markersize=2, markeredgecolor='k')
    ax.plot(nii_halpha_withcut_bridge, oiii_hbeta_withcut_bridge, 'x', color='r', markersize=6, markeredgecolor='r')
    ax.plot(nii_halpha_withcut_north, oiii_hbeta_withcut_north, 'o', color='g', markersize=2, markeredgecolor='g')
    ax.plot(nii_halpha_withcut_south, oiii_hbeta_withcut_south, 'o', color='b', markersize=2, markeredgecolor='b')
    ax.plot(np.arange(-1, 0, 0.01), y_agn_hii_line, '-', color='k')
    ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')

    ax.set_xlim(-1,0.3)
    ax.set_ylim(-1,1)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

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


    fig.savefig(ipac_taffy_figdir + 'bpt_nii.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # BPT with [OI]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{log\left( \frac{[OI]}{H\alpha} \right)}$', fontsize=15)
    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    y_agn_hii_line = 1.33 + 0.73 / (np.arange(-2.5, -0.8, 0.01) + 0.59)
    y_liner_seyfert_line = 1.30 + 1.18 * np.arange(-1.1, 0, 0.01)

    #ax.plot(oi_halpha_withcut, oiii_hbeta_withcut, 'o', color='k', markersize=2, markeredgecolor='k')
    ax.plot(oi_halpha_withcut_bridge, oiii_hbeta_withcut_bridge, 'x', color='r', markersize=6, markeredgecolor='r')
    ax.plot(oi_halpha_withcut_north, oiii_hbeta_withcut_north, 'o', color='g', markersize=2, markeredgecolor='g')
    ax.plot(oi_halpha_withcut_south, oiii_hbeta_withcut_south, 'o', color='b', markersize=2, markeredgecolor='b')
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
                                         bbox_to_anchor=(0.15, 0.2),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_hiibox)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(ipac_taffy_figdir + 'bpt_oi.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # BPT with [SII]
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{log\left( \frac{[SII]}{H\alpha} \right)}$', fontsize=15)
    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)

    y_agn_hii_line = 1.3 + 0.72 / (np.arange(-1, 0.1, 0.01) - 0.32)
    y_liner_seyfert_line = 0.76 + 1.89 * np.arange(-0.3, 1, 0.01)

    #ax.plot(sii_halpha_withcut, oiii_hbeta_withcut, 'o', color='k', markersize=2, markeredgecolor='k')
    ax.plot(sii_halpha_withcut_bridge, oiii_hbeta_withcut_bridge, 'x', color='r', markersize=6, markeredgecolor='r')
    ax.plot(sii_halpha_withcut_north, oiii_hbeta_withcut_north, 'o', color='g', markersize=2, markeredgecolor='g')
    ax.plot(sii_halpha_withcut_south, oiii_hbeta_withcut_south, 'o', color='b', markersize=2, markeredgecolor='b')
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

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(ipac_taffy_figdir + 'bpt_sii.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # map pixel by pixel maps for quantities that go into the BPT diagram
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(nii6583[0], halpha[0])), cmap='summer', origin='lower', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax)

    vdisp_line = vdisp_line2
    vel_comp = 2
    x = np.arange(vdisp_line.shape[0])
    y = np.arange(vdisp_line.shape[1])
    X, Y = np.meshgrid(x, y)
    #levels = np.array([35,45,55,65,75,85])  # uncomment this line when using velocity dispersion component 1
    levels = np.array([50,75,100,125,150])  # uncomment this line when using velocity dispersion component 2

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    fig.savefig(ipac_taffy_figdir + 'nii_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(oiii5007[0], hbeta[0])), cmap='summer', origin='lower', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    fig.savefig(ipac_taffy_figdir + 'oiii_hbeta_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(oi6300[0], halpha[0])), cmap='summer', origin='lower', vmin=-2.5, vmax=0)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    fig.savefig(ipac_taffy_figdir + 'oi_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(sii, halpha[0])), cmap='summer', origin='lower', vmin=-1, vmax=0.2)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    fig.savefig(ipac_taffy_figdir + 'sii_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)