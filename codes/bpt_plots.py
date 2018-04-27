from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg

import os
import sys
import glob
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
ipac_taffy_figdir = home + '/Desktop/ipac/taffy_lzifu/figures_stitched_cube/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'
em_line_dir = taffy_extdir + 'emission_line_ratios/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

def getregionmask(pn, final_arr_shape, regionname):
    """
    This function will return a Boolean mask defining the region of interest.

    Arguments:
    pn: Object of Polygon class. 
    This is the polygon defining the region of interest.

    final_arr_shape: Numpy array shape.
    This is two integers separated by a comma and enclosed in parentheses.
    e.g. for the taffy VIRUS-P IFU data this should be (58,58)
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

def zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, \
    halpha, hbeta, oiii5007, nii6583, oi6300, sii, \
    halpha_err, hbeta_err, oiii5007_err, nii6583_err, oi6300_err, sii_err, orig_shape):
    """
    The arrays that this function expects are of the same shape as orig_shape.
    """

    zipped_ha_ind = zip(ha_valid_ind[0], ha_valid_ind[1])
    zipped_hb_ind = zip(hb_valid_ind[0], hb_valid_ind[1])
    zipped_oiii_ind = zip(oiii_valid_ind[0], oiii_valid_ind[1])
    zipped_nii_ind = zip(nii_valid_ind[0], nii_valid_ind[1])
    zipped_oi_ind = zip(oi_valid_ind[0], oi_valid_ind[1])
    zipped_sii_ind = zip(sii_valid_ind[0], sii_valid_ind[1])

    # creating arrays with all values set to -9999.0 which I'm calling the default null value
    nii_halpha_withcut = np.ones(orig_shape) * -9999.0
    oi_halpha_withcut = np.ones(orig_shape) * -9999.0
    sii_halpha_withcut = np.ones(orig_shape) * -9999.0

    oiii_hbeta_for_nii_withcut = np.ones(orig_shape) * -9999.0
    oiii_hbeta_for_oi_withcut = np.ones(orig_shape) * -9999.0
    oiii_hbeta_for_sii_withcut = np.ones(orig_shape) * -9999.0

    halpha_withcut = np.ones(orig_shape) * -9999.0
    hbeta_withcut = np.ones(orig_shape) * -9999.0
    oiii5007_withcut = np.ones(orig_shape) * -9999.0
    oi6300_withcut = np.ones(orig_shape) * -9999.0
    nii6583_withcut = np.ones(orig_shape) * -9999.0
    sii_withcut = np.ones(orig_shape) * -9999.0
    
    #----Error arrays----#
    nii_halpha_err_withcut = np.ones(orig_shape) * -9999.0
    oi_halpha_err_withcut = np.ones(orig_shape) * -9999.0
    sii_halpha_err_withcut = np.ones(orig_shape) * -9999.0

    oiii_hbeta_for_nii_err_withcut = np.ones(orig_shape) * -9999.0
    oiii_hbeta_for_oi_err_withcut = np.ones(orig_shape) * -9999.0
    oiii_hbeta_for_sii_err_withcut = np.ones(orig_shape) * -9999.0

    # Put error arrays in the form you want them 
    halpha_err = halpha_err / halpha
    nii6583_err = nii6583_err / nii6583
    oi6300_err = oi6300_err / oi6300
    sii_err = sii_err / sii
    oiii5007_err = oiii5007_err / oiii5007
    hbeta_err = hbeta_err / hbeta

    # Loop for [NII]
    for i in range(len(zipped_ha_ind)):

        if zipped_ha_ind[i] in zipped_nii_ind:
            if zipped_ha_ind[i] in zipped_hb_ind:
                if zipped_ha_ind[i] in zipped_oiii_ind:

                    nii_halpha_withcut[zipped_ha_ind[i]] = \
                    np.log10(nii6583[zipped_ha_ind[i]] / halpha[zipped_ha_ind[i]])

                    nii_halpha_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(halpha_err[zipped_ha_ind[i]]**2 + nii6583_err[zipped_ha_ind[i]]**2) / np.log(10)

                    oiii_hbeta_for_nii_withcut[zipped_ha_ind[i]] = \
                    np.log10(oiii5007[zipped_ha_ind[i]] / hbeta[zipped_ha_ind[i]])
                    
                    oiii_hbeta_for_nii_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(oiii5007_err[zipped_ha_ind[i]]**2 + hbeta_err[zipped_ha_ind[i]]**2) / np.log(10)

                    halpha_withcut[zipped_ha_ind[i]] = halpha[zipped_ha_ind[i]]
                    hbeta_withcut[zipped_ha_ind[i]] = hbeta[zipped_ha_ind[i]]
                    oiii5007_withcut[zipped_ha_ind[i]] = oiii5007[zipped_ha_ind[i]]
                    nii6583_withcut[zipped_ha_ind[i]] = nii6583[zipped_ha_ind[i]]

    # Loop again for [OI]
    for i in range(len(zipped_ha_ind)):

        if zipped_ha_ind[i] in zipped_oi_ind:
            if zipped_ha_ind[i] in zipped_hb_ind:
                if zipped_ha_ind[i] in zipped_oiii_ind:

                    oi_halpha_withcut[zipped_ha_ind[i]] = \
                    np.log10(oi6300[zipped_ha_ind[i]] / halpha[zipped_ha_ind[i]])

                    oi_halpha_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(halpha_err[zipped_ha_ind[i]]**2 + oi6300_err[zipped_ha_ind[i]]**2) / np.log(10)

                    oiii_hbeta_for_oi_withcut[zipped_ha_ind[i]] = \
                    np.log10(oiii5007[zipped_ha_ind[i]] / hbeta[zipped_ha_ind[i]])

                    oiii_hbeta_for_oi_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(oiii5007_err[zipped_ha_ind[i]]**2 + hbeta_err[zipped_ha_ind[i]]**2) / np.log(10)

                    halpha_withcut[zipped_ha_ind[i]] = halpha[zipped_ha_ind[i]]
                    hbeta_withcut[zipped_ha_ind[i]] = hbeta[zipped_ha_ind[i]]
                    oiii5007_withcut[zipped_ha_ind[i]] = oiii5007[zipped_ha_ind[i]]
                    oi6300_withcut[zipped_ha_ind[i]] = oi6300[zipped_ha_ind[i]]

    # Loop again for [SII]
    for i in range(len(zipped_ha_ind)):

        if zipped_ha_ind[i] in zipped_hb_ind:
            if zipped_ha_ind[i] in zipped_oiii_ind:
                if zipped_ha_ind[i] in zipped_sii_ind:

                    sii_halpha_withcut[zipped_ha_ind[i]] = \
                    np.log10(sii[zipped_ha_ind[i]] / halpha[zipped_ha_ind[i]])

                    sii_halpha_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(halpha_err[zipped_ha_ind[i]]**2 + sii_err[zipped_ha_ind[i]]**2) / np.log(10)

                    oiii_hbeta_for_sii_withcut[zipped_ha_ind[i]] = \
                    np.log10(oiii5007[zipped_ha_ind[i]] / hbeta[zipped_ha_ind[i]])

                    oiii_hbeta_for_sii_err_withcut[zipped_ha_ind[i]] = \
                    np.sqrt(oiii5007_err[zipped_ha_ind[i]]**2 + hbeta_err[zipped_ha_ind[i]]**2) / np.log(10)

                    halpha_withcut[zipped_ha_ind[i]] = halpha[zipped_ha_ind[i]]
                    hbeta_withcut[zipped_ha_ind[i]] = hbeta[zipped_ha_ind[i]]
                    oiii5007_withcut[zipped_ha_ind[i]] = oiii5007[zipped_ha_ind[i]]
                    sii_withcut[zipped_ha_ind[i]] = sii[zipped_ha_ind[i]]

    return nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut 

def get_arr_withsigcut(sig_cut, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, orig_shape):
    """
    The arrays that this function expects are of the same shape as orig_shape.
    """

    # apply sigma cuts and create arrays for bpt diagram
    ha_valid_ind = np.where((halpha/halpha_err) > sig_cut)
    hb_valid_ind = np.where((hbeta/hbeta_err) > sig_cut)
    oiii_valid_ind = np.where((oiii5007/oiii5007_err) > sig_cut)
    nii_valid_ind = np.where((nii6583/nii6583_err) > sig_cut)
    oi_valid_ind = np.where((oi6300/oi6300_err) > sig_cut)
    sii_valid_ind = np.where((sii/sii_err) > sig_cut)

    nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut = \
    zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, \
    halpha, hbeta, oiii5007, nii6583, oi6300, sii, \
    halpha_err, hbeta_err, oiii5007_err, nii6583_err, oi6300_err, sii_err, \
    orig_shape)

    return nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut

def get_arr_basecut(halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, orig_shape):
    """
    The arrays that this function expects are of the same shape as orig_shape.
    """

    ha_valid_ind = np.where(halpha >= 10)
    hb_valid_ind = np.where(hbeta >= 4)
    oiii_valid_ind = np.where(oiii5007 >= 1)
    nii_valid_ind = np.where(nii6583 >= 8)
    oi_valid_ind = np.where(oi6300 >= 0.5)
    sii_valid_ind = np.where(sii >= 6)

    nii_halpha_withbasecut, oi_halpha_withbasecut, sii_halpha_withbasecut, \
    halpha_withbasecut, hbeta_withbasecut, oiii5007_withbasecut, oi6300_withbasecut, nii6583_withbasecut, sii_withbasecut, \
    oiii_hbeta_for_nii_withbasecut, oiii_hbeta_for_oi_withbasecut, oiii_hbeta_for_sii_withbasecut =\
    zip_and_append(ha_valid_ind, hb_valid_ind, oiii_valid_ind, nii_valid_ind, oi_valid_ind, sii_valid_ind, \
    halpha, hbeta, oiii5007, nii6583, oi6300, sii, orig_shape) 

    return nii_halpha_withbasecut, oi_halpha_withbasecut, sii_halpha_withbasecut, \
    halpha_withbasecut, hbeta_withbasecut, oiii5007_withbasecut, oi6300_withbasecut, nii6583_withbasecut, sii_withbasecut, \
    oiii_hbeta_for_nii_withbasecut, oiii_hbeta_for_oi_withbasecut, oiii_hbeta_for_sii_withbasecut

def check_mappings_field(model, fieldname):

    check_result = False

    if fieldname in model.dtype.names:
        check_result = True

    return check_result

def mappings_oi_nii_sii():

    all_n = ['M_n1_', 'T_n0_01_', 'U_n0_1_', 'V_n10_', 'L_n100_', 'S_n1000_']  # all solar abundances
    # the T model (n=0.01) doesn't have a B=5 microgauss model so it will not be considered,
    # even though it is part of this list.

    mappings_oi_halpha_v100 = []
    mappings_oi_halpha_v125 = []
    mappings_oi_halpha_v150 = []
    mappings_oi_halpha_v175 = []
    mappings_oi_halpha_v200 = []
    mappings_oi_halpha_v225 = []
    mappings_oi_halpha_v250 = []
    mappings_oi_halpha_v300 = []
    mappings_oi_halpha_v500 = []
    mappings_oi_halpha_v800 = []

    mappings_nii_halpha_v100 = []
    mappings_nii_halpha_v125 = []
    mappings_nii_halpha_v150 = []
    mappings_nii_halpha_v175 = []
    mappings_nii_halpha_v200 = []
    mappings_nii_halpha_v225 = []
    mappings_nii_halpha_v250 = []
    mappings_nii_halpha_v300 = []
    mappings_nii_halpha_v350 = []
    mappings_nii_halpha_v400 = []
    mappings_nii_halpha_v450 = []
    mappings_nii_halpha_v500 = []
    mappings_nii_halpha_v800 = []

    mappings_oiii_hbeta_v100 = []
    mappings_oiii_hbeta_v125 = []
    mappings_oiii_hbeta_v150 = []
    mappings_oiii_hbeta_v175 = []
    mappings_oiii_hbeta_v200 = []
    mappings_oiii_hbeta_v225 = []
    mappings_oiii_hbeta_v250 = []
    mappings_oiii_hbeta_v300 = []
    mappings_oiii_hbeta_v350 = []
    mappings_oiii_hbeta_v400 = []
    mappings_oiii_hbeta_v450 = []
    mappings_oiii_hbeta_v500 = []
    mappings_oiii_hbeta_v800 = []

    mappings_sii_halpha_v100 = []
    mappings_sii_halpha_v125 = []
    mappings_sii_halpha_v150 = []
    mappings_sii_halpha_v175 = []
    mappings_sii_halpha_v200 = []
    mappings_sii_halpha_v225 = []
    mappings_sii_halpha_v250 = []
    mappings_sii_halpha_v300 = []
    mappings_sii_halpha_v500 = []
    mappings_sii_halpha_v800 = []

    for curr_n in all_n:
        
        end_string = 'b5_s_lines.txt'
        if 'U_' in curr_n:
            end_string = 'b5_0_s_lines.txt' 

        for fl_name in glob.glob(em_line_dir + curr_n + end_string):

            mappings_model = np.genfromtxt(fl_name, dtype=None, names=True, skip_header=10)  

            mappings_oi_idx = np.where((mappings_model['Atom'] == 'O') & (mappings_model['Species'] == 'I')\
             & (mappings_model['Wavelength'] == '6300.2000'))[0]

            mappings_nii6583_idx = np.where((mappings_model['Atom'] == 'N') & (mappings_model['Species'] == 'II')\
             & (mappings_model['Wavelength'] == '6583.3400'))[0]

            mappings_nii6547_idx = np.where((mappings_model['Atom'] == 'N') & (mappings_model['Species'] == 'II')\
             & (mappings_model['Wavelength'] == '6547.9600'))[0]

            mappings_halpha_idx = np.where((mappings_model['Atom'] == 'H') & (mappings_model['Species'] == 'I')\
             & (mappings_model['Wavelength'] == '6562.8000'))[0]

            mappings_oiii_idx = np.where((mappings_model['Atom'] == 'O') & (mappings_model['Species'] == 'III')\
             & (mappings_model['Wavelength'] == '5006.7700'))[0]

            mappings_hbeta_idx = np.where((mappings_model['Atom'] == 'H') & (mappings_model['Species'] == 'I')\
             & (mappings_model['Wavelength'] == '4861.3200'))[0]

            mappings_sii6716_idx = np.where((mappings_model['Atom'] == 'S') & (mappings_model['Species'] == 'II')\
             & (mappings_model['Wavelength'] == '6716.3100'))[0]

            mappings_sii6731_idx = np.where((mappings_model['Atom'] == 'S') & (mappings_model['Species'] == 'II')\
             & (mappings_model['Wavelength'] == '6730.6800'))[0]

            mappings_oi_line = mappings_model[mappings_oi_idx]
            mappings_halpha_line = mappings_model[mappings_halpha_idx]
            mappings_oiii_line = mappings_model[mappings_oiii_idx]
            mappings_nii6583_line = mappings_model[mappings_nii6583_idx]
            mappings_nii6547_line = mappings_model[mappings_nii6547_idx]
            mappings_sii6716_line = mappings_model[mappings_sii6716_idx]
            mappings_sii6731_line = mappings_model[mappings_sii6731_idx]

            if check_mappings_field(mappings_model, '100'):
                mappings_oi_halpha_v100.append(np.log10(float(mappings_oi_line['100']) / float(mappings_halpha_line['100'])))
                mappings_nii_halpha_v100.append(np.log10((float(mappings_nii6583_line['100']) + float(mappings_nii6547_line['100'])) / float(mappings_halpha_line['100'])))
                mappings_oiii_hbeta_v100.append(np.log10(float(mappings_oiii_line['100'])))
                sii_sum = float(mappings_sii6716_line['100']) + float(mappings_sii6731_line['100'])
                mappings_sii_halpha_v100.append(np.log10(sii_sum / float(mappings_halpha_line['100'])))

            if check_mappings_field(mappings_model, '125'):
                mappings_oi_halpha_v125.append(np.log10(float(mappings_oi_line['125']) / float(mappings_halpha_line['125'])))
                mappings_nii_halpha_v125.append(np.log10((float(mappings_nii6583_line['125']) + float(mappings_nii6547_line['125'])) / float(mappings_halpha_line['125'])))
                mappings_oiii_hbeta_v125.append(np.log10(float(mappings_oiii_line['125'])))
                sii_sum = float(mappings_sii6716_line['125']) + float(mappings_sii6731_line['125'])
                mappings_sii_halpha_v125.append(np.log10(sii_sum / float(mappings_halpha_line['125'])))

            if check_mappings_field(mappings_model, '150'):
                mappings_oi_halpha_v150.append(np.log10(float(mappings_oi_line['150']) / float(mappings_halpha_line['150'])))
                mappings_nii_halpha_v150.append(np.log10((float(mappings_nii6583_line['150']) + float(mappings_nii6547_line['150'])) / float(mappings_halpha_line['150'])))
                mappings_oiii_hbeta_v150.append(np.log10(float(mappings_oiii_line['150'])))
                sii_sum = float(mappings_sii6716_line['150']) + float(mappings_sii6731_line['150'])
                mappings_sii_halpha_v150.append(np.log10(sii_sum / float(mappings_halpha_line['150'])))

            if check_mappings_field(mappings_model, '175'):
                mappings_oi_halpha_v175.append(np.log10(float(mappings_oi_line['175']) / float(mappings_halpha_line['175'])))
                mappings_nii_halpha_v175.append(np.log10((float(mappings_nii6583_line['175']) + float(mappings_nii6547_line['175'])) / float(mappings_halpha_line['175'])))
                mappings_oiii_hbeta_v175.append(np.log10(float(mappings_oiii_line['175'])))
                sii_sum = float(mappings_sii6716_line['175']) + float(mappings_sii6731_line['175'])
                mappings_sii_halpha_v175.append(np.log10(sii_sum / float(mappings_halpha_line['175'])))

            if check_mappings_field(mappings_model, '200'):
                mappings_oi_halpha_v200.append(np.log10(float(mappings_oi_line['200']) / float(mappings_halpha_line['200'])))
                mappings_nii_halpha_v200.append(np.log10((float(mappings_nii6583_line['200']) + float(mappings_nii6547_line['200'])) / float(mappings_halpha_line['200'])))
                mappings_oiii_hbeta_v200.append(np.log10(float(mappings_oiii_line['200'])))
                sii_sum = float(mappings_sii6716_line['200']) + float(mappings_sii6731_line['200'])
                mappings_sii_halpha_v200.append(np.log10(sii_sum / float(mappings_halpha_line['200'])))

            if check_mappings_field(mappings_model, '225'):
                mappings_oi_halpha_v225.append(np.log10(float(mappings_oi_line['225']) / float(mappings_halpha_line['225'])))
                mappings_nii_halpha_v225.append(np.log10((float(mappings_nii6583_line['225']) + float(mappings_nii6547_line['225'])) / float(mappings_halpha_line['225'])))
                mappings_oiii_hbeta_v225.append(np.log10(float(mappings_oiii_line['225'])))
                sii_sum = float(mappings_sii6716_line['225']) + float(mappings_sii6731_line['225'])
                mappings_sii_halpha_v225.append(np.log10(sii_sum / float(mappings_halpha_line['225'])))

            if check_mappings_field(mappings_model, '250'):
                mappings_oi_halpha_v250.append(np.log10(float(mappings_oi_line['250']) / float(mappings_halpha_line['250'])))
                mappings_nii_halpha_v250.append(np.log10((float(mappings_nii6583_line['250']) + float(mappings_nii6547_line['250'])) / float(mappings_halpha_line['250'])))
                mappings_oiii_hbeta_v250.append(np.log10(float(mappings_oiii_line['250'])))
                sii_sum = float(mappings_sii6716_line['250']) + float(mappings_sii6731_line['250'])
                mappings_sii_halpha_v250.append(np.log10(sii_sum / float(mappings_halpha_line['250'])))

            if check_mappings_field(mappings_model, '300'):
                mappings_oi_halpha_v300.append(np.log10(float(mappings_oi_line['300']) / float(mappings_halpha_line['300'])))
                mappings_nii_halpha_v300.append(np.log10((float(mappings_nii6583_line['300']) + float(mappings_nii6547_line['300'])) / float(mappings_halpha_line['300'])))
                mappings_oiii_hbeta_v300.append(np.log10(float(mappings_oiii_line['300'])))
                sii_sum = float(mappings_sii6716_line['300']) + float(mappings_sii6731_line['300'])
                mappings_sii_halpha_v300.append(np.log10(sii_sum / float(mappings_halpha_line['300'])))

            if check_mappings_field(mappings_model, '350'):
                mappings_nii_halpha_v350.append(np.log10((float(mappings_nii6583_line['350']) + float(mappings_nii6547_line['350'])) / float(mappings_halpha_line['350'])))
                mappings_oiii_hbeta_v350.append(np.log10(float(mappings_oiii_line['350'])))

            if check_mappings_field(mappings_model, '400'):
                mappings_nii_halpha_v400.append(np.log10((float(mappings_nii6583_line['400']) + float(mappings_nii6547_line['400'])) / float(mappings_halpha_line['400'])))
                mappings_oiii_hbeta_v400.append(np.log10(float(mappings_oiii_line['400'])))

            if check_mappings_field(mappings_model, '450'):
                mappings_nii_halpha_v450.append(np.log10((float(mappings_nii6583_line['450']) + float(mappings_nii6547_line['450'])) / float(mappings_halpha_line['450'])))
                mappings_oiii_hbeta_v450.append(np.log10(float(mappings_oiii_line['450'])))

            if check_mappings_field(mappings_model, '500'):
                mappings_oi_halpha_v500.append(np.log10(float(mappings_oi_line['500']) / float(mappings_halpha_line['500'])))
                mappings_nii_halpha_v500.append(np.log10((float(mappings_nii6583_line['500']) + float(mappings_nii6547_line['500'])) / float(mappings_halpha_line['500'])))
                mappings_oiii_hbeta_v500.append(np.log10(float(mappings_oiii_line['500'])))
                sii_sum = float(mappings_sii6716_line['500']) + float(mappings_sii6731_line['500'])
                mappings_sii_halpha_v500.append(np.log10(sii_sum / float(mappings_halpha_line['500'])))

            if check_mappings_field(mappings_model, '800'):
                mappings_oi_halpha_v800.append(np.log10(float(mappings_oi_line['800']) / float(mappings_halpha_line['800'])))
                mappings_nii_halpha_v800.append(np.log10((float(mappings_nii6583_line['800']) + float(mappings_nii6547_line['800'])) / float(mappings_halpha_line['800'])))
                mappings_oiii_hbeta_v800.append(np.log10(float(mappings_oiii_line['800'])))
                sii_sum = float(mappings_sii6716_line['800']) + float(mappings_sii6731_line['800'])
                mappings_sii_halpha_v800.append(np.log10(sii_sum / float(mappings_halpha_line['800'])))

    return mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
           mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
           mappings_oi_halpha_v250, mappings_oi_halpha_v300, mappings_oi_halpha_v500,\
           mappings_oi_halpha_v800,\
           mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
           mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
           mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
           mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
           mappings_nii_halpha_v800,\
        mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
        mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
        mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
        mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
        mappings_oiii_hbeta_v800,\
        mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
        mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
        mappings_sii_halpha_v250, mappings_sii_halpha_v300, mappings_sii_halpha_v500,\
        mappings_sii_halpha_v800

def bpt_to_spatial():

    # spatial mapping for interesting bpt points
    all_oiii_hbeta = [oiii_hbeta_withcut_north, oiii_hbeta_withcut_south, oiii_hbeta_withcut_bridge]
    all_oi_halpha  = [oi_halpha_withcut_north, oi_halpha_withcut_south, oi_halpha_withcut_bridge]
    all_nii_halpha  = [nii_halpha_withcut_north, nii_halpha_withcut_south, nii_halpha_withcut_bridge]
    all_sii_halpha  = [sii_halpha_withcut_north, sii_halpha_withcut_south, sii_halpha_withcut_bridge]

    current_bpt_array = np.zeros(oiii_hbeta_withcut.shape)

    oi_agn_hii = lambda x: 1.33 + 0.73 / (x + 0.59)
    oi_liner_seyfert = lambda x: 1.30 + 1.18 * x

    nii_hii = lambda x: 1.3 + 0.61 / (x - 0.05)
    nii_agn = lambda x: 1.19 + 0.61 / (x - 0.47)

    sii_hii = lambda x: 1.30 + 0.72 / (x - 0.32)
    sii_agn = lambda x: 0.76 + 1.89 * x

    current_plot = 'OI'

    for w in range(3):

        current_oiii_hbeta = all_oiii_hbeta[w]

        if current_plot == 'OI':
            current_hii_func = oi_agn_hii
            current_agn_func = oi_liner_seyfert
            current_X_halpha  = all_oi_halpha[w]

        elif current_plot == 'NII':
            current_hii_func = nii_hii
            current_agn_func = nii_agn
            current_X_halpha  = all_nii_halpha[w]

        elif current_plot == 'SII':
            current_hii_func = sii_hii
            current_agn_func = sii_agn
            current_X_halpha  = all_sii_halpha[w]

        for u in range(current_oiii_hbeta.shape[0]):

            for v in range(current_oiii_hbeta.shape[1]):

                current_y = current_oiii_hbeta[u,v]
                current_x = current_X_halpha[u,v]

                # add 1 to both u and v to convert to ds9 y and x respectively

                # for LINER region
                if current_plot == 'OI' or current_plot == 'SII':

                    if (current_y > current_hii_func(current_x)) and (current_y < current_agn_func(current_x)) and (current_x != 0.0) and (current_y != 0.0):
                        #print 'LINER', current_x, current_y, v, u
                        current_bpt_array[u,v] = 1.0

                    # for HII region
                    if (current_y <= current_hii_func(current_x)) and (current_y != 0.0):
                        #print 'HII', current_x, current_y
                        current_bpt_array[u,v] = 2.0

                    # for SEYFERT region
                    if (current_y > current_hii_func(current_x)) and (current_y > current_agn_func(current_x)):
                        #print 'SEYFERT', current_x, current_y
                        current_bpt_array[u,v] = 3.0

                if current_plot == 'NII':

                    if (current_y > nii_agn(current_x)) and (current_y != 0.0):
                        # AGN
                        #print 'AGN', current_x, current_y
                        current_bpt_array[u,v] = 1.0

                    if current_y < nii_hii(current_x):
                        # HII region
                        #print 'HII', current_x, current_y
                        current_bpt_array[u,v] = 2.0

                    if (current_y < nii_agn(current_x)) and (current_y > nii_hii(current_x)):
                        # HII-AGN composite
                        #print 'HII-AGN composite', current_x, current_y
                        current_bpt_array[u,v] = 3.0

    plt.imshow(current_bpt_array, origin='lower')
    plt.colorbar()
    #plt.show()
    #sys.exit(0)

    return None

def getallmasks(shape):

    # define regions and create corresponding masks
    # done by hand first in ds9 and copied over here
    # Use a polygon shaped region and save it in the image coord system
    # so that ds9 will write both x and y coords of the region which you
    # can copy and paste here from the .reg file.
    bridge_region = [32.281801,57.268842,38.538533,57.614395,45.57281,57.702324,47.595164,49.173263,\
    42.495313,39.325275,38.450604,30.884142,35.373108,26.839433,27.793224,19.116218,22.095459,19.749303,\
    15.870124,18.483133,12.599185,23.442299,10.058727,30.487856,14.788552,34.751095,24.463539,45.30257,\
    27.627401,47.586352,33.341603,48.641282,31.758745,54.356923]

    bridge_list = []
    for i in range(0,len(bridge_region),2):
        bridge_list.append([int(round(bridge_region[i])),int(round(bridge_region[i+1]))])

    bridge_pn = pg.Polygon(bridge_list)
    bridge_mask = getregionmask(bridge_pn, shape, "bridge region.")

    # -------------------------- #
    north_galaxy_region = [20.778985,57.420501,32.283659,57.25855,31.756076,54.356913,33.338764,48.641556,\
    27.623409,47.586438,24.457975,45.300312,14.785798,34.748937,9.2463044,34.573104,5.7291782,37.474758,\
    9.2463475,44.684877]

    north_galaxy_list = []
    for i in range(0,len(north_galaxy_region),2):
        north_galaxy_list.append([int(round(north_galaxy_region[i])),int(round(north_galaxy_region[i+1]))])

    north_galaxy_pn = pg.Polygon(north_galaxy_list)
    north_mask = getregionmask(north_galaxy_pn, shape, "north galaxy region.")

    # -------------------------- #
    south_galaxy_region = [47.596266,49.175706,52.888458,43.299641,52.255373,30.743457,47.507236,17.554187,\
    44.447325,11.856422,39.593674,8.4799692,31.785626,9.6406249,24.083093,10.379224,21.550753,15.021847,\
    22.103279,19.736942,27.808648,19.103012,35.374908,26.873766,38.491912,30.927514,42.500817,39.360307]

    south_galaxy_list = []
    for i in range(0,len(south_galaxy_region),2):
        south_galaxy_list.append([int(round(south_galaxy_region[i])),int(round(south_galaxy_region[i+1]))])

    south_galaxy_pn = pg.Polygon(south_galaxy_list)
    south_mask = getregionmask(south_galaxy_pn, shape, "south galaxy region.")

    return bridge_mask, north_mask, south_mask

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # Read in stitched cube
    stitched_cube = fits.open(taffy_extdir + 'stitched_cube.fits')

    # assign line arrays
    halpha = stitched_cube['HALPHA'].data[0]
    hbeta = stitched_cube['HBETA'].data[0]
    nii6583 = stitched_cube['NII6583'].data[0]
    oiii5007 = stitched_cube['OIII5007'].data[0]
    oi6300 = stitched_cube['OI6300'].data[0]
    oi6364 = stitched_cube['OI6364'].data[0]
    sii6716 = stitched_cube['SII6716'].data[0]
    sii6731 = stitched_cube['SII6731'].data[0]

    halpha_err = stitched_cube['HALPHA_ERR'].data[0]
    hbeta_err = stitched_cube['HBETA_ERR'].data[0]
    nii6583_err = stitched_cube['NII6583_ERR'].data[0]
    oiii5007_err = stitched_cube['OIII5007_ERR'].data[0]
    oi6300_err = stitched_cube['OI6300_ERR'].data[0]
    oi6364_err = stitched_cube['OI6364_ERR'].data[0]
    sii6716_err = stitched_cube['SII6716_ERR'].data[0]
    sii6731_err = stitched_cube['SII6731_ERR'].data[0]

    # add lines which are doublets
    sii = sii6716 + sii6731
    sii_err = np.sqrt((sii6716_err)**2 + (sii6731_err)**2)

    # get arrays with some baseline level cut off
    #nii_halpha_withcut, oiii_hbeta_withcut, oi_halpha_withcut, sii_halpha_withcut,\
    #halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut  =\
    # get_arr_basecut(halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    #nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha.shape)

    # get arrays with significance cut applied
    nii_halpha_withcut, oi_halpha_withcut, sii_halpha_withcut, \
    nii_halpha_err_withcut, oi_halpha_err_withcut, sii_halpha_err_withcut, \
    halpha_withcut, hbeta_withcut, oiii5007_withcut, oi6300_withcut, nii6583_withcut, sii_withcut, \
    oiii_hbeta_for_nii_withcut, oiii_hbeta_for_oi_withcut, oiii_hbeta_for_sii_withcut, \
    oiii_hbeta_for_nii_err_withcut, oiii_hbeta_for_oi_err_withcut, oiii_hbeta_for_sii_err_withcut =\
     get_arr_withsigcut(3, halpha, halpha_err, hbeta, hbeta_err, oiii5007, oiii5007_err,\
    nii6583, nii6583_err, oi6300, oi6300_err, sii, sii_err, halpha.shape)

    # get the region masks
    bridge_mask = vcm.get_region_mask('bridge_bpt_new')
    north_mask = vcm.get_region_mask('north_galaxy_bpt')
    south_mask = vcm.get_region_mask('south_galaxy_bpt')

    # apply bridge mask
    nii_halpha_withcut_bridge = ma.array(nii_halpha_withcut, mask=bridge_mask)
    oi_halpha_withcut_bridge = ma.array(oi_halpha_withcut, mask=bridge_mask)
    sii_halpha_withcut_bridge = ma.array(sii_halpha_withcut, mask=bridge_mask)

    oiii_hbeta_for_nii_withcut_bridge = ma.array(oiii_hbeta_for_nii_withcut, mask=bridge_mask)
    oiii_hbeta_for_oi_withcut_bridge = ma.array(oiii_hbeta_for_oi_withcut, mask=bridge_mask)
    oiii_hbeta_for_sii_withcut_bridge = ma.array(oiii_hbeta_for_sii_withcut, mask=bridge_mask)

    # on errors
    nii_halpha_err_withcut_bridge = ma.array(nii_halpha_err_withcut, mask=bridge_mask)
    oi_halpha_err_withcut_bridge = ma.array(oi_halpha_err_withcut, mask=bridge_mask)
    sii_halpha_err_withcut_bridge = ma.array(sii_halpha_err_withcut, mask=bridge_mask)

    oiii_hbeta_for_nii_err_withcut_bridge = ma.array(oiii_hbeta_for_nii_err_withcut, mask=bridge_mask)
    oiii_hbeta_for_oi_err_withcut_bridge = ma.array(oiii_hbeta_for_oi_err_withcut, mask=bridge_mask)
    oiii_hbeta_for_sii_err_withcut_bridge = ma.array(oiii_hbeta_for_sii_err_withcut, mask=bridge_mask)

    halpha_withcut_bridge = ma.array(halpha_withcut, mask=bridge_mask)
    hbeta_withcut_bridge = ma.array(hbeta_withcut, mask=bridge_mask)
    oiii5007_withcut_bridge = ma.array(oiii5007_withcut, mask=bridge_mask)
    oi6300_withcut_bridge = ma.array(oi6300_withcut, mask=bridge_mask)
    nii6583_withcut_bridge = ma.array(nii6583_withcut, mask=bridge_mask)
    sii_withcut_bridge = ma.array(sii_withcut, mask=bridge_mask)

    # apply north mask
    nii_halpha_withcut_north = ma.array(nii_halpha_withcut, mask=north_mask)
    oi_halpha_withcut_north = ma.array(oi_halpha_withcut, mask=north_mask)
    sii_halpha_withcut_north = ma.array(sii_halpha_withcut, mask=north_mask)

    oiii_hbeta_for_nii_withcut_north = ma.array(oiii_hbeta_for_nii_withcut, mask=north_mask)
    oiii_hbeta_for_oi_withcut_north = ma.array(oiii_hbeta_for_oi_withcut, mask=north_mask)
    oiii_hbeta_for_sii_withcut_north = ma.array(oiii_hbeta_for_sii_withcut, mask=north_mask)

    # on errors
    nii_halpha_err_withcut_north = ma.array(nii_halpha_err_withcut, mask=north_mask)
    oi_halpha_err_withcut_north = ma.array(oi_halpha_err_withcut, mask=north_mask)
    sii_halpha_err_withcut_north = ma.array(sii_halpha_err_withcut, mask=north_mask)

    oiii_hbeta_for_nii_err_withcut_north = ma.array(oiii_hbeta_for_nii_err_withcut, mask=north_mask)
    oiii_hbeta_for_oi_err_withcut_north = ma.array(oiii_hbeta_for_oi_err_withcut, mask=north_mask)
    oiii_hbeta_for_sii_err_withcut_north = ma.array(oiii_hbeta_for_sii_err_withcut, mask=north_mask)

    halpha_withcut_north = ma.array(halpha_withcut, mask=north_mask)
    hbeta_withcut_north = ma.array(hbeta_withcut, mask=north_mask)
    oiii5007_withcut_north = ma.array(oiii5007_withcut, mask=north_mask)
    oi6300_withcut_north = ma.array(oi6300_withcut, mask=north_mask)
    nii6583_withcut_north = ma.array(nii6583_withcut, mask=north_mask)
    sii_withcut_north = ma.array(sii_withcut, mask=north_mask)

    # apply south mask
    nii_halpha_withcut_south = ma.array(nii_halpha_withcut, mask=south_mask)
    oi_halpha_withcut_south = ma.array(oi_halpha_withcut, mask=south_mask)
    sii_halpha_withcut_south = ma.array(sii_halpha_withcut, mask=south_mask)

    oiii_hbeta_for_nii_withcut_south = ma.array(oiii_hbeta_for_nii_withcut, mask=south_mask)
    oiii_hbeta_for_oi_withcut_south = ma.array(oiii_hbeta_for_oi_withcut, mask=south_mask)
    oiii_hbeta_for_sii_withcut_south = ma.array(oiii_hbeta_for_sii_withcut, mask=south_mask)

    # on errors
    nii_halpha_err_withcut_south = ma.array(nii_halpha_err_withcut, mask=south_mask)
    oi_halpha_err_withcut_south = ma.array(oi_halpha_err_withcut, mask=south_mask)
    sii_halpha_err_withcut_south = ma.array(sii_halpha_err_withcut, mask=south_mask)

    oiii_hbeta_for_nii_err_withcut_south = ma.array(oiii_hbeta_for_nii_err_withcut, mask=south_mask)
    oiii_hbeta_for_oi_err_withcut_south = ma.array(oiii_hbeta_for_oi_err_withcut, mask=south_mask)
    oiii_hbeta_for_sii_err_withcut_south = ma.array(oiii_hbeta_for_sii_err_withcut, mask=south_mask)

    halpha_withcut_south = ma.array(halpha_withcut, mask=south_mask)
    hbeta_withcut_south = ma.array(hbeta_withcut, mask=south_mask)
    oiii5007_withcut_south = ma.array(oiii5007_withcut, mask=south_mask)
    oi6300_withcut_south = ma.array(oi6300_withcut, mask=south_mask)
    nii6583_withcut_south = ma.array(nii6583_withcut, mask=south_mask)
    sii_withcut_south = ma.array(sii_withcut, mask=south_mask)

    # read in Mappings III models and overplot
    mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
    mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
    mappings_oi_halpha_v250, mappings_oi_halpha_v300, mappings_oi_halpha_v500,\
    mappings_oi_halpha_v800,\
    mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
    mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
    mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
    mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
    mappings_nii_halpha_v800,\
    mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
    mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
    mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
    mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
    mappings_oiii_hbeta_v800,\
    mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
    mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
    mappings_sii_halpha_v250, mappings_sii_halpha_v300, mappings_sii_halpha_v500,\
    mappings_sii_halpha_v800 = mappings_oi_nii_sii()

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

    nii_nonzero = np.nonzero(nii_halpha_withcut)

    #ax.plot(nii_halpha_withcut_bridge[nii_nonzero], oiii_hbeta_for_nii_withcut_bridge[nii_nonzero], \
    #    'x', color='maroon', markersize=8, markeredgecolor='maroon')
    #ax.plot(nii_halpha_withcut_north[nii_nonzero], oiii_hbeta_for_nii_withcut_north[nii_nonzero], \
    #    'o', color='goldenrod', markersize=3, markeredgecolor='None')
    #ax.plot(nii_halpha_withcut_south[nii_nonzero], oiii_hbeta_for_nii_withcut_south[nii_nonzero], \
    #    'o', color='midnightblue', markersize=3, markeredgecolor='None')

    ax.errorbar(nii_halpha_withcut_bridge[nii_nonzero], oiii_hbeta_for_nii_withcut_bridge[nii_nonzero], \
        xerr=nii_halpha_err_withcut_bridge[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_bridge[nii_nonzero], \
        color='maroon', markersize=8, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(nii_halpha_withcut_north[nii_nonzero], oiii_hbeta_for_nii_withcut_north[nii_nonzero], \
        xerr=nii_halpha_err_withcut_north[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_north[nii_nonzero], \
        color='goldenrod', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(nii_halpha_withcut_south[nii_nonzero], oiii_hbeta_for_nii_withcut_south[nii_nonzero], \
        xerr=nii_halpha_err_withcut_south[nii_nonzero], yerr=oiii_hbeta_for_nii_err_withcut_south[nii_nonzero], \
        color='midnightblue', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

    ax.plot(np.arange(-1, 0, 0.01), y_agn_hii_line, '-', color='k')
    ax.plot(np.arange(-1, 0.4, 0.01), y_liner_seyfert_line, '--', color='k')

    ax.plot(mappings_nii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s', zorder=10)
    ax.plot(mappings_nii_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s', zorder=10)

    ax.legend(loc=0, prop={'size':10})

    ax.set_xlim(-1,0.3)
    ax.set_ylim(-1,1)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    # region labels
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
                                         bbox_to_anchor=(0.32, 0.3),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_hiibox)


    fig.savefig(ipac_taffy_figdir + 'bpt_nii_no_thresh_full_errbar.eps', dpi=300, bbox_inches='tight')
    #fig.savefig(ipac_taffy_figdir + 'bpt_nii_no_thresh.eps', dpi=300, bbox_inches='tight')

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

    oi_nonzero = np.nonzero(oi_halpha_withcut)

    #ax.plot(oi_halpha_withcut_bridge[oi_nonzero], oiii_hbeta_for_oi_withcut_bridge[oi_nonzero], \
    #    'x', color='maroon', markersize=8, markeredgecolor='maroon')
    #ax.plot(oi_halpha_withcut_north[oi_nonzero], oiii_hbeta_for_oi_withcut_north[oi_nonzero], \
    #    'o', color='goldenrod', markersize=3, markeredgecolor='None')
    #ax.plot(oi_halpha_withcut_south[oi_nonzero], oiii_hbeta_for_oi_withcut_south[oi_nonzero], \
    #    'o', color='midnightblue', markersize=3, markeredgecolor='None')

    ax.errorbar(oi_halpha_withcut_bridge[oi_nonzero], oiii_hbeta_for_oi_withcut_bridge[oi_nonzero], \
        xerr=oi_halpha_err_withcut_bridge[oi_nonzero], yerr=oiii_hbeta_for_oi_err_withcut_bridge[oi_nonzero], \
        color='maroon', markersize=8, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(oi_halpha_withcut_north[oi_nonzero], oiii_hbeta_for_oi_withcut_north[oi_nonzero], \
        xerr=oi_halpha_err_withcut_north[oi_nonzero], yerr=oiii_hbeta_for_oi_err_withcut_north[oi_nonzero], \
        color='goldenrod', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(oi_halpha_withcut_south[oi_nonzero], oiii_hbeta_for_oi_withcut_south[oi_nonzero], \
        xerr=oi_halpha_err_withcut_south[oi_nonzero], yerr=oiii_hbeta_for_oi_err_withcut_south[oi_nonzero], \
        color='midnightblue', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

    ax.plot(np.arange(-2.5, -0.8, 0.01), y_agn_hii_line, '-', color='k')
    ax.plot(np.arange(-1.1, 0, 0.01), y_liner_seyfert_line, '--', color='k')

    ax.plot(mappings_oi_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s', zorder=10)
    ax.plot(mappings_oi_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s', zorder=10)
    ax.plot(mappings_oi_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s', zorder=10)
    ax.plot(mappings_oi_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s', zorder=10)
    ax.plot(mappings_oi_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s', zorder=10)
    ax.plot(mappings_oi_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s', zorder=10)

    ax.legend(loc=0, prop={'size':10})

    ax.set_xlim(-2.0,0)
    ax.set_ylim(-1,1)

    # region labels
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

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(ipac_taffy_figdir + 'bpt_oi_no_thresh_full_errbar.eps', dpi=300, bbox_inches='tight')
    #fig.savefig(ipac_taffy_figdir + 'bpt_oi_no_thresh.eps', dpi=300, bbox_inches='tight')

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

    sii_nonzero = np.nonzero(sii_halpha_withcut)

    #ax.plot(sii_halpha_withcut_bridge[sii_nonzero], oiii_hbeta_for_sii_withcut_bridge[sii_nonzero], \
    #    'x', color='maroon', markersize=8, markeredgecolor='maroon')
    #ax.plot(sii_halpha_withcut_north[sii_nonzero], oiii_hbeta_for_sii_withcut_north[sii_nonzero], \
    #    'o', color='goldenrod', markersize=3, markeredgecolor='None')
    #ax.plot(sii_halpha_withcut_south[sii_nonzero], oiii_hbeta_for_sii_withcut_south[sii_nonzero], \
    #    'o', color='midnightblue', markersize=3, markeredgecolor='None')

    ax.errorbar(sii_halpha_withcut_bridge[sii_nonzero], oiii_hbeta_for_sii_withcut_bridge[sii_nonzero], \
        xerr=sii_halpha_err_withcut_bridge[sii_nonzero], yerr=oiii_hbeta_for_sii_err_withcut_bridge[sii_nonzero], \
        color='maroon', markersize=8, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.25)
    ax.errorbar(sii_halpha_withcut_north[sii_nonzero], oiii_hbeta_for_sii_withcut_north[sii_nonzero], \
        xerr=sii_halpha_err_withcut_north[sii_nonzero], yerr=oiii_hbeta_for_sii_err_withcut_north[sii_nonzero], \
        color='goldenrod', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)
    ax.errorbar(sii_halpha_withcut_south[sii_nonzero], oiii_hbeta_for_sii_withcut_south[sii_nonzero], \
        xerr=sii_halpha_err_withcut_south[sii_nonzero], yerr=oiii_hbeta_for_sii_err_withcut_south[sii_nonzero], \
        color='midnightblue', markersize=3, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.25)

    ax.plot(np.arange(-1, 0.1, 0.01), y_agn_hii_line, '-', color='k')
    ax.plot(np.arange(-0.3, 1, 0.01), y_liner_seyfert_line, '--', color='k')

    ax.plot(mappings_sii_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s', zorder=10)
    ax.plot(mappings_sii_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s', zorder=10)
    ax.plot(mappings_sii_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s', zorder=10)
    ax.plot(mappings_sii_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s', zorder=10)
    ax.plot(mappings_sii_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s', zorder=10)
    ax.plot(mappings_sii_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s', zorder=10)

    ax.legend(loc=0, prop={'size':10})

    ax.set_xlim(-1,0.5)
    ax.set_ylim(-1,1)

    # region labels
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

    fig.savefig(ipac_taffy_figdir + 'bpt_sii_no_thresh_full_errbar.eps', dpi=300, bbox_inches='tight')
    #fig.savefig(ipac_taffy_figdir + 'bpt_sii_no_thresh.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()
    sys.exit(0)

    # map pixel by pixel maps for quantities that go into the BPT diagram
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(nii6583, halpha)), cmap='summer', origin='lower', vmin=-1, vmax=1)
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

    #fig.savefig(ipac_taffy_figdir + 'nii_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(oiii5007, hbeta)), cmap='summer', origin='lower', vmin=-1, vmax=1)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    #fig.savefig(ipac_taffy_figdir + 'oiii_hbeta_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(oi6300, halpha)), cmap='summer', origin='lower', vmin=-2.5, vmax=0)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    #fig.savefig(ipac_taffy_figdir + 'oi_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------ #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    im = ax.imshow(np.log10(np.divide(sii, halpha)), cmap='summer', origin='lower', vmin=-1, vmax=0.2)
    plt.colorbar(im, ax=ax)

    c = ax.contour(X, Y, vdisp_line, levels=levels, cmap='inferno')
    ax.clabel(c, inline=True, inline_spacing=1, fontsize=4, fmt='%1.2f')

    #fig.savefig(ipac_taffy_figdir + 'sii_halpha_map_velodisp_comp' + str(vel_comp) + '.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)