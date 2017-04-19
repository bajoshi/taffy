from __future__ import division

import numpy as np
from astropy.io import fits
from scipy.integrate import simps

import os
import sys

import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/phil/TAFFY/products/'
taffy_dir = '/Users/bhavinjoshi/Desktop/ipac/taffy/'
hcg95c_dir = '/Users/bhavinjoshi/Desktop/ipac/hcg95c/'

sys.path.append(hcg95c_dir + 'codes/')
import imp
chk = imp.load_source('chk', '/Users/bhavinjoshi/Desktop/ipac/hcg95c/codes/check_lzifu_fits.py')

def get_line_array(hdu, line_name, line_side, line_width):

    #  Common lines:
    labels = np.array(['OII3726', 'OII3729', 'NeIII3869', 'Hepsilon', 'Hdelta', 'Hgamma', 'OIII4363', 'Hbeta',\
     'OIII4959', 'OIII5007', 'OI6300', 'OI6364 ', 'NII6548', 'Halpha', 'NII6583', 'SII6716', 'SII6731'])

    wavelengths = np.array([3726.032, 3728.815, 3869.060, 3970.072, 4101.734, 4340.464, 4363.210, 4861.325,\
     4958.911, 5006.843, 6300.304, 6363.7, 6548.040, 6562.800, 6583.460, 6716.440, 6730.810])

    # select line(s) to plot; you must also tell it which pixel(s) you'd like to plot
    redshift = 0.015
    line_idx = np.where(labels == line_name)[0]
    line_wav = wavelengths[line_idx] * (1+redshift)
    lam_low = line_wav - line_width
    lam_up = line_wav + line_width

    lam_b, lam_r = chk.create_lambda_arrays(h)

    if line_side == 'blue':
        lam_low_idx = np.argmin(abs(lam_b - lam_low))
        lam_up_idx = np.argmin(abs(lam_b - lam_up))
        line = hdu['B_LINE'].data
        cont = hdu['B_CONTINUUM'].data

    elif line_side == 'red':
        lam_low_idx = np.argmin(abs(lam_r - lam_low))
        lam_up_idx = np.argmin(abs(lam_r - lam_up))
        line = hdu['R_LINE'].data
        cont = hdu['R_CONTINUUM'].data

    spec = line + cont

    # pix_x and pix_y over here should correspond to x and y shown in ds9
    pix_x = 23
    pix_y = 38

    arr_x = pix_y - 1
    arr_y = pix_x - 1

    line_arr = spec[lam_low_idx:lam_up_idx, arr_x, arr_y]
    wav_arr = np.linspace(lam_low, lam_up, len(line_arr))

    """
    if you want the flux in the line (in the same units as the line) then 
    integrate the area under the line fit from lzifu 
    in this case make sure that hte given with **ONLY** covers the required line
    BE CAREFUL with the width!!
    the line flux varies quite a bit (and not intuitively)
    depending on where the end points of the integration are 
    """
    line_flux = simps(y=line[lam_low_idx:lam_up_idx, arr_x, arr_y], x=wav_arr)

    return line_arr, wav_arr, line_flux

if __name__ == '__main__':
    """
    This code overplots the optical lines and FIR lines and compares the total 
    fluxes within the lines.
    """

    # read in lzifu output file
    h = fits.open(taffy_products + 'TAFFY_2_comp.fits')

    # assign line arrays
    halpha = h['HALPHA'].data
    hbeta = h['HBETA'].data
    nii6583 = h['NII6583'].data
    oiii5007 = h['OIII5007'].data
    oi6300 = h['OI6300'].data
    oi6364 = h['OI6364'].data
    sii6716 = h['SII6716'].data
    sii6731 = h['SII6731'].data

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

    # get line and wavelength array to plot and also the flux in the line
    line_arr, wav_arr, line_flux = get_line_array(h, 'OI6300', 'red', 10)
    print line_flux

    # plot line
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(wav_arr, line_arr)

    plt.show()

    sys.exit(0)