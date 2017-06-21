from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from scipy.optimize import curve_fit
from astropy.modeling import models, fitting

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

import bpt_plots as bpt

def gauss(x, a, x0, sigma):
    return 1 - a*np.exp(-(x-x0)**2/(2*sigma**2))

if __name__ == '__main__':

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'big_cube_2_comp.fits')

    blue_cont = h['B_CONTINUUM'].data

    # mask elements where LZIFU gave NaNs
    region_file = open(taffy_extdir + 'all_possibly_notnan_pixels.reg')
    region_list = np.array(region_file.readlines()[-1].split('(')[1].split(')')[0].split(','))
    region_list = region_list.astype(np.float64)
    region_file.close()

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    region_pn = pg.Polygon(pn_list)
    all_mask = bpt.getregionmask(region_pn, (58,58), "region.")

    # create wavelength array
    # I read these data from the header
    delt = 0.3  # the wav axis is sampled at 0.3A
    blue_wav_start = 4662
    total_blue_res_elem = 2227

    blue_wav_arr = [blue_wav_start + delt*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # find hbeta index in wavelength array
    redshift = 0.0145  # average z
    hbeta_air_wav = 4861.363
    hbeta_wav = hbeta_air_wav*(1+redshift)
    hbeta_idx = np.argmin(abs(blue_wav_arr - hbeta_wav))

    # make empty ew and cont map
    # fit gaussian to abs dip for each spaxel
    ew_map = np.zeros((58,58))

    # get continuum value at hbeta
    for i in range(20,58):
        for j in range(40,58):

            if [i,j] in all_mask:
                ew_map[i,j] = np.nan
                continue

            if all(np.isnan(blue_cont[:,i,j])):
                ew_map[i,j] = np.nan
                continue

            print i, j

            #fig = plt.figure()
            #ax = fig.add_subplot(111)
            #ax.plot(np.arange(len(blue_cont[:,i,j])), blue_cont[:,i,j])
            #plt.show()
            #sys.exit(0)

            left_arr = blue_cont[hbeta_idx-150:hbeta_idx-100, i , j]
            right_arr = blue_cont[hbeta_idx+100:hbeta_idx+150, i , j]
            cont_mean = np.nanmean(np.concatenate([left_arr, right_arr]))
            print cont_mean, hbeta_idx

            blue_cont_fit_yarr = blue_cont[hbeta_idx-50:hbeta_idx+65, i , j]
            blue_cont_fit_xarr = np.linspace(hbeta_idx-50, hbeta_idx+65, len(blue_cont_fit_yarr))

            if any(np.isnan(blue_cont_fit_yarr)):
                ew_map[i,j] = np.nan
                continue

            #popt, pcov = curve_fit(gauss, blue_cont_fit_xarr, blue_cont_fit_yarr, p0=[20, hbeta_idx, 20])

            gauss_init = models.GaussianAbsorption1D(amplitude=0.25, mean=hbeta_idx, stddev=20.0)
            fit_gauss = fitting.LevMarLSQFitter()
            g = fit_gauss(gauss_init, blue_cont_fit_xarr/cont_mean, blue_cont_fit_yarr)

            print fit_gauss.fit_info['message']

            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(blue_cont_fit_xarr, blue_cont_fit_yarr/cont_mean, lw=2)
            ax.plot(blue_cont_fit_xarr, g(blue_cont_fit_xarr), ls='--', color='g', lw=2)
            #ax.plot(blue_cont_fit_xarr, gauss(blue_cont_fit_xarr, *popt), ls='--', color='lightgreen', lw=2)

            plt.show()

            sys.exit(0)

            abs_area = np.sqrt(2 * np.pi) * popt[0] * abs(popt[2])
            ew_map[i,j] = abs_area / cont_mean

    ew_map = ma.array(ew_map, mask=all_mask)

    plt.imshow(ew_map, vmin=0, vmax=5, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.show()
    
    sys.exit(0)