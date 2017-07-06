from __future__ import division

import numpy as np
import numpy.ma as ma
from astropy.io import fits
import Polygon as pg
from scipy.integrate import simps 
from astropy.modeling import models, fitting

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = '/Volumes/Bhavins_backup/ipac/TAFFY/'

import bpt_plots as bpt

def test_pix(blue_cont):

    pix_x = 18
    pix_y = 39
    arr_x = pix_y - 1
    arr_y = pix_x - 1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.arange(len(blue_cont[:,arr_x,arr_y])), blue_cont[:,arr_x,arr_y], color='k')
    plt.show()

    return None

def get_pn(region_list):

    pn_list = []
    for i in range(0,len(region_list),2):
        pn_list.append([int(round(region_list[i])),int(round(region_list[i+1]))])

    return pg.Polygon(pn_list)

if __name__ == '__main__':

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'products_big_cube_velsort/big_cube_2_comp_velsort.fits')

    blue_cont = h['B_CONTINUUM'].data

    # mask elements where LZIFU gave NaNs
    region_file = open(taffy_extdir + 'NS_2017new_xy.reg')

    for line in region_file.readlines()[3:]:

        line_reg = line.split('(')[1].split(')')[1].split(' ')[-1].rstrip()
        region_list = np.array(line.split('(')[1].split(')')[0].split(','))

        if line_reg == 'north':
            north_pn = get_pn(region_list.astype(np.float64))
            north_mask = bpt.getregionmask(north_pn, (58,58), "North galaxy region.")

        elif line_reg == 'south':
            south_pn = get_pn(region_list.astype(np.float64))
            south_mask = bpt.getregionmask(south_pn, (58,58), "South galaxy region.")

    region_file.close()

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
    for i in range(0,58):
        for j in range(0,58):

            pix_x = j + 1
            pix_y = i + 1

            #print "On ds9 pixel (x,y):", pix_x, pix_y

            # checks to see that the pixel is not masked and 
            # that no elements in the fitting array are NaN
            #if all_mask[i,j]:
            #    ew_map[i,j] = np.nan
            #    continue

            if all(np.isnan(blue_cont[:,i,j])):
                ew_map[i,j] = np.nan
                continue

            # initialize fitting arrays
            left_arr = blue_cont[hbeta_idx-175:hbeta_idx-125, i , j]
            right_arr = blue_cont[hbeta_idx+125:hbeta_idx+175, i , j]
            cont_mean = np.nanmean(np.concatenate([left_arr, right_arr]))

            blue_cont_fit_yarr = blue_cont[hbeta_idx-500:hbeta_idx+500, i , j]
            blue_cont_fit_xarr = np.linspace(hbeta_idx-500, hbeta_idx+500, len(blue_cont_fit_yarr))

            # normalize the y array by the continuum to get a continuum level of approx. 1
            blue_cont_fit_yarr_norm = blue_cont_fit_yarr / cont_mean

            # actual fitting
            # using astropy gaussian absorption and a straight line model
            # Sometimes there is no absorption and it still tries to fit a gaussian;
            # in these cases a line should provide a better fit.

            gauss_init = models.GaussianAbsorption1D(amplitude=0.25, mean=hbeta_idx, stddev=10.0)
            fit_gauss = fitting.LevMarLSQFitter()
            g = fit_gauss(gauss_init, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            linear_int = models.Linear1D(slope=0.1, intercept=1.0)
            fit_line = fitting.LevMarLSQFitter()
            l = fit_line(linear_int, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            comb_int = gauss_init + linear_int
            fit_comb = fitting.LevMarLSQFitter()
            gl = fit_comb(comb_int, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            gauss_diff2 = np.sum((g(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)          
            line_diff2 = np.sum((l(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)
            comb_diff2 = np.sum((gl(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)

            # get model parameters
            comb_amp = gl.parameters[0]
            comb_mean = gl.parameters[1]
            comb_stddev = gl.parameters[2]
            comb_slope = gl.parameters[3]
            comb_intercept = gl.parameters[4]

            # plot to check fit
            """
            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(np.arange(len(blue_cont[:,i,j])), blue_cont[:,i,j], color='k')
            ax.plot(blue_cont_fit_xarr, blue_cont_fit_yarr_norm, lw=2)
            ax.plot(blue_cont_fit_xarr, g(blue_cont_fit_xarr), ls='--', color='g', lw=2)
            ax.plot(blue_cont_fit_xarr, l(blue_cont_fit_xarr), ls='--', color='r', lw=2)
            ax.plot(blue_cont_fit_xarr, gl(blue_cont_fit_xarr), ls='--', color='orange', lw=2)

            plt.show()
            """

            if gauss_diff2 < line_diff2:

                if comb_diff2 < gauss_diff2:

                    blue_cont_area_yarr = blue_cont[hbeta_idx-70:hbeta_idx+70, i, j]
                    blue_cont_area_xarr = np.linspace(hbeta_idx-70, hbeta_idx+70, len(blue_cont_area_yarr))
                    a1 = simps(y=gl(blue_cont_area_xarr), x=blue_cont_area_xarr)
                    a2 = cont_mean * (blue_cont_area_xarr[-1] - blue_cont_area_xarr[0])
                    abs_area = a2 - a1

                    abs_area_analytic = comb_amp * comb_stddev * np.sqrt(2 * np.pi)

            else:
                abs_area_analytic = 0.0

            ew_map[i,j] = abs_area_analytic #/ cont_mean

            #print abs_area_analytic, abs_area
            #print cont_mean
            #print ew_map[i,j]
            #sys.exit(0)

    #ew_map = ma.array(ew_map, mask=all_mask)

    ew_map_north = ma.array(ew_map, mask=north_mask)
    ew_map_south = ma.array(ew_map, mask=south_mask)

    plt.imshow(ew_map_north, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.show()

    plt.imshow(ew_map_south, cmap='viridis', origin='lower')
    plt.colorbar()
    plt.show()

    sys.exit(0)