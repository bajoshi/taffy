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
taffy_extdir = home + "/Desktop/ipac/taffy_lzifu/"

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import bpt_plots as bpt

def plot_pix(blue_cont):

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

def plot_map(ew_map, ew_map_north, ew_map_south):

    # NaN out small EW
    nan_idx = np.where(ew_map < 0.1)
    ew_map[nan_idx] = np.nan

    # combine masks and plot
    comb_mask = (north_mask == 1) & (south_mask == 1)

    ew_map = ma.array(ew_map, mask=comb_mask)

    # make figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    cax = ax.imshow(ew_map, vmin=0.0, vmax=20, cmap='viridis', origin='lower', interpolation='None')
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel(r'$\mathrm{H \beta\ EW\, [\AA]}$')

    ax.minorticks_on()

    plt.show()

    fig.savefig(taffy_extdir + 'figures_stitched_cube/hbeta_ew_map.png', \
        dpi=300, bbox_inches='tight')

    return None

def read_map_and_plot(north_mask, south_mask):

    # load previously saved ew map
    ew_map = np.load(taffy_extdir + 'ew_map.npy')

    # plot
    plot_map(ew_map, north_mask, south_mask)

    return None

def save_map_as_fits(blue_cont_hdu, north_mask, south_mask):

    ew_map = np.load(taffy_extdir + 'ew_map.npy')

    nan_idx = np.where(ew_map < 0.1)
    ew_map[nan_idx] = np.nan

    comb_mask = (north_mask == 1) & (south_mask == 1)
    ew_map[comb_mask] = np.nan

    hdu = fits.PrimaryHDU(data=ew_map, header=blue_cont_hdu.header)
    hdu.writeto(taffy_extdir + 'ew_map.fits', overwrite=True)

    return None

if __name__ == '__main__':
    """
    TO-DO:
    1. You need to somehow be able to measure the 
    average error on the continuum on+near Hbeta.
    This error will then be used to remove spaxels 
    which have low significance on the Hbeta abs measurement.
    1a. Is there a way to get a corresponding error 
    cube for the data cube? Ask Phil.
    """

    # read in lzifu output file
    h = fits.open(taffy_extdir + 'stitched_cube.fits')

    blue_cont = h['B_CONTINUUM'].data

    # mask elements where LZIFU gave NaNs
    region_file = open(taffy_extdir + 'NS_EW.reg')

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

    read_map_and_plot(north_mask, south_mask)
    sys.exit(0)

    # create wavelength array
    # I read these data from the header
    delt = 0.3  # the wav axis is sampled at 0.3A
    blue_wav_start = 4662
    total_blue_res_elem = 2227

    blue_wav_arr = [blue_wav_start + delt*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # find hbeta index in wavelength array
    redshift = 0.0145  # average z
    hbeta_air_wav = 4861.363  # air wavelength
    hbeta_wav = hbeta_air_wav*(1+redshift)
    hbeta_idx = np.argmin(abs(blue_wav_arr - hbeta_wav))

    # make empty ew and cont map
    # fit gaussian to abs dip for each spaxel
    ew_map = np.zeros((58,58))
    snr_map = np.zeros((58,58))

    """
    Try
    for i in range(32, 37):#(0,58):
        for j in range(8, 13):#(0,58):
    
    # To check the high EW region in the north galaxy east.
    """

    # get continuum value at hbeta
    for i in range(0,58):
        for j in range(0,58):

            pix_x = j + 1
            pix_y = i + 1

            #print "On ds9 pixel (x,y):", pix_x, pix_y

            if all(np.isnan(blue_cont[:,i,j])):
                ew_map[i,j] = np.nan
                continue

            # initialize fitting arrays
            left_arr = blue_cont[hbeta_idx-175:hbeta_idx-125, i , j]
            right_arr = blue_cont[hbeta_idx+125:hbeta_idx+175, i , j]
            cont_mean = np.nanmean(np.concatenate([left_arr, right_arr]))

            blue_cont_fit_yarr = blue_cont[hbeta_idx-500:hbeta_idx+500, i , j]
            blue_cont_fit_xarr = np.linspace(hbeta_idx-500, hbeta_idx+500, len(blue_cont_fit_yarr))

            # Meaure SNR for hbeta abs in spaxel
            # I'm doing this before the normalization because
            # teh normalization divides my continuum signal but
            # for the SNR to be accurate the error should 
            # also have been normalized. Since I don't have 
            # the error yet I can't do that.
            snr = np.nanmean(np.concatenate([left_arr, right_arr])) / np.nanstd(np.concatenate([left_arr, right_arr]))
            #print "SNR at ds9 pix X,Y:", snr, pix_x, pix_y
            snr_map[i,j] = snr

            if snr < 5.0:
                # i.e. if the flux around the hbeta line does not have 
                # a significant measurement then do not attempt to measure 
                # the EW. The EW is 0 by default.
                continue

            # normalize the y array by the continuum to get a continuum level of approx. 1
            blue_cont_fit_yarr_norm = blue_cont_fit_yarr / cont_mean

            # actual fitting
            # using astropy gaussian absorption and a straight line model
            # Sometimes there is no absorption and it still tries to fit a gaussian;
            # in these cases a line should provide a better fit.
            #gauss_init = models.GaussianAbsorption1D(amplitude=0.25, mean=hbeta_idx, stddev=10.0)
            # Not using the GaussianAbsorption1D model becuase it is deprecated
            # Instead I'm subtracting a Const1D model from the Gaussian1D model
            # as suggested in the Astropy documentation.

            gauss_abs_init1 = models.Const1D(amplitude=1.0)
            gauss_abs_init2 = models.Gaussian1D(amplitude=0.25, mean=hbeta_idx, stddev=10.0)
            gauss_init = gauss_abs_init1 - gauss_abs_init2
            fit_gauss = fitting.LevMarLSQFitter()
            g = fit_gauss(gauss_init, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            linear_int = models.Linear1D(slope=0.1, intercept=1.0)
            fit_line = fitting.LinearLSQFitter()
            l = fit_line(linear_int, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            comb_int = gauss_init + linear_int
            fit_comb = fitting.LevMarLSQFitter()
            gl = fit_comb(comb_int, blue_cont_fit_xarr, blue_cont_fit_yarr_norm)

            gauss_diff2 = np.sum((g(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)
            line_diff2 = np.sum((l(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)
            comb_diff2 = np.sum((gl(blue_cont_fit_xarr) - blue_cont_fit_yarr_norm)**2)

            # get model parameters
            #print gl.param_names
            comb_amp = gl.parameters[1]
            comb_mean = gl.parameters[2]
            comb_stddev = gl.parameters[3]
            comb_slope = gl.parameters[4]
            comb_intercept = gl.parameters[5]

            if gauss_diff2 < line_diff2:
                # i.e. the chi2 for fitting a gaussian has to be 
                # better than that for a fitting a straight line.
                if comb_diff2 < gauss_diff2:
                    # i.e. the chi2 for fitting a gaussian + straight line
                    # has to be better than fitting just a gaussian by itself.
                    abs_area_analytic = comb_amp * comb_stddev * np.sqrt(2 * np.pi) * delt
                    # Multiplying by delt because the wavelength axis is not sampled
                    # at every 1 Angstrom. It is sampled every 0.3 Angstroms.

            else:
                abs_area_analytic = 0.0

            ew_map[i,j] = abs_area_analytic

            # plot to check fit
            """
            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(np.arange(len(blue_cont[:,i,j])), blue_cont[:,i,j] / cont_mean, color='k')
            ax.plot(blue_cont_fit_xarr, blue_cont_fit_yarr_norm, lw=2)
            ax.plot(blue_cont_fit_xarr, g(blue_cont_fit_xarr), ls='--', color='g', lw=2)
            ax.plot(blue_cont_fit_xarr, l(blue_cont_fit_xarr), ls='--', color='r', lw=2)
            ax.plot(blue_cont_fit_xarr, gl(blue_cont_fit_xarr), ls='--', color='orange', lw=2)

            plt.show()
            plt.clf()
            plt.cla()
            plt.close()

            # Print info
            print "H-beta EW:", "{:.2f}".format(ew_map[i,j]), "at DS9 pixel (x,y)", pix_x, pix_y,
            print "    Continuum mean:", cont_mean,
            print "Amplitude:", "{:.2f}".format(comb_amp),
            print "Std Dev:", "{:.2f}".format(comb_stddev)
            """

    # save map as numpy array
    np.save(taffy_extdir + 'ew_map.npy', ew_map)

    # close file
    h.close()

    # Mask spaxels with low SNR
    val_idx = np.where(snr_map >= 5.0)
    snr_mask = np.ones(snr_map.shape, dtype=bool)
    snr_mask[val_idx] = False
    snr_map = ma.array(snr_map, mask=snr_mask)

    plt.imshow(snr_map, vmin=5, vmax=100, origin='lower', interpolation='None')
    plt.show()

    plot_map(ew_map, north_mask, south_mask)
    save_map_as_fits(h['B_CONTINUUM'], north_mask, south_mask)

    sys.exit(0)