from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = '/Volumes/Bhavins_backup/ipac/TAFFY/products_big_cube_velsort/'
taffy_data = '/Volumes/Bhavins_backup/ipac/TAFFY/data/'
taffydir = home + '/Desktop/ipac/taffy/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in observed data
    obs_b = fits.open(taffy_data + 'Taffy_B.fits')
    obs_r = fits.open(taffy_data + 'Taffy_R.fits')

    obs_noise_b = obs_b[1].data
    obs_data_b = obs_b[0].data
    obs_noise_r = obs_r[1].data
    obs_data_r = obs_r[0].data

    # read in lzifu output file
    h = fits.open(taffy_products + 'big_cube_2_comp_velsort.fits')

    # assign line arrays
    # Each one of these cubes has the shape 
    # Blue: (2227,58,58)
    # Red: (2350,58,58)
    b_cont = h['B_CONTINUUM'].data
    b_resid = h['B_RESIDFIT'].data
    b_line = h['B_LINE'].data
    b_line_comp1 = h['B_LINE_COMP1'].data
    b_line_comp2 = h['B_LINE_COMP2'].data

    r_cont = h['R_CONTINUUM'].data
    r_resid = h['R_RESIDFIT'].data
    r_line = h['R_LINE'].data
    r_line_comp1 = h['R_LINE_COMP1'].data
    r_line_comp2 = h['R_LINE_COMP2'].data

    # I'm choosing to fit a gaussian to just hte [OIII]5007 line.
    # This is because I think this line has high enough SNR and 
    # also is not contaminated by any other lines. I would've 
    # chosen H-alpha but that has the problem of contamination 
    # from the [NII] satellite lines.

    # For the [OIII]5007 line, this is what I found while looking 
    # at the data in ds9 --
    # It appears that both the northern and southern galaxies have 
    # double peaked profiles. The southern galaxy more so than the 
    # northern one. The southern galaxy has double peaked lines all
    # above its nucleus and not so much below the nucleus. The  
    # southern nucleus itself seems to have really broad lines.
    # The northern galaxy mostly has double peaked lines closer to 
    # the HII region. It might also have double peaked lines near 
    # the nucleus and in the NW region but this isn't obvious.
    # The HII region definitely has very obvious double peaked lines.
    # The rest of the bridge doesn't show two peaks in the [OIII] line.
    # I think the rest of the bridge can be fit with a single Gaussian.

    # create wavelength array
    # I read these data from the header
    delt = 0.3  # the wav axis is sampled at 0.3A
    blue_wav_start = 4662
    total_blue_res_elem = 2227

    blue_wav_arr = [blue_wav_start + delt*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # find hbeta index in wavelength array
    redshift = 0.0145  # average z
    linename = 'oiii5007'
    line_air_wav = 5006.84  # 5006.84 for oiii5007  # 4861.363 for hbeta
    line_wav = line_air_wav*(1+redshift)
    line_idx = np.argmin(abs(blue_wav_arr - line_wav))

    # loop over all spaxels and fit a gaussian to the individual line fits
    # initialize arrays for saving fit parameters
    amp_comp1 = np.ones((58,58)) * -9999.0
    vel_comp1 = np.ones((58,58)) * -9999.0
    std_comp1 = np.ones((58,58)) * -9999.0

    amp_comp2 = np.ones((58,58)) * -9999.0
    vel_comp2 = np.ones((58,58)) * -9999.0
    std_comp2 = np.ones((58,58)) * -9999.0

    # conv ds9 coords to array coords 
    # to be able to check with ds9
    pix_x = 24
    pix_y = 39
    arr_x = pix_y - 1
    arr_y = pix_x - 1

    # start looping
    count = 0
    for i in range(58):  #(arr_x, arr_x + 3):  # If you want to analyze a 3x3 block enter the pix coords of the low left corner above
        for j in range(58):  #(arr_y, arr_y + 3):

            # slice array to get only region containing line
            linepad = 50
            line_y_arr_comp1 = b_line_comp1[line_idx-linepad:line_idx+linepad, i, j]
            line_x_arr_comp1 = np.linspace(line_idx-linepad, line_idx+linepad, len(line_y_arr_comp1))

            line_y_arr_comp2 = b_line_comp2[line_idx-linepad:line_idx+linepad, i, j]
            line_x_arr_comp2 = np.linspace(line_idx-linepad, line_idx+linepad, len(line_y_arr_comp2))

            # fitting
            gauss_init = models.Gaussian1D(amplitude=5.0, mean=line_idx, stddev=5.0)
            fit_gauss = fitting.LevMarLSQFitter()

            g1 = fit_gauss(gauss_init, line_x_arr_comp1, line_y_arr_comp1)
            g2 = fit_gauss(gauss_init, line_x_arr_comp2, line_y_arr_comp2)

            # also fit raw data by a single gaussian
            line_y_arr_single_gaussfit = obs_data_b[line_idx-linepad:line_idx+linepad, i, j]
            g = fit_gauss(gauss_init, line_x_arr_comp1, line_y_arr_single_gaussfit)

            # save in arrays
            amp_comp1[i,j] = g1.parameters[0]
            vel_comp1[i,j] = g1.parameters[1]
            std_comp1[i,j] = g1.parameters[2]

            #print amp_comp1[i,j], vel_comp1[i,j], std_comp1[i,j]

            amp_comp2[i,j] = g2.parameters[0]
            vel_comp2[i,j] = g2.parameters[1]
            std_comp2[i,j] = g2.parameters[2]

            #print amp_comp2[i,j], vel_comp2[i,j], std_comp2[i,j]

            #print "amp diff", amp_comp2[i,j] - amp_comp1[i,j]
            #print "mean diff", vel_comp2[i,j] - vel_comp1[i,j]
            #print "std diff", std_comp2[i,j] - std_comp1[i,j]

            # plot to check
            # comp 1
            """
            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(line_x_arr_comp1, line_y_arr_comp1, '.', color='b')
            ax.plot(line_x_arr_comp1, g1(line_x_arr_comp1), ls='--', color='b', lw=2)

            ax.plot(line_x_arr_comp2, line_y_arr_comp2, '.', color='r')
            ax.plot(line_x_arr_comp2, g2(line_x_arr_comp2), ls='--', color='r', lw=2)

            ax.plot(line_x_arr_comp1, b_line[line_idx-linepad:line_idx+linepad, i, j], color='k')

            # also showing the raw data and *MY* single gaussian fit to the raw data
            ax.plot(line_x_arr_comp1, line_y_arr_single_gaussfit, color='gray')
            ax.plot(line_x_arr_comp1, g(line_x_arr_comp1), ls='--', color='g', lw=2)

            plt.show()
            plt.clf()
            plt.cla()
            plt.close()
            """

            #if count == 2: sys.exit(0)

            count += 1

    # save fit parameters
    savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'
    np.save(savedir + 'amp_' + linename + '_comp1.npy', amp_comp1)
    np.save(savedir + 'vel_' + linename + '_comp1.npy', vel_comp1)
    np.save(savedir + 'std_' + linename + '_comp1.npy', std_comp1)

    np.save(savedir + 'amp_' + linename + '_comp2.npy', amp_comp2)
    np.save(savedir + 'vel_' + linename + '_comp2.npy', vel_comp2)
    np.save(savedir + 'std_' + linename + '_comp2.npy', std_comp2)

    # get mask and set all masked elements to np.nan
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels')

    amp_diff = np.absolute(amp_comp2 - amp_comp1)
    mean_diff = np.absolute(vel_comp2 - vel_comp1)
    std_diff = np.absolute(std_comp2 - std_comp1)

    amp_diff = np.ma.array(amp_diff, mask=all_mask)
    mean_diff = np.ma.array(mean_diff, mask=all_mask)
    std_diff = np.ma.array(std_diff, mask=all_mask)

    amp_diff = np.ma.filled(amp_diff, fill_value=np.nan)
    mean_diff = np.ma.filled(mean_diff, fill_value=np.nan)
    std_diff = np.ma.filled(std_diff, fill_value=np.nan)

    # save differences as fits file
    new_hdu = fits.PrimaryHDU(data=amp_diff)
    new_hdu.writeto(savedir + 'amp_diff_' + linename + '_2m1.fits', clobber=True)
    del new_hdu

    new_hdu = fits.PrimaryHDU(data=mean_diff)
    new_hdu.writeto(savedir + 'mean_diff_' + linename + '_2m1.fits', clobber=True)
    del new_hdu
    
    new_hdu = fits.PrimaryHDU(data=std_diff)
    new_hdu.writeto(savedir + 'std_diff_' + linename + '_2m1.fits', clobber=True)
    del new_hdu

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)