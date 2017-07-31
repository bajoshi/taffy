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
savedir = '/Volumes/Bhavins_backup/ipac/TAFFY/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm

speed_of_light = 299792.458  # km/s

def plot_indices(idx):

    dummy = np.zeros((58,58))
    dummy[idx] = 1.0

    plt.imshow(dummy, origin='lower', cmap='Greys')
    plt.show()

    return None

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
    # above its nucleus and not below the nucleus. The  
    # southern nucleus itself seems to have really broad lines.
    # The northern galaxy mostly has double peaked lines closer to 
    # the HII region. It might also have double peaked lines near 
    # the nucleus and in the NW region but this isn't obvious.
    # The HII region definitely has very obvious double peaked lines.
    # The rest of the bridge doesn't show two peaks in the [OIII] line.
    # I think the rest of the bridge can be fit with a single Gaussian.

    # create wavelength array
    # I read these data from the corresponding headers
    delt_b = 0.3  # i.e. the wav axis is sampled at 0.3A
    blue_wav_start = 4662.0
    total_blue_res_elem = 2227
    blue_res = 1.6

    blue_wav_arr = [blue_wav_start + delt_b*i for i in range(total_blue_res_elem)]
    blue_wav_arr = np.asarray(blue_wav_arr)

    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350
    red_res = 1.5

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # find hbeta index in wavelength array
    redshift = 0.0145  # average z
    linename = 'halpha'
    line_air_wav = 6562.80  # 6562.80 for halpha  # 5006.84 for oiii5007  # 4861.363 for hbeta
    channel = 'red'

    if channel == 'blue':
        wav_arr = blue_wav_arr
        line_comp1 = b_line_comp1
        line_comp2 = b_line_comp2
        line_total = b_line
        obs_data = obs_data_b
        data_res = blue_res
        delt = delt_b

    elif channel == 'red':
        wav_arr = red_wav_arr
        line_comp1 = r_line_comp1
        line_comp2 = r_line_comp2
        line_total = r_line
        obs_data = obs_data_r
        data_res = red_res
        delt = delt_r

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
    pix_x = 50
    pix_y = 33
    arr_x = pix_y - 1
    arr_y = pix_x - 1
    box_size = 3

    comp1_inv_idx = np.zeros((58,58))
    comp2_inv_idx = np.zeros((58,58))

    # start looping
    count = 0
    for i in range(58):  # (arr_x, arr_x + box_size):  # If you want to analyze a 3x3 block enter the pix coords of the low left corner above
        for j in range(58):  # (arr_y, arr_y + box_size):

            # slice array to get only region containing line
            linepad_left = 35
            linepad_right = 50
            # find the center of the biggest peak and call that the line idx
            line_wav = line_air_wav*(1+redshift)
            line_idx = np.argmin(abs(wav_arr - line_wav))
            max_val_idx = np.argmax(obs_data[line_idx-80: line_idx+80, i, j])
            line_idx = line_idx-80+max_val_idx

            line_y_arr_comp1 = line_comp1[line_idx-linepad_left:line_idx+linepad_right, i, j]
            line_x_arr_comp1 = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_comp1))

            line_y_arr_comp2 = line_comp2[line_idx-linepad_left:line_idx+linepad_right, i, j]
            line_x_arr_comp2 = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_comp2))

            # fitting
            gauss_init_lowcomp = models.Gaussian1D(amplitude=5.0, mean=line_idx-10, stddev=5.0)
            gauss_init_highcomp = models.Gaussian1D(amplitude=5.0, mean=line_idx+10, stddev=5.0)
            fit_gauss = fitting.LevMarLSQFitter()

            g1 = fit_gauss(gauss_init_lowcomp, line_x_arr_comp1, line_y_arr_comp1)
            g2 = fit_gauss(gauss_init_highcomp, line_x_arr_comp2, line_y_arr_comp2)

            # also fit raw data by a single gaussian
            line_y_arr_data = obs_data[line_idx-linepad_left:line_idx+linepad_right, i, j]
            g = fit_gauss(gauss_init_lowcomp, line_x_arr_comp1, line_y_arr_data)

            # save in arrays
            amp_comp1[i,j] = g1.parameters[0]
            vel_comp1[i,j] = g1.parameters[1]
            std_comp1[i,j] = g1.parameters[2]

            #print amp_comp1[i,j], vel_comp1[i,j], std_comp1[i,j]

            amp_comp2[i,j] = g2.parameters[0]
            vel_comp2[i,j] = g2.parameters[1]
            std_comp2[i,j] = g2.parameters[2]

            isinvalid_comp1_fit = np.allclose(g1(line_x_arr_comp1), np.zeros(len(g1(line_x_arr_comp1))))
            isinvalid_comp2_fit = np.allclose(g2(line_x_arr_comp2), np.zeros(len(g2(line_x_arr_comp2))))

            if isinvalid_comp1_fit:
                comp1_inv_idx[i,j] = 1.0

            if isinvalid_comp2_fit:
                comp2_inv_idx[i,j] = 1.0

            #print amp_comp2[i,j], vel_comp2[i,j], std_comp2[i,j]

            """
            #print "amp diff", amp_comp2[i,j] - amp_comp1[i,j]
            print "at pixel", i, j
            print "mean diff", (((vel_comp2[i,j] - vel_comp1[i,j]) * 0.3) / line_air_wav) * speed_of_light
            print "std devs", std_comp2[i,j], std_comp1[i,j]
            print "All zeros in comp1", np.allclose(g1(line_x_arr_comp1), np.zeros(len(g1(line_x_arr_comp1))))
            print "All zeros in comp2", np.allclose(g2(line_x_arr_comp2), np.zeros(len(g2(line_x_arr_comp2))))
            print '\n' 

            # plot to check
            fig = plt.figure()
            ax = fig.add_subplot(111)

            ax.plot(line_x_arr_comp1, line_y_arr_comp1, '.', color='b')
            ax.plot(line_x_arr_comp1, g1(line_x_arr_comp1), ls='--', color='b', lw=2)

            ax.plot(line_x_arr_comp2, line_y_arr_comp2, '.', color='r')
            ax.plot(line_x_arr_comp2, g2(line_x_arr_comp2), ls='--', color='r', lw=2)

            ax.plot(line_x_arr_comp1, line_total[line_idx-linepad_left:line_idx+linepad_right, i, j], color='k')

            # also showing the raw data and *MY* single gaussian fit to the raw data
            ax.plot(line_x_arr_comp1, line_y_arr_data, color='gray')
            ax.plot(line_x_arr_comp1, g(line_x_arr_comp1), ls='--', color='g', lw=2)

            plt.show()
            plt.clf()
            plt.cla()
            plt.close()
            """

            count += 1

    # save fit parameters
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

    # now save the different cases separately (see action_items.txt for more details)
    # first convert mean diff arr to units of km/s
    # currently it is sampled at 0.3A so multiplying by that gives the diff in Angstroms
    # convert angstroms to km/s by using line wav and speed of light
    mean_diff *= 0.3
    mean_diff = (mean_diff / line_air_wav) * speed_of_light

    # I'm puting the std values in abs() because I noticed that some std 
    # values are negative. I have no idea why, for now.
    #std_comp1 = np.absolute(std_comp1)
    #std_comp2 = np.absolute(std_comp2)

    # 1. if mean and std are not too different ==> there is only a single comp
    single_idx = np.where((mean_diff < 45) & (std_comp2 < 1.5 * std_comp1))
    plot_indices(single_idx)
    single_idx_arr = np.zeros((58,58))
    single_idx_arr[single_idx] = 1.0

    # 2. Different mean but same std ==> There are two components
    diffmean_idx = np.where((mean_diff >= 45) & (std_comp2 < 1.5 * std_comp1))
    #large_stddiff_idx = np.where(np.absolute(std_comp1 - std_comp2) * delt >= data_res)

    #for k in range(len(large_stddiff_idx[0])):
    #    print large_stddiff_idx[0][k], large_stddiff_idx[1][k]

    #    if (large_stddiff_idx[0][k], large_stddiff_idx[1][k]) in zip(diffmean_idx[0], diffmean_idx[1]):
    #        print "   ", large_stddiff_idx[0][k], large_stddiff_idx[1][k]

    #sys.exit(0)
    plot_indices(diffmean_idx)
    diffmean_idx_arr = np.zeros((58,58))
    diffmean_idx_arr[diffmean_idx] = 1.0

    # 3. Different std but same mean ==> There are two components
    diffstd_idx = np.where((mean_diff < 45) & (std_comp2 >= 1.5 * std_comp1))
    plot_indices(diffstd_idx)
    diffstd_idx_arr = np.zeros((58,58))
    diffstd_idx_arr[diffstd_idx] = 1.0

    # 4. Different mean and std ==> There are two components
    diffboth_idx = np.where((mean_diff >= 45) & (std_comp2 >= 1.5 * std_comp1))
    plot_indices(diffboth_idx)
    diffboth_idx_arr = np.zeros((58,58))
    diffboth_idx_arr[diffboth_idx] = 1.0

    # save all cases
    all_cases_hdu = fits.PrimaryHDU()
    all_cases_hdulist = fits.HDUList(all_cases_hdu)
    all_cases_hdulist.append(fits.ImageHDU(data=comp1_inv_idx))
    all_cases_hdulist.append(fits.ImageHDU(data=comp2_inv_idx))
    all_cases_hdulist.append(fits.ImageHDU(data=single_idx_arr))
    all_cases_hdulist.append(fits.ImageHDU(data=diffmean_idx_arr))
    all_cases_hdulist.append(fits.ImageHDU(data=diffstd_idx_arr))
    all_cases_hdulist.append(fits.ImageHDU(data=diffboth_idx_arr))
    all_cases_hdulist.writeto(savedir + 'all_cases_indices.fits', clobber=True)

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