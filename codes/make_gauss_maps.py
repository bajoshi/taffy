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
taffy_products = home + '/Desktop/ipac/taffy_lzifu/products_work/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
taffydir = home + '/Desktop/ipac/taffy/'
savedir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

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
    h = fits.open(taffy_extdir + 'Taffy_2_comp_patched.fits')

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

    # find line index in wavelength array
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

    amp_onecomp = np.ones((58,58)) * -9999.0
    vel_onecomp = np.ones((58,58)) * -9999.0
    std_onecomp = np.ones((58,58)) * -9999.0

    # conv ds9 coords to array coords 
    # to be able to check with ds9
    pix_x = 43
    pix_y = 20
    arr_x = pix_y - 1
    arr_y = pix_x - 1
    box_size = 3

    # I'm forcing some low SNR spaxels to get the one comp fit
    force_onecomp_list = [[28, 20],[29, 20],[30, 20],[29, 21],[30, 21],[29, 22],[30, 22],\
    [38,11]]
    force_onecomp_arr_x = []
    force_onecomp_arr_y = []
    for u in range(len(force_onecomp_list)):
        force_onecomp_arr_x.append(force_onecomp_list[u][1] - 1)
        force_onecomp_arr_y.append(force_onecomp_list[u][0] - 1)

    force_onecomp_arr = zip(force_onecomp_arr_x, force_onecomp_arr_y)

    # set up arrays to flag invalid fits
    comp1_inv_idx = np.zeros((58,58))
    comp2_inv_idx = np.zeros((58,58))

    # start looping
    count = 0
    for i in range(arr_x, arr_x + box_size):  # If you want to analyze a block enter the pix coords of the low left corner above
        for j in range(arr_y, arr_y + box_size):

            # find the center of the biggest peak and call that the line idx
            line_wav = line_air_wav*(1+redshift)
            line_idx = np.argmin(abs(wav_arr - line_wav))
            max_val_idx = np.argmax(obs_data[line_idx-45: line_idx+45, i, j])
            # it searches 45 spectral elements on each side of line idx 
            # for the peak because 45 spectral elements at Halpha corresponds to
            # about 600 km/s which is the max velocity I'm willing to believe.
            # Also, another issue with giving it a larger width to search for the
            # peak runs into the risk of confusing NII6584 with Halpha and 
            # consequently getting a wrong high velocity.
            line_idx = line_idx-45+max_val_idx

            # generate arrays that will be given to fitting function
            # slice array to get only region containing line
            linepad_left = 50
            linepad_right = 50

            # some spaxels need special treatment
            if i == 35 and j == 30:
                linepad_left = 25
            if i == 47 and j == 26:
                linepad_left = 75
            if i == 48 and j == 26:
                linepad_left = 75

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

            # save lzifu total fit to array for plotting
            line_y_arr_total = line_total[line_idx-linepad_left:line_idx+linepad_right, i, j]

            # also fit raw data by a single gaussian
            if (i,j) in force_onecomp_arr:
                linepad_left = 25
                linepad_right = 25

            line_y_arr_data = obs_data[line_idx-linepad_left:line_idx+linepad_right, i, j]
            line_x_arr_data = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_data))

            # find pseudo continuum and subtract
            pseudo_cont_arr_left = obs_data[line_idx-10-(linepad_left-10):line_idx-10, i, j]
            pseudo_cont_arr_right = obs_data[line_idx+10:line_idx+10+(linepad_right-10), i, j]
            cont_level = np.nanmean(np.concatenate((pseudo_cont_arr_left, pseudo_cont_arr_right)))

            line_y_arr_data -= cont_level  # subtract the continuum level before trying to fit the single comp gaussian directly to the data
            gauss_init_onecomp = models.Gaussian1D(amplitude=np.nanmean(line_y_arr_data), mean=line_idx-10, stddev=5.0)
            g = fit_gauss(gauss_init_onecomp, line_x_arr_data, line_y_arr_data)

            # make sure that the component you are calling the second component
            # is indeed at a higher velocity that the first component. This is 
            # how LZIFU sorts them i.e. by increasing order of velocity.
            # This has to be done before the numbers are saved in arrays.
            phys_vel_comp1 = float(format(((g1.parameters[1] * 0.3 + red_wav_start - line_air_wav) / line_air_wav) * speed_of_light - speed_of_light*redshift, '.2f'))
            phys_vel_comp2 = float(format(((g2.parameters[1] * 0.3 + red_wav_start - line_air_wav) / line_air_wav) * speed_of_light - speed_of_light*redshift, '.2f'))
            if phys_vel_comp2 < phys_vel_comp1:
                # if it found a low vel fit for the fit for comp 2
                # then switch the fits for 1 and 2 components
                amp_comp1[i,j] = g2.parameters[0]
                vel_comp1[i,j] = g2.parameters[1]
                std_comp1[i,j] = abs(g2.parameters[2])

                amp_comp2[i,j] = g1.parameters[0]
                vel_comp2[i,j] = g1.parameters[1]
                std_comp2[i,j] = abs(g1.parameters[2])

            elif phys_vel_comp2 >= phys_vel_comp1:
                amp_comp1[i,j] = g1.parameters[0]
                vel_comp1[i,j] = g1.parameters[1]
                std_comp1[i,j] = abs(g1.parameters[2])

                amp_comp2[i,j] = g2.parameters[0]
                vel_comp2[i,j] = g2.parameters[1]
                std_comp2[i,j] = abs(g2.parameters[2])

            amp_onecomp[i,j] = g.parameters[0]
            vel_onecomp[i,j] = g.parameters[1]
            std_onecomp[i,j] = abs(g.parameters[2])

            # the fit is invalid if either one of hte line fits gave all zeros
            isinvalid_comp1_fit = np.allclose(g1(line_x_arr_comp1), np.zeros(len(g1(line_x_arr_comp1))), rtol=1e-5, atol=1e-3)
            isinvalid_comp2_fit = np.allclose(g2(line_x_arr_comp2), np.zeros(len(g2(line_x_arr_comp2))), rtol=1e-5, atol=1e-3)

            if isinvalid_comp1_fit:
                comp1_inv_idx[i,j] = 1.0

            if isinvalid_comp2_fit:
                comp2_inv_idx[i,j] = 1.0

            # uncomment the follwing block to run checks and 
            # add sys.exit(0) right after for loop is done
            # also uncomment the for loop range
            print "amp diff", amp_comp2[i,j] - amp_comp1[i,j]
            print "at pixel", j+1, i+1
            print "line idx and center", line_idx, line_idx * 0.3 + red_wav_start
            print "mean vel 1 and 2 [km/s]", phys_vel_comp1, phys_vel_comp2
            print "mean diff", format((((vel_comp2[i,j] - vel_comp1[i,j]) * 0.3) / line_air_wav) * speed_of_light, '.2f'), "km/s"
            print "std devs", format(std_comp2[i,j], '.2f'), format(std_comp1[i,j], '.2f')
            print "continuum level, amp and vdisp for onecomp fit", format(cont_level, '.2f'), format(amp_onecomp[i,j], '.2f'), format(std_onecomp[i,j], '.2f')
            print "All zeros in comp1", np.allclose(g1(line_x_arr_comp1), np.zeros(len(g1(line_x_arr_comp1))), rtol=1e-5, atol=1e-3)
            print "All zeros in comp2", np.allclose(g2(line_x_arr_comp2), np.zeros(len(g2(line_x_arr_comp2))), rtol=1e-5, atol=1e-3)
            print '\n' 

            # plot to check
            fig = plt.figure()
            ax = fig.add_subplot(111)

            # plot first comp fit
            ax.plot(line_x_arr_comp1, line_y_arr_comp1, '.', color='b')
            ax.plot(line_x_arr_comp1, g1(line_x_arr_comp1), ls='--', color='b', lw=2)

            # plot second comp fit
            ax.plot(line_x_arr_comp2, line_y_arr_comp2, '.', color='r')
            ax.plot(line_x_arr_comp2, g2(line_x_arr_comp2), ls='--', color='r', lw=2)

            # plot lzifu total fit
            ax.plot(line_x_arr_comp1, line_y_arr_total, color='k')

            # also show the raw data and *MY* single gaussian fit to the raw data
            ax.plot(line_x_arr_data, line_y_arr_data, color='gray')
            ax.plot(line_x_arr_data, g(line_x_arr_data), ls='--', color='g', lw=2)

            plt.show()
            plt.clf()
            plt.cla()
            plt.close()

            count += 1

    sys.exit(0)

    # save fit parameters
    np.save(savedir + 'amp_' + linename + '_comp1.npy', amp_comp1)
    np.save(savedir + 'vel_' + linename + '_comp1.npy', vel_comp1)
    np.save(savedir + 'std_' + linename + '_comp1.npy', std_comp1)

    np.save(savedir + 'amp_' + linename + '_comp2.npy', amp_comp2)
    np.save(savedir + 'vel_' + linename + '_comp2.npy', vel_comp2)
    np.save(savedir + 'std_' + linename + '_comp2.npy', std_comp2)

    np.save(savedir + 'amp_' + linename + '_onecomp.npy', amp_onecomp)
    np.save(savedir + 'vel_' + linename + '_onecomp.npy', vel_onecomp)
    np.save(savedir + 'std_' + linename + '_onecomp.npy', std_onecomp)

    # get mask and set all masked elements to np.nan
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

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
    single_idx = np.where((mean_diff < 35) & (std_comp2 < 1.5 * std_comp1) & (std_comp1 < 1.5 * std_comp2))
    #plot_indices(single_idx)
    single_idx_arr = np.zeros((58,58))
    single_idx_arr[single_idx] = 1.0

    # 2. Different mean but same std ==> There are two components
    diffmean_idx = np.where((mean_diff >= 35) & (std_comp2 < 1.5 * std_comp1) & (std_comp1 < 1.5 * std_comp2))
    #large_stddiff_idx = np.where(np.absolute(std_comp1 - std_comp2) * delt >= data_res)

    #for k in range(len(large_stddiff_idx[0])):
    #    print large_stddiff_idx[0][k], large_stddiff_idx[1][k]

    #    if (large_stddiff_idx[0][k], large_stddiff_idx[1][k]) in zip(diffmean_idx[0], diffmean_idx[1]):
    #        print "   ", large_stddiff_idx[0][k], large_stddiff_idx[1][k]

    #sys.exit(0)
    #plot_indices(diffmean_idx)
    diffmean_idx_arr = np.zeros((58,58))
    diffmean_idx_arr[diffmean_idx] = 1.0

    # 3. Different std but same mean ==> There are two components
    diffstd_idx = np.where((mean_diff < 35) & ((std_comp2 >= 1.5 * std_comp1) | (std_comp1 >= 1.5 * std_comp2)))
    #plot_indices(diffstd_idx)
    diffstd_idx_arr = np.zeros((58,58))
    diffstd_idx_arr[diffstd_idx] = 1.0

    # 4. Different mean and std ==> There are two components
    diffboth_idx = np.where((mean_diff >= 35) & ((std_comp2 >= 1.5 * std_comp1) | (std_comp1 >= 1.5 * std_comp2)))
    #plot_indices(diffboth_idx)
    diffboth_idx_arr = np.zeros((58,58))
    diffboth_idx_arr[diffboth_idx] = 1.0

    # save all cases
    all_cases_hdu = fits.PrimaryHDU()
    hdr = fits.Header()
    all_cases_hdulist = fits.HDUList(all_cases_hdu)
    hdr['EXTNAME'] = 'COMP1_INV'
    all_cases_hdulist.append(fits.ImageHDU(data=comp1_inv_idx, header=hdr))
    hdr['EXTNAME'] = 'COMP2_INV'
    all_cases_hdulist.append(fits.ImageHDU(data=comp2_inv_idx, header=hdr))
    hdr['EXTNAME'] = 'SINGLE_IDX'
    all_cases_hdulist.append(fits.ImageHDU(data=single_idx_arr, header=hdr))
    hdr['EXTNAME'] = 'DIFFMEAN_IDX'
    all_cases_hdulist.append(fits.ImageHDU(data=diffmean_idx_arr, header=hdr))
    hdr['EXTNAME'] = 'DIFFSTD_IDX'
    all_cases_hdulist.append(fits.ImageHDU(data=diffstd_idx_arr, header=hdr))
    hdr['EXTNAME'] = 'DIFFBOTH_IDX'
    all_cases_hdulist.append(fits.ImageHDU(data=diffboth_idx_arr, header=hdr))
    all_cases_hdulist.writeto(savedir + 'all_cases_indices.fits', overwrite=True)

    h.close()
    sys.exit(0)

    # save differences as fits file
    new_hdu = fits.PrimaryHDU(data=amp_diff)
    new_hdu.writeto(savedir + 'amp_diff_' + linename + '_2m1.fits', overwrite=True)
    del new_hdu

    new_hdu = fits.PrimaryHDU(data=mean_diff)
    new_hdu.writeto(savedir + 'mean_diff_' + linename + '_2m1.fits', overwrite=True)
    del new_hdu
    
    new_hdu = fits.PrimaryHDU(data=std_diff)
    new_hdu.writeto(savedir + 'std_diff_' + linename + '_2m1.fits', overwrite=True)
    del new_hdu

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)