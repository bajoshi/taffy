from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.modeling import models, fitting

import sys
import os

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + '/Desktop/ipac/taffy/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

speed_of_light = 299792.458  # km/s
halpha_air_wav = 6562.80  # Angstroms

if __name__ == '__main__':

    # read in LZIFU line fitting results
    total_spec = np.genfromtxt(taffy_extdir + 'taffys_nuc_avg_R_line.dat', dtype=None, names=['channel', 'lineflux'])
    comp1_spec = np.genfromtxt(taffy_extdir + 'taffys_nuc_avg_R_line_comp1.dat', dtype=None, names=['channel', 'lineflux'])
    comp2_spec = np.genfromtxt(taffy_extdir + 'taffys_nuc_avg_R_line_comp2.dat', dtype=None, names=['channel', 'lineflux'])

    # read in observed data
    obs_data_r_total = np.genfromtxt(taffy_extdir + 'taffys_nuc_avg_R_data.dat', dtype=None, names=['channel', 'lineflux'])

    # create wavelength array
    # I read these data from the corresponding headers
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

    wav_arr = red_wav_arr
    line_comp1 = comp1_spec['lineflux']
    line_comp2 = comp2_spec['lineflux']
    line_total = total_spec['lineflux']
    obs_data = obs_data_r_total['lineflux']
    data_res = red_res
    delt = delt_r

    # Start fitting
    # find the center of the biggest peak and call that the line idx
    line_wav = line_air_wav*(1+redshift)
    line_idx = np.argmin(abs(wav_arr - line_wav))
    max_val_idx = np.argmax(obs_data[line_idx-45: line_idx+45])
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

    line_y_arr_comp1 = line_comp1[line_idx-linepad_left:line_idx+linepad_right]
    line_x_arr_comp1 = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_comp1))

    line_y_arr_comp2 = line_comp2[line_idx-linepad_left:line_idx+linepad_right]
    line_x_arr_comp2 = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_comp2))

    # fitting
    gauss_init_lowcomp = models.Gaussian1D(amplitude=5.0, mean=line_idx-10, stddev=5.0)
    gauss_init_highcomp = models.Gaussian1D(amplitude=5.0, mean=line_idx+10, stddev=5.0)
    fit_gauss = fitting.LevMarLSQFitter()

    g1 = fit_gauss(gauss_init_lowcomp, line_x_arr_comp1, line_y_arr_comp1)
    g2 = fit_gauss(gauss_init_highcomp, line_x_arr_comp2, line_y_arr_comp2)

    # save lzifu total fit to array for plotting
    line_y_arr_total = line_total[line_idx-linepad_left:line_idx+linepad_right]

    # also fit raw data by a single gaussian
    line_y_arr_data = obs_data[line_idx-linepad_left:line_idx+linepad_right]
    line_x_arr_data = np.linspace(line_idx-linepad_left, line_idx+linepad_right, len(line_y_arr_data))

    # find pseudo continuum and subtract
    pseudo_cont_arr_left = obs_data[line_idx-10-(linepad_left-10):line_idx-10]
    pseudo_cont_arr_right = obs_data[line_idx+10:line_idx+10+(linepad_right-10)]
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
        amp_comp1 = g2.parameters[0]
        vel_comp1 = g2.parameters[1]
        std_comp1 = abs(g2.parameters[2])

        amp_comp2 = g1.parameters[0]
        vel_comp2 = g1.parameters[1]
        std_comp2 = abs(g1.parameters[2])

    elif phys_vel_comp2 >= phys_vel_comp1:
        amp_comp1 = g1.parameters[0]
        vel_comp1 = g1.parameters[1]
        std_comp1 = abs(g1.parameters[2])

        amp_comp2 = g2.parameters[0]
        vel_comp2 = g2.parameters[1]
        std_comp2 = abs(g2.parameters[2])

    amp_onecomp = g.parameters[0]
    vel_onecomp = g.parameters[1]
    std_onecomp = abs(g.parameters[2])

    # print useful info
    print "mean vel 1 and 2 [km/s]", phys_vel_comp1, phys_vel_comp2
    print "mean diff", format((((vel_comp2 - vel_comp1) * 0.3) / line_air_wav) * speed_of_light, '.2f'), "km/s"
    print "std devs comp2 and comp1", format(std_comp2, '.2f'), format(std_comp1, '.2f')
    print "continuum level, amp and vdisp for onecomp fit", \
    format(cont_level, '.2f'), format(amp_onecomp, '.2f'), format(std_onecomp, '.2f')

    # plot to check
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel('Channel', fontsize=12)
    ax.set_ylabel(r'$\mathrm{Emission\ line\ flux}\ [10^{-18}\, \mathrm{erg\, s^{-1}\, cm^{-2}\, \AA^{-1}}]$', fontsize=12)

    # plot first comp fit
    #ax.plot(line_x_arr_comp1, line_y_arr_comp1, '.', color='b')
    ax.plot(line_x_arr_comp1, g1(line_x_arr_comp1), ls='--', color='b', lw=2, label='Low vel component fit')

    # plot second comp fit
    #ax.plot(line_x_arr_comp2, line_y_arr_comp2, '.', color='r')
    ax.plot(line_x_arr_comp2, g2(line_x_arr_comp2), ls='--', color='r', lw=2, label='High vel component fit')

    # plot lzifu total fit
    ax.plot(line_x_arr_comp1, line_y_arr_total, color='k', label='LZIFU fit (two comp)')

    # also show the raw data and *MY* single gaussian fit to the raw data
    ax.plot(line_x_arr_data, line_y_arr_data, color='gray', label='Raw data')
    ax.plot(line_x_arr_data, g(line_x_arr_data), ls='--', color='g', lw=2, label='Single Gaussian fit')

    # Text for line FWHM
    f = FontProperties()
    f.set_weight('bold')

    fwhm_low = (std_comp1 * 2.354 * 0.3 / halpha_air_wav) * speed_of_light
    fwhm_high = (std_comp2 * 2.354 * 0.3 / halpha_air_wav) * speed_of_light

    ax.text(0.02, 0.95, r'$\mathrm{FWHM_{Low\ vel} = }$' + str(int(fwhm_low)) + r'$\mathrm{\, km\, s^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)
    ax.text(0.02, 0.9, r'$\mathrm{FWHM_{High\ vel} = }$' + str(int(fwhm_high)) + r'$\mathrm{\, km\, s^{-1}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', fontproperties=f, size=13)

    ax.legend()
    ax.minorticks_on()

    fig.savefig(taffy_extdir + 'figures_stitched_cube/gaussfit_example_snuc.png', dpi=300, bbox_inches='tight')

    sys.exit(0)