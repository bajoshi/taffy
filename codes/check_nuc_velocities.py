from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import line_fluxes as lf

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + '/Desktop/ipac/taffy/'

if __name__ == '__main__':
    
    # constants:
    lightspeed = 299792.458  # km/s

    line_arr_northnuc, wav_arr_northnuc, line_flux_northnuc = lf.get_line_array('northnuc_sum_R_line', 'Halpha', 50, 50)
    line_arr_southnuc, wav_arr_southnuc, line_flux_southnuc = lf.get_line_array('southnuc_sum_R_line', 'Halpha', 50, 50)

    # make figure and grid to plot and set up axes
    gs = gridspec.GridSpec(1,3, height_ratios=[1], width_ratios=[1,0.2,1])
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure()
    ax_northnuc = fig.add_subplot(gs[0,0])
    ax_southnuc = fig.add_subplot(gs[0,2])

    halpha_air_wav = 6562.801
    helio_vel_arr_northnuc = ((wav_arr_northnuc - halpha_air_wav) / halpha_air_wav) * lightspeed
    helio_vel_arr_southnuc = ((wav_arr_southnuc - halpha_air_wav) / halpha_air_wav) * lightspeed

    ax_northnuc.plot(helio_vel_arr_northnuc, line_arr_northnuc, color='r')
    ax_southnuc.plot(helio_vel_arr_southnuc, line_arr_southnuc, color='r')

    # find peak and draw a vertical line and print out velocity at peak value
    halpha_north_peak_idx = np.argmax(line_arr_northnuc)
    halpha_south_peak_idx = np.argmax(line_arr_southnuc)

    print "North galaxy velocity at peak H-alpha", helio_vel_arr_northnuc[halpha_north_peak_idx]
    print "South galaxy velocity at peak H-alpha", helio_vel_arr_southnuc[halpha_south_peak_idx]

    ax_northnuc.axvline(x=helio_vel_arr_northnuc[halpha_north_peak_idx], ls='--')
    ax_southnuc.axvline(x=helio_vel_arr_southnuc[halpha_south_peak_idx], ls='--')

    plt.show()

    sys.exit(0)