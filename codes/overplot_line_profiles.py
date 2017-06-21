from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

import line_fluxes as lf

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_dir = home + '/Desktop/ipac/taffy/'

if __name__ == '__main__':

    # constants:
    lightspeed = 299792.458  # km/s
    
    #  Common lines:
    labels = np.array(['OII3726', 'OII3729', 'NeIII3869', 'Hepsilon', 'Hdelta', 'Hgamma', 'OIII4363', 'Hbeta',\
     'OIII4959', 'OIII5007', 'OI6300', 'OI6364', 'NII6548', 'Halpha', 'NII6583', 'SII6716', 'SII6731'])

    wavelengths = np.array([3726.032, 3728.815, 3869.060, 3970.072, 4101.734, 4340.464, 4363.210, 4861.325,\
     4958.911, 5006.843, 6300.304, 6363.7, 6548.040, 6562.800, 6583.460, 6716.440, 6730.810])

    all_blue_regions = ['taffy_n_sum_B_line', 'taffy_s_sum_B_line', 'bridge_b1_sum_B_line', 'bridge_b2_sum_B_line', \
    'bridge_b3_sum_B_line', 'bridge_b4_sum_B_line', 'exgalac_hii_sum_B_line']

    all_red_regions = ['taffy_n_sum_R_line', 'taffy_s_sum_R_line', 'bridge_b1_sum_R_line', 'bridge_b2_sum_R_line', \
    'bridge_b3_sum_R_line', 'bridge_b4_sum_R_line', 'exgalac_hii_sum_R_line']

    bridge_b1_cp, bridge_b1_oi, bridge_b2_cp, bridge_b2_oi, bridge_b3_cp, bridge_b3_oi, bridge_b4_cp, bridge_b4_oi,\
    taffy_n_cp, taffy_n_oi, taffy_s_cp, taffy_s_oi, exgalac_hii_cp, exgalac_hii_oi = lf.get_herschel_spectra(spectype='helio')

    # make figure and grid to plot and set up axes
    gs1 = gridspec.GridSpec(3,5, height_ratios=[1,0.05,1], width_ratios=[1,0.05,1,0.05,1])
    gs1.update(left=0.1, right=0.9, bottom=0.35, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure()

    ax_north = fig.add_subplot(gs1[0,0])
    ax_south = fig.add_subplot(gs1[0,2])

    ax_b1 = fig.add_subplot(gs1[0,4])
    ax_b2 = fig.add_subplot(gs1[2,0])
    ax_b3 = fig.add_subplot(gs1[2,2])
    ax_b4 = fig.add_subplot(gs1[2,4])

    gs2 = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
    gs2.update(left=0.37, right=0.63, bottom=0.06, top=0.32, wspace=0.02, hspace=0.02)

    ax_xgal = fig.add_subplot(gs2[0,0])

    all_axes = [ax_north, ax_south, ax_b1, ax_b2, ax_b3, ax_b4, ax_xgal]
    all_herschel_cp_spectra = [taffy_n_cp, taffy_s_cp, bridge_b1_cp, bridge_b2_cp, bridge_b3_cp, bridge_b4_cp, exgalac_hii_cp]

    # actual plotting
    line_name = 'OI6300'

    count = 0
    for region_name in all_red_regions:

        line_arr, wav_arr, line_flux = lf.get_line_array(region_name, line_name, 50, 50)

        halpha_air_wav = 6562.801
        hbeta_air_wav = 4861.363
        oi_air_wav = 6300.30

        if line_name == 'Halpha':
            helio_vel_arr = ((wav_arr - halpha_air_wav) / halpha_air_wav) * lightspeed
        elif line_name == 'Hbeta':
            helio_vel_arr = ((wav_arr - hbeta_air_wav) / hbeta_air_wav) * lightspeed
        elif line_name == 'OIII5007':
            helio_vel_arr = ((wav_arr - 5008.239) / 5008.239) * lightspeed
        elif line_name == 'OIII4959':
            helio_vel_arr = ((wav_arr - 4960.295) / 4960.295) * lightspeed
        elif line_name == 'OI6300':
            helio_vel_arr = ((wav_arr - oi_air_wav) / oi_air_wav) * lightspeed

        all_axes[count].plot(helio_vel_arr, line_arr, color='r')

        yticklabels = all_axes[count].get_yticks().tolist()
        xticklabels = all_axes[count].get_xticks().tolist()
        if count == 3: xticklabels[-1] = ''
        if count == 5: xticklabels[0] = ''

        #if count in [0,1]:
        #    all_axes[count].set_ylim(-2000,16000)
        #    all_axes[count].set_yticklabels(['-0.2', '0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6'],\
        #     size=8, rotation=30)
        #    all_axes[count].tick_params('both', which='major', pad=1)

        #if count == 2:
        #    all_axes[count].set_ylim(-100,900)

        #if count == 3:
        #    all_axes[count].set_ylim(-200,1500)

        #if count == 4:
        #    all_axes[count].set_ylim(-200,1200)

        #if count == 5:
        #    all_axes[count].set_ylim(-50,300)

        #if count == 6:
        #    all_axes[count].set_ylim(-500,3500)

        all_axes[count].set_xticklabels(xticklabels, size=8, rotation=30)
        all_axes[count].set_yticklabels(yticklabels, size=8, rotation=30)

        if count in [0,1,2,4]:
            all_axes[count].set_xticklabels([])

        # ----------- twin axes plotting Herschel C+ data ----------- #
        ax_cp = all_axes[count].twinx()

        ax_cp.plot(all_herschel_cp_spectra[count]['wav'], all_herschel_cp_spectra[count]['flux'], color='k')
        #ax_cp.axhline(y=0, color='b', linestyle='--')

        # set limits before setting tick labels
        if count == 0:
            ax_cp.set_ylim(-10, 60)

        if count == 1:
            ax_cp.set_ylim(-5, 35)

        if count == 2:
            ax_cp.set_ylim(-0.5,2.5)

        ax_cp.set_yticklabels(ax_cp.get_yticks().tolist(), size=8, rotation=30)

        count += 1

    fig.savefig(taffy_dir + 'figures/' + line_name + '_cplus_profile.eps', dpi=300, bbox_inches='tight')

    plt.show()
    sys.exit(0)