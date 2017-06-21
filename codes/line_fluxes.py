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
taffy_dir = home + '/Desktop/ipac/taffy/'
hcg95c_dir = home + '/Desktop/ipac/hcg95c/'

def get_line_array(region_name, line_name, line_width_low, line_width_up):

    #  Common lines:
    labels = np.array(['OII3726', 'OII3729', 'NeIII3869', 'Hepsilon', 'Hdelta', 'Hgamma', 'OIII4363', 'Hbeta',\
     'OIII4959', 'OIII5007', 'OI6300', 'OI6364', 'NII6548', 'Halpha', 'NII6583', 'SII6716', 'SII6731'])

    wavelengths = np.array([3726.032, 3728.815, 3869.060, 3970.072, 4101.734, 4340.464, 4363.210, 4861.325,\
     4958.911, 5006.843, 6300.304, 6363.7, 6548.040, 6562.800, 6583.460, 6716.440, 6730.810])

    # select line(s) to plot
    redshift = 0.0145
    line_idx = np.where(labels == line_name)[0]
    line_wav = wavelengths[line_idx] * (1+redshift)
    print line_wav
    lam_low = line_wav - line_width_low
    lam_up = line_wav + line_width_up

    # lower limits for [NII]6583: 6672, 6668.4, 6667, 6666, 6668, 6664.7, 6670
    # upper limit for [NII]6583: 6700

    # lower limit for [NII]6548: 6625
    # upper limits for [NII]6548: 6652, 6648, 6646, 6647, 6649, 6647.7, 6650.5

    # lower limit for [SII]6716: 6785
    # upper limits for [SII]6716: 6824.5, 6820.8, 6818.5, 6818, 6822, 6820.6, 6822.9 # these are also lower limits for [SII]6731
    # upper limit for [SII]6731: 6855

    specfile = np.genfromtxt(taffy_dir + 'taffy_cube_summed_fits/' + region_name + '.dat', dtype=None, names=['wav', 'flux'])

    lam = specfile['wav']
    line_flux = specfile['flux']

    lam_low_idx = np.argmin(abs(lam - lam_low))
    lam_up_idx = np.argmin(abs(lam - lam_up))

    line_arr = line_flux[lam_low_idx:lam_up_idx+1]
    wav_arr = np.linspace(lam_low, lam_up, len(line_arr))

    """
    if you want the flux in the line (in the same units as the line) then 
    integrate the area under the line in the observed data
    in this case make sure that hte given with **ONLY** covers the required line
    BE CAREFUL with the width!!
    the line flux varies quite a bit (and not intuitively)
    depending on where the end points of the integration are 
    """
    line_flux_intg = simps(y=line_arr, x=wav_arr)

    return line_arr, wav_arr, line_flux_intg

def get_herschel_spectra(spectype='wav'):

    bridge_b1_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b1_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    bridge_b1_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b1_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    bridge_b2_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b2_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    bridge_b2_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b2_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    bridge_b3_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b3_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    bridge_b3_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b3_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    bridge_b4_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b4_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    bridge_b4_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b4_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    taffy_n_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/taffy_n_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    taffy_n_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/taffy_n_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    taffy_s_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/taffy_s_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    taffy_s_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/taffy_s_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    exgalac_hii_cp = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/exgalac_hii_sum_c+.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])
    exgalac_hii_oi = np.genfromtxt(home + '/Desktop/ipac/taffy/herschel_spectra_from_phil/exgalac_hii_sum_OI.dat-' + spectype + '.txt',\
     dtype=None, names=['wav','flux'])

    return bridge_b1_cp, bridge_b1_oi, bridge_b2_cp, bridge_b2_oi, bridge_b3_cp, bridge_b3_oi, bridge_b4_cp, bridge_b4_oi,\
    taffy_n_cp, taffy_n_oi, taffy_s_cp, taffy_s_oi, exgalac_hii_cp, exgalac_hii_oi

def get_herschel_flux(specfile, lam_intg_low, lam_intg_high):

    lam = specfile['wav']
    flux = specfile['flux']

    low_idx = np.argmin(abs(lam - lam_intg_low))
    high_idx = np.argmin(abs(lam - lam_intg_high))

    flux_intg = flux[low_idx:high_idx+1]
    wav_intg = lam[low_idx:high_idx+1]

    line_flux = simps(y=flux_intg, x=wav_intg)

    return line_flux

def plot_herschel_spectrum(specfile):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(specfile['wav'], specfile['flux'], color='k')
    ax.axhline(y=0, color='b', linestyle='--')

    plt.show()
    plt.clf()
    plt.cla()
    plt.close()

    return None

def plot_all_herschel(allspec, allspecnames):

    for i in range(len(allspec)):    
        print allspecnames[i]
        plot_herschel_spectrum(allspec[i])

    return None

if __name__ == '__main__':
    """
    This code overplots the optical lines and FIR lines and compares the total 
    fluxes within the lines.
    """

    # use these lines if you want the line flux in just a single region
    # if it is a line that LZIFU provides then just overplot the region in ds9
    # over the total line map and get the region properties by clicking on analysis 
    # in region information.
    """
    single_region_name = 'snuc_2017new_R_line'
    single_line_name = 'Halpha'
    line_arr, wav_arr, line_flux = get_line_array(single_region_name, single_line_name, 8.16, 11.54)
    print "Flux in", single_line_name, "{:.3f}".format(line_flux * 1e-18 * 1e-3 / 1e-16), "for", single_region_name

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wav_arr, line_arr)
    ax.axhline(y=0, ls='--', color='k')
    plt.show()

    sys.exit(0)
    """

    # get line and wavelength array to plot and also the flux in the line

    alllines = ['Hbeta', 'OIII4959', 'OIII5007', 'OI6300', 'OI6364 ', \
    'NII6548', 'Halpha', 'NII6583', 'SII6716', 'SII6731']

    all_blue_regions = ['taffy_n_sum_B_line', 'taffy_s_sum_B_line', 'bridge_b1_sum_B_line', 'bridge_b2_sum_B_line', \
    'bridge_b3_sum_B_line', 'bridge_b4_sum_B_line', 'exgalac_hii_sum_B_line']

    all_red_regions = ['taffy_n_sum_R_line', 'taffy_s_sum_R_line', 'bridge_b1_sum_R_line', 'bridge_b2_sum_R_line', \
    'bridge_b3_sum_R_line', 'bridge_b4_sum_R_line', 'exgalac_hii_sum_R_line']

    line_name = 'SII6716'  # the current line who's integrated flux you want

    for region_name in all_red_regions:
        # make sure that this says all_blue or all_red regions depending on which channel the line actually is in

        line_arr, wav_arr, line_flux = get_line_array(region_name, line_name, 50, 50)
        print line_name, '& $', "{:.3f}".format(line_flux * 1e-18 * 1e-3 / 1e-16), '$ &'
        # the 1e-18 converts it to erg s^-1 cm^-2 and the 1e-3 converts that to W m^-2
        # the 1e-16 further will write it in units of [W/m^2]x1e-16

        # plot line
        # check if the limits for integrating are actually correct
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(wav_arr, line_arr)
        plt.show()

    sys.exit(0)

    # get herschel spectra and their fluxes
    bridge_b1_cp, bridge_b1_oi, bridge_b2_cp, bridge_b2_oi, bridge_b3_cp, bridge_b3_oi, bridge_b4_cp, bridge_b4_oi,\
    taffy_n_cp, taffy_n_oi, taffy_s_cp, taffy_s_oi, exgalac_hii_cp, exgalac_hii_oi = get_herschel_spectra()

    allspec = [bridge_b1_cp, bridge_b1_oi, bridge_b2_cp, bridge_b2_oi, bridge_b3_cp, bridge_b3_oi, bridge_b4_cp, bridge_b4_oi,\
     taffy_n_cp, taffy_n_oi, taffy_s_cp, taffy_s_oi, exgalac_hii_cp, exgalac_hii_oi]

    allspecnames = ['bridge_b1_cp', 'bridge_b1_oi', 'bridge_b2_cp', 'bridge_b2_oi', 'bridge_b3_cp', 'bridge_b3_oi',\
     'bridge_b4_cp', 'bridge_b4_oi', 'taffy_n_cp', 'taffy_n_oi', 'taffy_s_cp', 'taffy_s_oi', 'exgalac_hii_cp', 'exgalac_hii_oi']

    #plot_all_herschel(allspec, allspecnames)

    # eyeballed limits for integration
    # we only want to integrate where the line is actually there
    print 'bridge_b1_cp_flux [W/m^2]', get_herschel_flux(bridge_b1_cp, 159.58, 160.5) * 1e-26 * 1.2e10
    print 'bridge_b2_cp_flux [W/m^2]', get_herschel_flux(bridge_b2_cp, 159.59, 160.4) * 1e-26 * 1.2e10
    print 'bridge_b3_cp_flux [W/m^2]', get_herschel_flux(bridge_b3_cp, 159.64, 160.28) * 1e-26 * 1.2e10
    print 'bridge_b4_cp_flux [W/m^2]', get_herschel_flux(bridge_b4_cp, 159.92, 160.34) * 1e-26 * 1.2e10
    print 'taffy_n_cp_flux [W/m^2]', get_herschel_flux(taffy_n_cp, 159.66, 160.5) * 1e-26 * 1.2e10
    print 'taffy_s_cp_flux [W/m^2]', get_herschel_flux(taffy_s_cp, 159.57, 160.41) * 1e-26 * 1.2e10
    print 'exgalac_hii_cp_flux [W/m^2]', get_herschel_flux(exgalac_hii_cp, 159.66, 160.42) * 1e-26 * 1.2e10

    # by eye I did not see any OI in b1 and b4
    print 'bridge_b1_oi_flux [W/m^2]', get_herschel_flux(bridge_b1_oi, 64.03, 64.215) * 1e-26 * 7.515e10
    print 'bridge_b2_oi_flux [W/m^2]', get_herschel_flux(bridge_b2_oi, 64.035, 64.15) * 1e-26 * 7.515e10
    print 'bridge_b3_oi_flux [W/m^2]', get_herschel_flux(bridge_b3_oi, 64.06, 64.167) * 1e-26 * 7.515e10
    print 'bridge_b4_oi_flux [W/m^2]', get_herschel_flux(bridge_b4_oi, 64.03, 64.215) * 1e-26 * 7.515e10
    print 'taffy_n_oi_flux [W/m^2]', get_herschel_flux(taffy_n_oi, 64.03, 64.215) * 1e-26 * 7.515e10
    print 'taffy_s_oi_flux [W/m^2]', get_herschel_flux(taffy_s_oi, 64.025, 64.172) * 1e-26 * 7.515e10
    print 'exgalac_hii_oi_flux [W/m^2]', get_herschel_flux(exgalac_hii_oi, 64.0467, 64.1839) * 1e-26 * 7.515e10

    sys.exit(0)