from __future__ import division

import numpy as np
from astropy.io import fits

import sys
import os
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea
import matplotlib.transforms as mtransforms

home = os.getenv('HOME')  # does not have a trailing slash
desktop = home + '/Desktop/'
stacking_analysis_dir = home + "/Desktop/FIGS/stacking-analysis-pears/"

taffy_products = '/Users/baj/Desktop/ipac/taffy_lzifu/products/'
taffy_data = '/Users/baj/Desktop/ipac/taffy_lzifu/data/'
ipac_taffy_dir = home + '/Desktop/ipac/taffy/'
ipac_taffy_figdir = ipac_taffy_dir + 'figures/'
taffy_extdir = '/Users/baj/Desktop/ipac/taffy_lzifu/'

sys.path.append(stacking_analysis_dir + 'codes/')
import fast_chi2_jackknife as fcj

if __name__ == '__main__':
	
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in products fits file
    h = fits.open(taffy_products + 'Taffy_1_comp.fits')

    #total_ext = fcj.get_total_extensions(h)
    #bpt.print_extnames(h, total_ext)

    # read in observed data
    obs_b = fits.open(taffy_data + 'Taffy_B.fits')
    obs_r = fits.open(taffy_data + 'Taffy_R.fits')

    obs_noise_b = obs_b[1].data
    obs_data_b = obs_b[0].data
    obs_noise_r = obs_r[1].data
    obs_data_r = obs_r[0].data

    # read in required data for fit from red and blue channels
    b_cont = h['B_CONTINUUM'].data
    b_resid = h['B_RESIDFIT'].data
    b_line = h['B_LINE'].data
    b_mask = h['B_CONT_MASK'].data

    r_cont = h['R_CONTINUUM'].data
    r_resid = h['R_RESIDFIT'].data
    r_line = h['R_LINE'].data
    r_mask = h['R_CONT_MASK'].data

    # find ranges of fit
    #print "y and x array coords of fitted spaxels for blue and red cube respectively; continuum -- ",\
    # np.unique(np.where(np.isfinite(b_cont))[1]), np.unique(np.where(np.isfinite(b_cont))[2])
    #print "y and x array coords of fitted spaxels for blue and red cube respectively; continuum -- ",\
    # np.unique(np.where(np.isfinite(r_cont))[1]), np.unique(np.where(np.isfinite(r_cont))[2])
    #print "y and x array coords of fitted spaxels for blue and red cube respectively; line -- ",\
    # np.unique(np.where(np.isfinite(b_line))[1]), np.unique(np.where(np.isfinite(b_line))[2])
    #print "y and x array coords of fitted spaxels for blue and red cube respectively; line -- ",\
    # np.unique(np.where(np.isfinite(r_line))[1]), np.unique(np.where(np.isfinite(r_line))[2])

    # create lambda arrays
    b_cont_hdr = h['B_CONTINUUM'].header
    tot_spec_elem = b_cont_hdr['NAXIS3']
    spec_start = b_cont_hdr['CRVAL3']
    spec_delta = b_cont_hdr['CDELT3']
    lam_end = spec_start + tot_spec_elem*spec_delta - spec_delta
    lam_b = np.linspace(spec_start, lam_end, tot_spec_elem)

    r_cont_hdr = h['R_CONTINUUM'].header
    tot_spec_elem = r_cont_hdr['NAXIS3']
    spec_start = r_cont_hdr['CRVAL3']
    spec_delta = r_cont_hdr['CDELT3']
    lam_end = spec_start + tot_spec_elem*spec_delta - spec_delta
    lam_r = np.linspace(spec_start, lam_end, tot_spec_elem)

    # plot the spectrum from a pixel
    # these pixel coordinates correspond to the ones shown in ds9
    # teh array coordinates are tranformations to the numpy array coords
    # if you only fit a small region by editing hte i and j ranges in 
    # lzifu_loop_spaxel.pro it includes both limits of the range.
    # example: 
    # for i=9,14 do begin
    #     for j=35,40 do begin
    # will include both x pixel values from 9 (+1) to 14 (+1) (ds9 values) including both numbers
    # and similarly for y.
    pix_x = 11
    pix_y = 39
    arr_x = pix_y - 1
    arr_y = pix_x - 1

    print 'Current pixel (ds9 coords)', pix_x, pix_y

    # actual plot using gridspec 
    gs = gridspec.GridSpec(20,20)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.00, hspace=3.0)

    fig = plt.figure()
    ax1 = fig.add_subplot(gs[:15,:10])
    ax2 = fig.add_subplot(gs[15:,:10])
    ax3 = fig.add_subplot(gs[:15,10:])
    ax4 = fig.add_subplot(gs[15:,10:])

    # blue subplots
    flam_obs_b = obs_data_b[:,arr_x,arr_y]
    ferr_obs_b = np.sqrt(obs_noise_b[:,arr_x,arr_y])
    cont_fit_b = b_cont[:,arr_x,arr_y]
    line_fit_b = b_line[:,arr_x,arr_y]
    mask_fit_b = b_mask[:,arr_x,arr_y]
    residfit_b = (flam_obs_b - (cont_fit_b + line_fit_b)) / ferr_obs_b

    lam_b_masked = np.ma.array(lam_b, mask=mask_fit_b)

    if np.any(np.isfinite(line_fit_b+cont_fit_b)):
        ax1.plot(lam_b, line_fit_b+cont_fit_b, color='b', lw=1.3, zorder=10)
    else:
        ax1.plot(lam_b, cont_fit_b, color='b', zorder=10)

    # to show line and continuum fits by themselves
    #ax1.plot(lam_b, cont_fit_b, color='g', zorder=10)
    #ax1.plot(lam_b, line_fit_b, color='g', zorder=10)
    
    ax1.plot(lam_b, flam_obs_b, '-', color='gray')
    #ax1.plot(lam_b, ferr_obs_b, '--', color='lightgray', linewidth=1, zorder=5)
    #ax1.axhline(y=ferr_obs_b[100], linestyle='--', linewidth=3, color='lightgray', zorder=5)
    #ax1.fill_between(lam_b, flam_obs_b + ferr_obs_b, flam_obs_b - ferr_obs_b, color='lightgray')

    print "valid wavlength indices for B, continuum:", np.unique(np.where(np.isfinite(cont_fit_b)))
    print "valid wavlength indices for B, line:", np.unique(np.where(np.isfinite(line_fit_b)))
    
    trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    #ax1.fill_between(lam_b_masked, 0, 1, facecolor='b', alpha=0.3, transform=trans)

    ax1.set_xlim(4800, 5300)  # old lims (4800, 5300)
    #ax1.set_ylim(35,90)

    ax1.get_xaxis().set_ticklabels([])
    #ax1.get_xaxis().set_ticklabels(['4600', '4700', '4800', '4900', '5000', '5100', '5200', '5300', ''],\
    # fontsize=10, rotation=35)
    ax1.tick_params(axis='both', which='major', labelsize=10)

    ax1.set_ylabel(r'$\mathrm{f_\lambda\ [erg\,s^{-1}\,cm^{-2}\,\AA^{-1}] / 1\times10^{-18}}$')
    ax1.yaxis.set_label_coords(-0.1, 0.5)

    ax1.minorticks_on()
    ax1.tick_params('both', width=1, length=3, which='minor')
    ax1.tick_params('both', width=1, length=4.7, which='major')
    ax1.grid(True)

    # second blue subplot
    ax2.plot(lam_b, residfit_b, color='gray')

    ax2.set_xlim(4800, 5300)
    ax2.get_xaxis().set_ticklabels(['4800', '4900', '5000', '5100', '5200', '5300', ''],\
     fontsize=10, rotation=35)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    ax2.set_ylabel(r'$\mathrm{(f^{data}_\lambda - f^{model}_\lambda)/\sigma}$')

    ax2.minorticks_on()
    ax2.tick_params('both', width=1, length=3, which='minor')
    ax2.tick_params('both', width=1, length=4.7, which='major')
    ax2.grid(True)

    # red subplots
    flam_obs_r = obs_data_r[:,arr_x,arr_y]
    ferr_obs_r = obs_noise_r[:,arr_x,arr_y]
    cont_fit_r = r_cont[:,arr_x,arr_y]
    line_fit_r = r_line[:,arr_x,arr_y]
    mask_fit_r = r_mask[:,arr_x,arr_y]
    residfit_r = (flam_obs_r - (cont_fit_r + line_fit_r)) / ferr_obs_r

    lam_r_masked = np.ma.array(lam_r, mask=mask_fit_r)

    print "valid wavlength indices for R, continuum:", np.unique(np.where(np.isfinite(cont_fit_r)))
    print "valid wavlength indices for R, line:", np.unique(np.where(np.isfinite(line_fit_r)))

    if np.any(np.isfinite(line_fit_r+cont_fit_r)):
        ax3.plot(lam_r, line_fit_r+cont_fit_r, color='r', lw=1.3, zorder=10)
    else:
        ax3.plot(lam_r, cont_fit_r, color='r', zorder=10)

    #ax3.plot(lam_r, cont_fit_r, color='g', zorder=10)
    #ax3.plot(lam_r, line_fit_r, color='g', zorder=10)

    ax3.plot(lam_r, flam_obs_r, '-', color='gray')
    #ax3.plot(lam_r, ferr_obs_r, '--', color='lightgray', linewidth=1, zorder=5)
    #ax3.axhline(y=ferr_obs_r[100], linestyle='--', linewidth=3, color='lightgray', zorder=5)
    #ax3.fill_between(lam_r, flam_obs_r + ferr_obs_r, flam_obs_r - ferr_obs_r, color='lightgray')

    trans = mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
    #ax3.fill_between(lam_r_masked, 0, 1, facecolor='r', alpha=0.3, transform=trans)

    ax3.set_xlim(6370, 6700)  # old lims (6370, 6700)
    #ax3.set_ylim(0,300)

    ax3.get_xaxis().set_ticklabels([])
    #ax3.get_xaxis().set_ticklabels(['6100', '6200', '6300', '6400', '6500', '6600', '6700', '6800', '6900'],\
    # fontsize=10, rotation=35)
    ax3.tick_params(axis='both', which='major', labelsize=10)

    ax3.minorticks_on()
    ax3.tick_params('both', width=1, length=3, which='minor')
    ax3.tick_params('both', width=1, length=4.7, which='major')
    ax3.grid(True)

    # second red subplot
    ax4.plot(lam_r, residfit_r, color='gray')

    ax4.set_xlabel(r'$\mathrm{Wavelength\ [\AA]}$')
    ax4.xaxis.set_label_coords(0.00, -0.3)

    ax4.set_xlim(6370, 6700)
    ax4.get_xaxis().set_ticklabels(['', '6400', '6450', '6500', '6500', '6600', '6650', '6700'],\
     fontsize=10, rotation=35)
    ax4.tick_params(axis='both', which='major', labelsize=10)

    ax4.minorticks_on()
    ax4.tick_params('both', width=1, length=3, which='minor')
    ax4.tick_params('both', width=1, length=4.7, which='major')
    ax4.grid(True)

    plt.show()

    #fig.savefig(ipac_taffy_figdir + 'lzifu_fit_pix_ds9xy_' + str(pix_x) + '_' + str(pix_y) + '_gray.eps',\
    # dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # overplot two pixels
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    flam_obs_b1 = obs_data_b[:,arr_x,arr_y]
    flam_obs_b2 = obs_data_b[:,arr_x+1,arr_y]

    ax.plot(lam_b, flam_obs_b1, '-', color='gray')
    ax.plot(lam_b, flam_obs_b2-50, '-', color='r')

    ax.set_xlim(4600, 5400)

    plt.show()
    """

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)