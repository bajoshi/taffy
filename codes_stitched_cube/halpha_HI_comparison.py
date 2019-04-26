from __future__ import division

import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Gaussian2DKernel
import Polygon as pg

import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import FontProperties

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffy_products = home + '/Desktop/ipac/taffy_lzifu/products_work/'
taffy_data = home + '/Desktop/ipac/taffy_lzifu/data/'
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'
taffydir = home + '/Desktop/ipac/taffy/'
gaussfits_dir = home + '/Desktop/ipac/taffy_lzifu/baj_gauss_fits_to_lzifu_linefits/'

sys.path.append(taffydir + 'codes/')
import vel_channel_map as vcm
import stitch_map as sm

speed_of_light = 299792.458  # km/s
redshift = 0.0145

def extract_spec(reg_pn, cube):
    """
    This function will extract and return the spectrum 
    corresponding to the given region from the given cube.
    Arguments:
    --- reg: ds9 region polygon in image coordinates.
    i.e., convert the region from the ds9 file to a polygon object.

    --- cube: data cube from which spectrum is to be extracted.
    """

    # Use polygon to generate mask which can then 
    # be used on cube which is anumpy array 
    # Code based on region mask making code in bpt_plots.py
    # that one is optimzed for the Taffy IFU data 
    # Here, I need it to work for both the optical and the HI cubes.
    hi_full_shape = (31, 128, 128)
    hi_cube_shape2d = (128,128)
    #ifu_cube_shape = (58,58)
    regionmask = np.zeros(hi_cube_shape2d, dtype=np.int)

    bbox = reg_pn.boundingBox()
    # Slightly padding all sides
    bbox_xmin = int(bbox[0] - 0.1)
    bbox_xmax = int(bbox[1] + 0.1)
    bbox_ymin = int(bbox[2] - 0.1)
    bbox_ymax = int(bbox[3] + 0.1)

    # Test python (row,col) to ds9 (x,y) conversion
    # Check spaxel spectrum in ds9 and the one plotted here
    # ds9 coords to test
    """
    ds9_x = 57
    ds9_y = 54
    py_row = ds9_y - 1
    py_col = ds9_x - 1

    print cube.shape
    print "Value at DS9 pixel (x,y):", ds9_x, ds9_y, "is", np.nanmean(cube[:, :, py_row, py_col])
    print cube[:, :, py_row, py_col].reshape((31))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(31), cube[:, :, py_row, py_col].reshape((31)), lw=0.5)
    ax.grid(True, lw=0.1)
    plt.show()
    sys.exit()
    """

    # loop over all possible pixels inside the bounding box
    # i.e., you need to loop over all possible pixel centers
    # if the pixel center is within the polygon then it gets
    # selected otherwise skip.
    for i in range(bbox_xmin, bbox_xmax+1):
        for j in range(bbox_ymin, bbox_ymax+1):

            # convert from ds9 coords to array coords
            arr_x = int(j) - 1
            arr_y = int(i) - 1

            regionmask[arr_x,arr_y] = reg_pn.isInside(i,j)

    """
    # plot to test
    # all the -1.0 here are to show everything in python frame
    # when you compare to ds9 the x and y are now in the correct
    # direction but you need to add 1.0 mentally.
    # Extract x and y coords for plotting
    all_regpn_x, all_regpn_y = zip(*reg_pn[0])
    all_regpn_x = np.asarray(all_regpn_x)
    all_regpn_y = np.asarray(all_regpn_y)

    # Make sure it is closed
    all_regpn_x = np.append(all_regpn_x, all_regpn_x[0])
    all_regpn_y = np.append(all_regpn_y, all_regpn_y[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot region
    ax.plot(all_regpn_x-1.0, all_regpn_y-1.0, color='r')

    # PLot mask
    ax.imshow(regionmask, origin='lower')

    ax.grid(True, lw=0.1)

    plt.show()
    sys.exit(0)
    """

    # ---------------- Extract spectrum ---------------- #
    cube = cube.reshape(hi_full_shape)
    print cube.shape

    # Get cube for flux
    spec_flux_cube = []
    # This will end up being of shape = (num_pix, num_wav)
    # wehre num_pix is pixels within the region
    # num_wav is numer of spectral elements

    for u in range(128):
        for v in range(128):

            if regionmask[u,v]:
                spec_flux_cube.append(cube[:, u, v])

                """
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(np.arange(31), cube[:, u, v])
                plt.show()
                """

    spec_flux_cube = np.array(spec_flux_cube)
    

    sys.exit(0)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.arange(len(spec_flux)), spec_flux)

    plt.show()

    sys.exit(0)

    return spec_wav, spec_flux

def main():

    # ------------------------- Read in data ------------------------- #
    # Read in HI data from Condon et al. 1993
    hi_cube_hdu = fits.open(home + '/Dropbox/Taffy/HI_Condon_et_al/taffy_HIcube.fits')
    hi_cube = hi_cube_hdu[0].data

    # Now read in optical IFU data
    # Only need red channel for Halpha
    obs_r_hdu = fits.open(taffy_data + 'Taffy_R.fits')
    obs_r = obs_r_hdu[0].data

    # Using my own line fits instead of lzifu
    r_line_total = np.load(taffy_extdir + 'halpha_profile_mylinefits.npy')

    # ------------------------- Other preliminaries ------------------------- #
    # create wavelength array
    # I read these data from the corresponding headers
    # red wav arr
    delt_r = 0.3  # i.e. the wav axis is sampled at 0.3A
    red_wav_start = 6165.0
    total_red_res_elem = 2350

    red_wav_arr = [red_wav_start + delt_r*i for i in range(total_red_res_elem)]
    red_wav_arr = np.asarray(red_wav_arr)

    # get mask of all possible not NaN pixels
    all_mask = vcm.get_region_mask('all_possibly_notnan_pixels_new')

    # Make sure that rest wavelength is Halpha wavelength in air!! the IFU data was taken on the ground
    # find line index in wavelength array
    redshift = 0.0145  # average z 
    sys_vel = redshift * speed_of_light  # avg systemic velocity is z*c = 4350 km/s
    halpha_air_wav = 6562.80
    halpha_wav = halpha_air_wav*(1+redshift)
    halpha_idx = np.argmin(abs(red_wav_arr - halpha_wav))

    # ------------------------- Plotting ------------------------- #
    # read in i band SDSS image
    sdss_i, wcs_sdss = vcm.get_sdss('i')

    # read in lzifu output file
    lzifu_hdulist, wcs_lzifu = vcm.get_lzifu_products()

    # plot sdss image
    fig, ax = vcm.plot_sdss_image(sdss_i, wcs_sdss) 

    # Now read in regions file and convert to polygons
    #with fl as open():


    regstr = 'polygon(52.75971,52.679517,55.559628,52.679492,55.559611,55.479493,52.759775,55.479517) # color=red width=2'
    regstr = regstr.split('(')[1].split(')')[0]
    regstr = regstr.split(',')
    reglist = list(regstr)  # this is now a list of floats represented as strings

    reg_pn_list = []
    for k in range(0,8,2):
        reg_pn_list.append([float(reglist[k]), float(reglist[k+1])])

    reg_pn = pg.Polygon(reg_pn_list)

    extract_spec(reg_pn, hi_cube)
    #extract_spec(reg_pn, obs_r)

    #plt.show()

    return None

if __name__ == '__main__':
    main()
    sys.exit(0)