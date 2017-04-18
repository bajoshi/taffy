from __future__ import division

import numpy as np
from astropy.io import fits

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # does not have a trailing slash
desktop = home + '/Desktop/'

taffy_products = '/Volumes/Bhavins_backup/phil/TAFFY/products/'
taffy_data = '/Volumes/Bhavins_backup/phil/TAFFY/data/'
ipac_taffy_dir = home + '/Desktop/ipac/taffy/'

if __name__ == '__main__':
	
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    #filename = '/Volumes/Bhavins_backup/phil/TAFFY/products/TAFFY_2_comp_OI6364.fits'
    #map_hdulist = fits.open(filename)
    #taffy_r = fits.open('/Volumes/Bhavins_backup/phil/TAFFY/data/Taffy_R.fits')

    filename = home + '/Desktop/ipac/Taffy_2/infrared/paper_recent/HIPE15_CALIBRATED_CUBES/TaffyC+_sub.fits' 
    mapname = 'Herschel_C+'
    map_hdulist = fits.open(filename)

    #map_total = map_hdulist[0].data[0]
    #map_vel1 = map_hdulist[0].data[1]
    #map_vel2 = map_hdulist[0].data[2]

    #map_for_cont = map_total
    summed = np.nansum(map_hdulist[1].data, axis=0)
    map_for_cont = summed

    # make contour plot
    x = np.arange(map_for_cont.shape[0])
    y = np.arange(map_for_cont.shape[1])
    X, Y = np.meshgrid(x, y)
    X, Y = X.T, Y.T

    fig = plt.figure()
    ax = fig.add_subplot(111)

    mapmin = np.nanmin(map_for_cont, axis=None)
    mapmax = np.nanmax(map_for_cont, axis=None)
    print mapmin, mapmax
    #levels = np.logspace(mapmin, mapmax , 5)
    # use max=95 or so for OI6300 and OI6364
    levels = np.array([1,3,5,8,9,12,15])

    c = ax.contour(X, Y, map_for_cont, levels=levels, cmap='gist_earth')
    ax.clabel(c, inline=True, inline_spacing=2, fontsize=5, fmt='%1.2f')
    ax.imshow(map_for_cont, origin='lower', cmap='gray_r')

    filename = os.path.basename(filename)
    #mapname = filename.split('.')[0].split('_')[-1]
    fig.savefig(ipac_taffy_dir + '/figures/contour_' + mapname + '.png', dpi=300, bbox_inches='tight')
    #plt.show()

    # total run time
    print "Total time taken --", time.time() - start, "seconds."
    sys.exit(0)