# simpsons rule test

from __future__ import division

import numpy as np
from scipy.integrate import simps 

import matplotlib.pyplot as plt

import sys

if __name__ == '__main__':
    
    bridge_b1_cp = np.genfromtxt('/Users/bhavinjoshi/Desktop/ipac/taffy/herschel_spectra_from_phil/bridge_b1_sum_c+.dat-wav.txt',\
     dtype=None, names=['wav','flux'])

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(bridge_b1_cp['wav'], bridge_b1_cp['flux'])
    ax.axhline(y=0, color='k', linestyle='--')

    plt.show()

    # eyeballed limits for integration
    # we only want to integrate where the line is actually there
    lam_intg_low = 159.6 
    lam_intg_high = 160.4

    lam = bridge_b1_cp['wav']
    flux = bridge_b1_cp['flux']

    low_idx = np.argmin(abs(lam - lam_intg_low))
    high_idx = np.argmin(abs(lam - lam_intg_high))

    flux_intg = flux[low_idx:high_idx+1]
    wav_intg = lam[low_idx:high_idx+1]

    print simps(y=flux_intg, x=wav_intg)

    sys.exit(0)