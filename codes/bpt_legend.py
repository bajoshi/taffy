from __future__ import division

import numpy as np

import sys
import os

import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import pylab

home = os.getenv('HOME')  # Does not have a trailing slash at the end
taffydir = home + "/Desktop/ipac/taffy/"
taffy_extdir = home + '/Desktop/ipac/taffy_lzifu/'

# Solution came from StackOverflow
# See (check all solutions): https://stackoverflow.com/questions/4534480/get-legend-as-a-separate-picture-in-matplotlib
if __name__ == '__main__':

    # create a figure for the data
    figData = pylab.figure()
    ax = pylab.gca()

    # Plot fake data
    # Only the labels, colors, markersize, and symbols matter here
    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='maroon', markersize=6, markeredgecolor='maroon', fmt='x', capsize=0, elinewidth=0.2, \
        label='East Bridge')
    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='darkorange', markersize=4, markeredgecolor='darkorange', fmt='o', zorder=5, capsize=0, elinewidth=0.2, \
        label='West Bridge')
    pylab.errorbar(0.5, 0.5, xerr=0.1, yerr=0.1, \
        color='None', markersize=7, markeredgecolor='darkorange', fmt='o', capsize=0, elinewidth=0.4, \
        label=r'$\left< \mathrm{West\ Bridge} \right>$')

    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='darkgreen', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.2, \
        label='Taffy-N')
    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='darkgreen', markersize=6, markeredgecolor='darkgreen', fmt='+', zorder=5, capsize=0, elinewidth=0.2, \
        label='West Taffy-N')

    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='midnightblue', markersize=3.5, markeredgecolor='None', fmt='o', capsize=0, elinewidth=0.2, \
        label='Taffy-S')
    pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
        color='midnightblue', markersize=4.5, markeredgecolor='midnightblue', fmt='d', zorder=5, capsize=0, elinewidth=0.2, \
        label='Taffy-S nuclear region')
    #pylab.errorbar(np.arange(10), np.arange(10), xerr=np.arange(10), yerr=np.arange(10), \
    #    color='None', markersize=5, markeredgecolor='limegreen', fmt='o', zorder=5, capsize=0, elinewidth=0.2, \
    #    label='Taffy-S nuclear minor axis')

    # create a second figure for the legend
    figLegend = pylab.figure()
    
    # produce a legend for the objects in the other figure
    pylab.figlegend(*ax.get_legend_handles_labels(), loc='center', frameon=False)
    
    # save the two figures to files
    figLegend.savefig(taffy_extdir + "emissionline_diagnostic_legend.png", dpi=300)

    sys.exit(0)