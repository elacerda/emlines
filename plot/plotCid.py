#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')
ALL_SFR__Tg = get_attrib_h5(h5, 'ALL_SFR__Tg')
ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
ALL_aSFRSD__Trg = get_attrib_h5(h5, 'ALL_aSFRSD__Trg')
ALL_aSFRSD_Ha__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha__rg')

h5.close()
    
for iT, tSF in enumerate(tSF__T):
    x1 = np.ma.log10(ALL_SFR__Tg[iT])
    x1label = r'$\log\ \mathrm{SFR}_\star\ [M_\odot yr^{-1}]$' 
    y1 = np.ma.log10(ALL_SFR_Ha__g)
    y1label = r'$\log\ \mathrm{SFR}_{neb}\ [M_\odot yr^{-1}]$' 
    x2 = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
    x2label = r'$\log\ \Sigma_{\mathrm{SFR}}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    y2 = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten())
    y2label = r'$\log\ \Sigma_{\mathrm{SFR}}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 

    fname = 'SFReSFRSD_%sMyr.png' % str(tSF / 1.e6)
    
    mask1 = ~(x1.mask | y1.mask)
    x1m = x1[mask1]
    y1m = y1[mask1]
    mask2 = ~(x2.mask | y2.mask)
    x2m = x2[mask2]
    y2m = y2[mask2]

    f, axArr = plt.subplots(1, 2)
    f.set_dpi(96)
    f.set_size_inches(20, 10)    
    ax = axArr[0]
    scat = ax.scatter(x1, y1, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.set_xlabel(x1label)
    ax.set_ylabel(y1label)
    ax = axArr[1]
    scat = ax.scatter(x2, y2, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.set_xlabel(x2label)
    ax.set_ylabel(y2label)
    plt.tight_layout()
    f.savefig(fname)
    plt.close(f)
