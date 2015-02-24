#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from CALIFAUtils.scripts import H5SFRData 

def iTGen(tSF__T, iT_values = [ 10, 11, 17 ]):
    for iT in iT_values:
        yield iT, tSF__T[iT]

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE CALIFAID' % (sys.argv[0])
    exit(1)

H = H5SFRData(h5file)
tSF__T = H.tSF__T
iT_values = [ 10, 11, 17 ]

for iT, tSF in iTGen(H.tSF__T, iT_values):
    x1 = np.ma.log10(H.tau_V__Trg[iT])
    y1 = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6) 
    x1label = r'$\log\ \tau_V^{\star}(R)$'
    y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'

    x2 = np.ma.log10(H.tau_V__Tg[iT])
    y2 = np.ma.log10(H.SFRSD__Tg[iT] * 1e6) 
    x2label = r'$\log\ \tau_V^{\star}$'
    y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'

    fname = 'SKeSKradius_%sMyr.png' % str(tSF / 1.e6)
    
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
