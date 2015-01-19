#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys
from plot_aux import H5SFRData, plot_text_ax

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    
    tSF__T = H.tSF__T[0:18]
    
    # zones
    tau_V__Tg = H.get_data_h5('tau_V__Tg')
    tau_V_neb__g = H.get_data_h5('tau_V_neb__g')

    ###########################################################################
    ###########################################################################
    ###########################################################################
    for iT,tSF in enumerate(tSF__T):
        xlabel = r'$\tau_V^\star$'
        ylabel = r'$\tau_V^{neb}$'
        fname = 'tauV_tauVNeb_age_%sMyr.png' % str(tSF / 1.e6)
        f = plt.figure()
        f.set_size_inches(10,10)
        ax = f.gca()
        x = tau_V__Tg[iT]
        y = tau_V_neb__g
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        xlim = [0, 2.5]
        ylim = [0, 5]
        scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 4, edgecolor = 'none', alpha = 0.6)
        a = ym.sum() / xm.sum()
        Y = a * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = 'b', ls = '--', lw = 0.5)
        txt = r'y = (%.2f)x $Y_{rms}$:%.2f' %  (a, Yrms)
        plot_text_ax(ax, txt, 0.92, 0.15, 16, 'bottom', 'right', color = 'b')
        txt = '%.2f Myr' % (tSF / 1e6)
        plot_text_ax(ax, txt, 0.05, 0.92, 16, 'top', 'left')

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.125))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.125))
        
        ax.grid(which = 'major')

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.savefig(fname)