#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys
from califa_scripts import H5SFRData
from plot_aux import plotOLSbisectorAxis, plot_text_ax


mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

iT_values = [4, 10, 11, 17]

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    tSF__T = np.asarray([H.tSF__T[i] for i in iT_values])
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    
    # zones
    tau_V__Tg = H.get_data_h5('tau_V__Tg')
    tau_V_neb__g = H.get_data_h5('tau_V_neb__g')

    ###########################################################################
    ###########################################################################
    ###########################################################################
    for iT,tSF in enumerate(tSF__T):        
        x = tau_V__Tg[iT]
        y = tau_V_neb__g
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xlabel = r'$\tau_V^\star$'
        ylabel = r'$\tau_V^{neb}$'
        fname = 'tauV_tauVNeb_age_%sMyr.png' % str(tSF / 1.e6)
        xlim = [0, 1.5]
        ylim = [0, 2.5]
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 5., edgecolor = 'none')
        a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = 'b', ls = '--', lw = 0.5)
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = 'b')
        #####################
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(which = 'major')
        txt = r'%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
        f.suptitle(txt, fontsize = 14)
        #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        f.savefig(fname)
        plt.close(f)
