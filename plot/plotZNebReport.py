#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 15/Jan/2015
#
import numpy as np
#import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import sys
from plot_aux import H5SFRData, plot_text_ax, plot_linreg_params, \
                     OLS_bisector, plotOLSbisectorAxis
from matplotlib.ticker import MultipleLocator
#from matplotlib.ticker import MultipleLocator

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
    
if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    K = H5SFRData(h5file)

    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69,8.27) 
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
       
    xlabel = r'$\log\ Z_\star(t)$' 
    ylabel = r'$\log\ Z_{neb}$'
       
    iT = 0
    
    tSF__T = K.tSF__T
    alogZ_mass__Tg = K.get_data_h5('alogZ_mass__Tg')
    logZ_neb_S06__g = K.get_data_h5('logZ_neb_S06__g')
       
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = alogZ_mass__Tg[iT] 
            y = logZ_neb_S06__g
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            age = tSF__T[iT]
            print 'Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
            #xran = [-3.5, 1]
            yran = [-0.75, 0.25]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            #ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")

            A, B = plotOLSbisectorAxis(ax, xm, ym, 0.92, 0.05, 8)                          
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.92, 8, 'top', 'left')
      
            #ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.2))
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            
            ax.grid(which = 'major')
               
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
                   
            if i == NRows - 1 and j == 3:
                ax.set_xlabel(xlabel)
                   
            if i == 2 and j == 0:
                ax.set_ylabel(ylabel)
               
            iT += 1
       
    f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
    f.savefig('logZ_report.png')
    plt.close(f)
