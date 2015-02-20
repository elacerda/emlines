#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import numpy as np
#import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from califa_scripts import H5SFRData
from plot_aux import plotLinRegAxis, plot_text_ax

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
        try:
            iT = np.int(sys.argv[2])
        except IndexError:
            iT = -1
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    
    SFR__Tg = H.SFR__Tg
    SFR_Ha__g = H.SFR_Ha__g
    SFRSD__Tg = H.SFRSD__Tg
    SFRSD_Ha__g = H.SFRSD_Ha__g
    aSFRSD__Trg = H.aSFRSD__Trg
    aSFRSD_Ha__rg = H.aSFRSD_Ha__rg
    
    if iT < 0:
        tSF__T = H.tSF__T
        ind = range(H.N_T)
    else:
        ind = [ iT ]
        tSF__T = [ H.tSF__T[iT] ]  
        
    ###################################################################################
    for i, age in enumerate(tSF__T):
        j = ind[i]
        f, axArr = plt.subplots(1, 3)
        f.set_size_inches(15,5) 
        ax = axArr[0]
        x = np.ma.log10(SFR__Tg[j])
        y = np.ma.log10(SFR_Ha__g)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \overline{SFR_\star}(t_{SF})\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
        #fname = 'SFR_SFRHa_%.2fMyr.png' % (age / 1e6)
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)

        ax = axArr[1]
        x = np.ma.log10(SFRSD__Tg[j] * 1e6)
        y = np.ma.log10(SFRSD_Ha__g * 1e6)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_{SF})\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
        #fname = 'SFRSD_SFRSDHa_%.2fMyr.png' % (age / 1e6)
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran) 
    
        ax = axArr[2]
        x = np.ma.log10(aSFRSD__Trg[j].flatten() * 1e6)
        y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_{SF}, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
        txt = r'$t_{SF}: %s$ Myr' % str(age / 1.e6)
        plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')

        
        fname = 'SFR_%.2fMyr.png' % (age / 1e6)
        f.subplots_adjust(wspace=0.22, hspace=0, left=0.1, bottom=0.15, right=0.95, top=0.95)
        f.savefig(fname)
        plt.close(f)
 