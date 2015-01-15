#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import numpy as np
#import h5py
import matplotlib as mpl
import sys
from plot_aux import H5SFRData, plotLinRegAge

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
    
    K = H5SFRData(h5file)
    
    SFR__Tg = K.get_data_h5('SFR__Tg')
    SFR_Ha__g = K.get_data_h5('SFR_Ha__g')
    aSFRSD_kpc__Trg = K.get_data_h5('aSFRSD_kpc__Trg')
    aSFRSD_Ha_kpc__rg = K.get_data_h5('aSFRSD_Ha_kpc__rg')

    min_pixel_to_plot = 5
    ###################################################################################
    if iT < 0:
        for iT, age in enumerate(K.tSF__T):
            x = np.ma.log10(SFR__Tg[iT])
            y = np.ma.log10(SFR_Ha__g)
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            xran = [-5, 0]
            yran = [-5, 0]
            xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
            ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
            fname = 'SFR_SFRHa_%.2fMyr.png' % (age / 1e6)
            plotLinRegAge(xm, ym, xlabel, ylabel, xran, yran, age, fname)
        
            x = np.ma.log10(aSFRSD_kpc__Trg[iT].flatten())
            y = np.ma.log10(aSFRSD_Ha_kpc__rg.flatten())
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            xran = [-3.5, 1.]
            yran = [-3.5, 1.]
            xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
            fname = 'aSFRSD_aSFRSDHa_%.2fMyr.png' % (age / 1e6)
            plotLinRegAge(xm, ym, xlabel, ylabel, xran, yran, age, fname) 
    else:
        age = K.tSF__T[iT]
        x = np.ma.log10(SFR__Tg[iT])
        y = np.ma.log10(SFR_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        xran = [-5, 0]
        yran = [-5, 0]
        xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
        fname = 'SFR_SFRHa_%.2fMyr.png' % (age / 1e6)
        plotLinRegAge(xm, ym, xlabel, ylabel, xran, yran, age, fname)
    
        x = np.ma.log10(aSFRSD_kpc__Trg[iT].flatten())
        y = np.ma.log10(aSFRSD_Ha_kpc__rg.flatten())
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        xran = [-3.5, 1.]
        yran = [-3.5, 1.]
        xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        fname = 'aSFRSD_aSFRSDHa_%.2fMyr.png' % (age / 1e6)
        plotLinRegAge(xm, ym, xlabel, ylabel, xran, yran, age, fname) 