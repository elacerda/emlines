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
        #califaID = sys.argv[2]
        try:
            _iT = np.int(sys.argv[2])
        except IndexError:
            _iT = -1
    except IndexError:
        print 'usage: %s HDF5FILE CALIFAID [IT]' % (sys.argv[0])
        exit(1)
    
    K = H5SFRData(h5file)
    
    for califaID in K.califaIDs:
        SFR__Tz = K.get_prop_gal('SFR__Tg', califaID)
        SFR_Ha__z = K.get_prop_gal('SFR_Ha__g', califaID)
        aSFRSD_kpc__Tr = K.get_prop_gal('aSFRSD_kpc__Trg', califaID)
        aSFRSD_Ha_kpc__r = K.get_prop_gal('aSFRSD_Ha_kpc__rg', califaID)
        
        min_pixel_to_plot = 5
        ###################################################################################
        if _iT < 0:
            print califaID
            for iT, age in enumerate(K.tSF__T[0:18]):
                print 'idade: %d' % iT
                x = np.ma.log10(SFR__Tz[iT])
                y = np.ma.log10(SFR_Ha__z)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    print len(xm), len(ym)
                    #xran = [-5, 0]
                    #yran = [-5, 0]
                    xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
                    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
                    fname = '%s_SFR_SFRHa_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotLinRegAge(xm, ym, xlabel, ylabel, None, None, age, fname)

                x = np.ma.log10(aSFRSD_kpc__Tr[iT])
                y = np.ma.log10(aSFRSD_Ha_kpc__r)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    print len(xm), len(ym)
                    #xran = [-3.5, 1.]
                    #yran = [-3.5, 1.]
                    xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                    fname = '%s_aSFRSD_aSFRSDHa_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotLinRegAge(xm, ym, xlabel, ylabel, None, None, age, fname) 
        else:
            print 'XX', califaID
            iT = _iT
            age = K.tSF__T[iT]
            x = np.ma.log10(SFR__Tz[iT])
            y = np.ma.log10(SFR_Ha__z)
            mask = ~(x.mask | y.mask)
            not_masked = mask.sum()
            
            if not_masked >= min_pixel_to_plot:
                xm = x[mask]
                ym = y[mask]
                #xran = [-5, 0]
                #yran = [-5, 0]
                print len(xm), len(ym)
                xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
                ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
                fname = '%s_SFR_SFRHa_%.2fMyr.png' % (califaID, (age / 1e6))
                plotLinRegAge(xm, ym, xlabel, ylabel, None, None, age, fname)
            
            x = np.ma.log10(aSFRSD_kpc__Tr[iT])
            y = np.ma.log10(aSFRSD_Ha_kpc__r)
            mask = ~(x.mask | y.mask)
            not_masked = mask.sum()
            
            if not_masked >= min_pixel_to_plot:
                xm = x[mask]
                ym = y[mask]
                #xran = [-3.5, 1.]
                #yran = [-3.5, 1.]
                xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                fname = '%s_aSFRSD_aSFRSDHa_%.2fMyr.png' % (califaID, (age / 1e6))
                plotLinRegAge(xm, ym, xlabel, ylabel, None, None, age, fname) 