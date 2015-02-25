#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import sys
import numpy as np
import matplotlib as mpl
from scipy import stats as st
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from CALIFAUtils.plots import plotOLSbisectorAxis, plot_text_ax
from CALIFAUtils.scripts import H5SFRData

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

iTmax = 18

def plotCorrAges(x__T, y, tSF__T, min_N, fname):
    Rs__T = np.ones_like(tSF__T, dtype = np.float)
    
    for iT, age in enumerate(tSF__T):
        x = x__T[iT]
        mask = ~(x.mask | y.mask)
        not_masked = mask.sum()
        
        if not_masked >= min_N:
            xm = x[mask]
            ym = y[mask]
            Rs__T[iT] = st.spearmanr(xm, ym)[0]
        
    f = plt.figure()
    f.set_dpi(100)
    f.set_size_inches(11.69,8.27) 
    ax = f.gca()
    xlabel = r'$\log\ t_\star$ [yr]'
    ylabel = r'coef. spearmann'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim([6.,8.])
    ax.set_ylim([0.,1.])
    
    plt.plot(np.log10(tSF__T), Rs__T, 'k-*')

    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which = 'major')
    
    ax.legend()
    f.savefig(fname)
    plt.close(f)
    
    
def plotSFRgal(x, y, xlabel, ylabel, xlim, ylim, age, N, fname):
    f = plt.figure()
    f.set_dpi(100)
    f.set_size_inches(11.69,8.27) 
    plot_suptitle = '%.2f Myr' % (age/1e6)
    f.suptitle(plot_suptitle)
    ax = f.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if xlim != None:
        ax.set_xlim(xlim)
        
    if ylim != None:
        ax.set_ylim(ylim)

    ax.scatter(x, y, c = 'black', marker = 'o', s = 5, edgecolor = 'none', alpha = 0.5, label='')
        
    txt = 'N: %d' % N
    plot_text_ax(ax, txt, 0.98, 0.14, 14, 'bottom', 'right')

    plotOLSbisectorAxis(ax, x, y, 0.98, 0.07, 14)

    ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3", label = '')

    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which = 'major')
    
    ax.legend()
    f.savefig(fname)
    plt.close(f)
    
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
    
    H = H5SFRData(h5file)
    tSF__T = H.tSF__T[0:iTmax]
    
    for califaID in H.califaIDs:
        SFR__Tz = H.get_prop_gal('SFR__Tg', califaID)
        SFR_Ha__z = H.get_prop_gal('SFR_Ha__g', califaID)
        SFRSD__Tz = H.get_prop_gal('SFRSD__Tg', califaID)
        SFRSD_Ha__z = H.get_prop_gal('SFRSD_Ha__g', califaID)
        aSFRSD__Tr = H.get_prop_gal('aSFRSD__Trg', califaID)
        aSFRSD_Ha__r = H.get_prop_gal('aSFRSD_Ha__rg', califaID)
        
        min_pixel_to_plot = 5
        ###################################################################################
        if _iT < 0:
            print califaID
            
            fname = '%s_SFR_Rs.png' % califaID
            plotCorrAges(SFR__Tz, SFR_Ha__z, tSF__T, min_pixel_to_plot, fname)
            fname = '%s_SFRSD_Rs.png' % califaID
            plotCorrAges(SFRSD__Tz, SFRSD_Ha__z, tSF__T, min_pixel_to_plot, fname)
            fname = '%s_aSFRSD_Rs.png' % califaID
            plotCorrAges(aSFRSD__Tr, aSFRSD_Ha__r, tSF__T, min_pixel_to_plot, fname)
            
            for iT, age in enumerate(tSF__T):
                x = np.ma.log10(SFR__Tz[iT])
                y = np.ma.log10(SFR_Ha__z)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    #xran = [-5, 0]
                    #yran = [-5, 0]
                    xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
                    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
                    fname = '%s_SFR_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname)

                x = np.ma.log10(SFRSD__Tz[iT] * 1e6)
                y = np.ma.log10(SFRSD_Ha__z * 1e6)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    #xran = [-3.5, 1.]
                    #yran = [-3.5, 1.]
                    xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
                    fname = '%s_SFRSD_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname) 

                x = np.ma.log10(aSFRSD__Tr[iT] * 1e6)
                y = np.ma.log10(aSFRSD_Ha__r * 1e6)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    #xran = [-3.5, 1.]
                    #yran = [-3.5, 1.]
                    xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                    fname = '%s_aSFRSD_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname) 
        else:
            print 'XX', califaID

            fname = '%s_SFR_Rs.png' % califaID
            plotCorrAges(SFR__Tz, SFR_Ha__z, tSF__T, min_pixel_to_plot, fname)
            fname = '%s_SFRSD_Rs.png' % califaID
            plotCorrAges(SFRSD__Tz, SFRSD_Ha__z, tSF__T, min_pixel_to_plot, fname)
            fname = '%s_aSFRSD_Rs.png' % califaID
            plotCorrAges(aSFRSD__Tr, aSFRSD_Ha__r, tSF__T, min_pixel_to_plot, fname)
            
            iT = _iT
            age = H.tSF__T[iT]
            x = np.ma.log10(SFR__Tz[iT])
            y = np.ma.log10(SFR_Ha__z)
            mask = ~(x.mask | y.mask)
            not_masked = mask.sum()
            
            if not_masked >= min_pixel_to_plot:
                xm = x[mask]
                ym = y[mask]
                #xran = [-5, 0]
                #yran = [-5, 0]
                xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
                ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
                fname = '%s_SFR_%.2fMyr.png' % (califaID, (age / 1e6))
                plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname)

                x = np.ma.log10(SFRSD__Tz[iT] * 1e6)
                y = np.ma.log10(SFRSD_Ha__z * 1e6)
                mask = ~(x.mask | y.mask)
                not_masked = mask.sum()
                
                if not_masked >= min_pixel_to_plot:
                    xm = x[mask]
                    ym = y[mask]
                    #xran = [-3.5, 1.]
                    #yran = [-3.5, 1.]
                    xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
                    fname = '%s_SFRSD_%.2fMyr.png' % (califaID, (age / 1e6))
                    plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname) 
            
            x = np.ma.log10(aSFRSD__Tr[iT] * 1e6)
            y = np.ma.log10(aSFRSD_Ha__r * 1e6)
            mask = ~(x.mask | y.mask)
            not_masked = mask.sum()
            
            if not_masked >= min_pixel_to_plot:
                xm = x[mask]
                ym = y[mask]
                #xran = [-3.5, 1.]
                #yran = [-3.5, 1.]
                xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
                ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
                fname = '%s_aSFRSD_%.2fMyr.png' % (califaID, (age / 1e6))
                plotSFRgal(xm, ym, xlabel, ylabel, None, None, age, not_masked, fname) 