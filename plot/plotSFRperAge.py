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
import CALIFAUtils as C
from CALIFAUtils.plots import plotLinRegAxis, plot_text_ax

#plot_zbins

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
    
    H = C.H5SFRData(h5file)
    
    SFR__Tg = H.SFR__Tg
    SFR_Ha__g = H.SFR_Ha__g
    SFRSD__Tg = H.SFRSD__Tg
    SFRSD_Ha__g = H.SFRSD_Ha__g
    aSFRSD__Trg = H.aSFRSD__Trg
    aSFRSD_Ha__rg = H.aSFRSD_Ha__rg
    integrated_SFR__Tg = H.integrated_SFR__Tg
    integrated_SFR_Ha__g = H.integrated_SFR_Ha__g
    integrated_SFRSD__Tg = H.integrated_SFRSD__Tg
    integrated_SFRSD_Ha__g = H.integrated_SFRSD_Ha__g
    
    if iT < 0:
        tSF__T = H.tSF__T
        ind = range(H.N_T)
    else:
        ind = [ iT ]
        tSF__T = [ H.tSF__T[iT] ]  
        
    ###################################################################################
    for i, tSF in enumerate(tSF__T):
        j = ind[i]
        f, axArr = plt.subplots(2, 3)
        f.set_size_inches(15,12) 
        ax = axArr[0,0]
        x = np.ma.log10(SFR__Tg[j])
        y = np.ma.log10(SFR_Ha__g)
        xm, ym = C.ma_mask_xy(x, y)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ SFR_\star(t_{SF})\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)

        ax = axArr[0,1]
        x = np.ma.log10(SFRSD__Tg[j] * 1e6)
        y = np.ma.log10(SFRSD_Ha__g * 1e6)
        xm, ym = C.ma_mask_xy(x, y)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \Sigma_{SFR}^\star(t_{SF})\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran) 
    
        ax = axArr[0,2]
        x = np.ma.log10(aSFRSD__Trg[j].flatten() * 1e6)
        y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
        xm, ym = C.ma_mask_xy(x, y)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \Sigma_{SFR}^\star(t_{SF}, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)

        ax = axArr[1,0]
        x = np.ma.log10(integrated_SFR__Tg[j])
        y = np.ma.log10(integrated_SFR_Ha__g)
        xm, ym = C.ma_mask_xy(x, y)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ SFR_\star^{int}(t_{SF})\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}^{int}\ [M_\odot yr^{-1}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)

        ax = axArr[1,1]
        x = np.ma.log10(integrated_SFRSD__Tg[j] * 1e6)
        y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
        xm, ym = C.ma_mask_xy(x, y)
        xran = [-5, 1]
        yran = [-5, 1]
        xlabel = r'$\log\ \Sigma_{SFR}^\star(int, t_{SF})\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(int)\ [M_\odot yr^{-1} kpc^{-2}]$'
        plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
        
        ax = axArr[1,2] 
        ax.set_axis_off()
        
        fname = 'SFR_%.2fMyr.png' % (tSF / 1e6)
        suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        f.suptitle(suptitle)
        f.subplots_adjust(wspace=0.22, hspace=0.22, left=0.1, bottom=0.15, right=0.95, top=0.90)
        f.savefig(fname)
        plt.close(f)
 