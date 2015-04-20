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
#from CALIFAUtils.plots import plotLinRegAxis, plot_text_ax

#plot_zbins
mask_radius = True
RNuc = 0.5

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
            iTs = np.int(sys.argv[2])
        except IndexError:
            iTs = -1
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = C.H5SFRData(h5file)
        
    if iTs < 0:
        tSF__T = H.tSF__T
        ind = range(H.N_T)
    else:
        ind = [ iTs ]
        tSF__T = [ H.tSF__T[iTs] ]  
        
    ###################################################################################
    for i, tSF in enumerate(tSF__T):
        iT = ind[i]
        
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][10, :] + H.aSFRSD__Trg[iT][9, :] / 2.)
        # SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[10, :] + H.aSFRSD_Ha__rg[9, :] / 2.)
        # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # # SFRSD_norm_GAL__g = H.integrated_SFRSD__Tg[iT]
        # # SFRSD_Ha_norm_GAL__g = H.integrated_SFRSD_Ha__g
        # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
        # aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

        x1 = H.SFR__Tg[iT]
        y1 = H.SFR_Ha__g
        x1label = r'$\log\ SFR_\star(t_{SF})\ [M_\odot yr^{-1}]$' 
        y1label = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
        
        x2 = H.SFRSD__Tg[iT]
        y2 = H.SFRSD_Ha__g
        x2label = r'$\log\ \Sigma_{SFR}^\star(t_{SF})\ [M_\odot yr^{-1} kpc^{-2}]$' 
        y2label = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'
        
        x3 = H.aSFRSD__Trg[iT]
        y3 = H.aSFRSD_Ha__rg
        x3label = r'$\log\ \Sigma_{SFR}^\star(t_{SF}, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        y3label = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'

        if mask_radius is True:
            x1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
            y1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
            x2[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
            y2[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
            x3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
            y3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
            suptitle = r'NGals:%d  R > %.2fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            fname = 'SFR_maskradius_%.2fMyr.png' % (tSF / 1e6)
        else:
            suptitle = r'NGals:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            fname = 'SFR_%.2fMyr.png' % (tSF / 1e6)
        
        #f, axArr = plt.subplots(2, 3)
        #f.set_size_inches(15,12) 
        f, axArr = plt.subplots(1, 3)
        f.set_size_inches(15,8) 
        kw = C.plot_zbins(
            return_kwargs = True,
            f = f,
            ax = axArr[0],
            debug = True,
            x = np.ma.log10(x1),
            y = np.ma.log10(y1),
            xlabel = x1label,
            ylabel = y1label,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # z = H.zone_dist_HLR__g,
            # zmask = None,
            # zlim = [RNuc, 3],
            # zlabel = r'R (HLR)',
            # cb = None,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            xlim = [-5, 1],
            ylim = [-5, 1],
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            y_major_locator = 1.,
            y_minor_locator = 0.2,
            ols = True,
            spearmanr = True,
            kwargs_figure = dict(figsize = (10, 8), dpi = 100),
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
            kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
            kwargs_ols_plot = dict(c = 'r', ls = '-', lw = 2, label = 'OLS'),
        )
        kw['ax'].plot(kw['ax'].get_xlim(), kw['ax'].get_xlim(), ls = '--', label = '', c = 'k')    
        kw = C.plot_zbins(
            return_kwargs = True,
            f = f,
            ax = axArr[1],
            debug = True,
            x = np.ma.log10(x2 * 1e6),
            y = np.ma.log10(y2 * 1e6),
            xlabel = x2label,
            ylabel = y2label,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # z = H.zone_dist_HLR__g,
            # zmask = None,
            # zlabel = r'R (HLR)',
            # zlim = [RNuc, 3],
            # cb = None,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            xlim = [-5, 1],
            ylim = [-5, 1],
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            y_major_locator = 1.,
            y_minor_locator = 0.2,
            ols = True,
            spearmanr = True,
            kwargs_figure = dict(figsize = (10, 8), dpi = 100),
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
            kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
            kwargs_ols_plot = dict(c = 'r', ls = '-', lw = 2, label = 'OLS'),
        )    
        kw['ax'].plot(kw['ax'].get_xlim(), kw['ax'].get_xlim(), ls = '--', label = '', c = 'k')
        kw = C.plot_zbins(
            return_kwargs = True,
            f = f,
            ax = axArr[2],
            debug = True,
            x = np.ma.log10(x3 * 1e6).flatten(),
            y = np.ma.log10(y3 * 1e6).flatten(),
            xlabel = x3label,
            ylabel = y3label,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # z = (H.RbinCenter__r[..., np.newaxis] * np.ones_like(H.aSFRSD_Ha__rg)).flatten(),
            # zmask = None,
            # zlabel = r'R (HLR)',
            # zlim = [RNuc, 3],
            # cb = True,
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            xlim = [-5, 1],
            ylim = [-5, 1],
            x_major_locator = 1.,
            x_minor_locator = 0.2,
            y_major_locator = 1.,
            y_minor_locator = 0.2,
            ols = True,
            spearmanr = True,
            kwargs_figure = dict(figsize = (10, 8), dpi = 100),
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8, label = ''),
            kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
            kwargs_ols_plot = dict(c = 'r', ls = '-', lw = 2, label = 'OLS'),
        )    
        kw['ax'].plot(kw['ax'].get_xlim(), kw['ax'].get_xlim(), ls = '--', label = '', c = 'k')
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         x = np.ma.log10(x1)
#         y = np.ma.log10(y1)
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-5, 1]
#         yran = [-5, 1]
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
# 
#         ax = axArr[0,1]
#         x = np.ma.log10(x2 * 1e6)
#         y = np.ma.log10(y2 * 1e6)
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-5, 1]
#         yran = [-5, 1]
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran) 
#     
#         ax = axArr[0,2]
#         x = np.ma.log10(x3.flatten() * 1e6)
#         y = np.ma.log10(y3.flatten() * 1e6)
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-5, 1]
#         yran = [-5, 1]
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         ax = axArr[1,0]
#         x = np.ma.log10(integrated_SFR__Tg[iT])
#         y = np.ma.log10(integrated_SFR_Ha__g)
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-5, 1]
#         yran = [-5, 1]
#         xlabel = r'$\log\ SFR_\star^{int}(t_{SF})\ [M_\odot yr^{-1}]$' 
#         ylabel = r'$\log\ SFR_{neb}^{int}\ [M_\odot yr^{-1}]$'
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
# 
#         ax = axArr[1,1]
#         x = np.ma.log10(integrated_SFRSD__Tg[iT] * 1e6)
#         y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-5, 1]
#         yran = [-5, 1]
#         xlabel = r'$\log\ \Sigma_{SFR}^\star(int, t_{SF})\ [M_\odot yr^{-1} kpc^{-2}]$' 
#         ylabel = r'$\log\ \Sigma_{SFR}^{neb}(int)\ [M_\odot yr^{-1} kpc^{-2}]$'
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
#         
#         ax = axArr[1,2] 
#         x = np.ma.log10(aSFRSD_norm__rg.flatten())
#         y = np.ma.log10(aSFRSD_Ha_norm__rg.flatten())
#         xm, ym = C.ma_mask_xy(x, y)
#         xran = [-2, 2]
#         yran = [-2, 2]
#         xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(1HLR)}$'
#         ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}(R)}{\Sigma_{SFR}^{neb}(1HLR)}$' 
#         plotLinRegAxis(ax, xm, ym, xlabel, ylabel, xran, yran)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

        f.suptitle(suptitle)
        f.subplots_adjust(wspace=0.22, hspace=0.22, left=0.1, bottom=0.15, right=0.95, top=0.90)
        f.savefig(fname)
        plt.close(f)
 