#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5, plotStatCorreAxis, \
                     density_contour, plotRunningStatsAxis
from scipy import stats as st


mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')

ALL_SFR__Tg = get_attrib_h5(h5, 'ALL_SFR__Tg')
ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
ALL_aSFRSD_kpc__Trg = get_attrib_h5(h5, 'ALL_aSFRSD_kpc__Trg')
ALL_aSFRSD_Ha_kpc__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha_kpc__rg')
ALL_tauV__Trg = get_attrib_h5(h5, 'ALL_tauV__Trg')
ALL_tau_V_neb__rg = get_attrib_h5(h5, 'ALL_tau_V_neb__rg')
ALL_dist_zone__g = get_attrib_h5(h5, 'ALL_dist_zone__g')
ALL_tauV__Tg = get_attrib_h5(h5, 'ALL_tauV__Tg')
ALL_tau_V_neb__g = get_attrib_h5(h5, 'ALL_tauVNeb__g')
ALL_L_int_Ha__g = get_attrib_h5(h5, 'ALL_L_int_Ha__g')
 
h5.close()

for iT,tSF in enumerate(tSF__T):
    x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[:, :].flatten())
    mask = ~(x.mask | y.mask)
    xm = x[mask]
    ym = y[mask]
    xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
    ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
    xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr.png' % str(tSF / 1.e6)
    plotSFR(xm,ym,xlabel,ylabel,xlim,ylim,tSF,fname)
      
    x = np.ma.log10(ALL_SFR__Tg[iT])
    y = np.ma.log10(ALL_SFR_Ha__g)
    xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
    fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
    ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
    plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
          
    x = ALL_tauV__Trg[iT, :, :].flatten()
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    xlabel = r'$\tau_V^\star(R)$'
    ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
      
    x = ALL_tau_V_neb__rg.flatten()
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    xlabel = r'$\tau_V^{neb}(R)$'
    ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    plotTau(x,y,xlabel,ylabel,None,None,tSF,fname)
     
    x = ALL_tauV__Tg[iT]
    y = ALL_tau_V_neb__g
    z = ALL_dist_zone__g
    mask = ~(x.mask | y.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\tau_V^\star$'
    ylabel = r'$\tau_V^{neb}$'
    zlabel = 'pixel distance (HLR)'
    fname = 'tauV_tauVNeb_pixDistHLR_age_%sMyr.png' % str(tSF / 1.e6)
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, tSF, fname)
    
    x = ALL_tauV__Tg[iT]
    y = ALL_tau_V_neb__g
    z = np.ma.log10(ALL_L_int_Ha__g)
    mask = ~(x.mask | y.mask | z.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\tau_V^\star$'
    ylabel = r'$\tau_V^{neb}$'
    zlabel = r'$L_{H\alpha}\ [L_\odot]$'
    fname = 'tauV_tauVNeb_LHa_age_%sMyr.png' % str(tSF / 1.e6)
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, tSF, fname)
    
 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 #    f = plt.figure()
 #    f.set_size_inches(10,10)
 #    ax = f.gca()
 #    sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', vmin = 4., vmax = 6.,  marker = 'o', s = 5., edgecolor = 'none')
 #    binsx = np.linspace(min(xm), max(xm), 21)
 #    binsy = np.linspace(min(ym), max(ym), 21)
 #    density_contour(xm, ym, binsx, binsy, ax = ax, color = 'k')
 #    ax.set_xlabel(xlabel)
 #    ax.set_ylabel(ylabel)
 #    ax.set_xlim(0, 2.5)
 #    ax.set_ylim(0, 2.5)
 #    #plotStatCorreAxis(ax, x, y, 0.03, 0.97, 16)
 #    plotRunningStatsAxis(ax, x, y, 'k')    
 #    cb = f.colorbar(sc)
 #    cb.set_label(zlabel)
 #    ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
 #     
 #    for iT2 in range(0, 8):
 #        x = ALL_tauV__Tg[iT2]
 #        y = ALL_tau_V_neb__g
 #        z = np.ma.log10(ALL_L_int_Ha__g)
 #        mask = ~(x.mask | y.mask | z.mask)
 #        xm = x[mask]
 #        ym = y[mask]
 #        zm = z[mask]
 #        plotRunningStatsAxis(ax, xm, ym, 'k')
 #         
 #    f.savefig(fname)
 #    plt.close(f)
 # 
 #    #iT = [7, 13, 17, 27]    
 #    f, axArr = plt.subplots(2,2)
 #    ax = axArr[0,0]
 #    x = ALL_dist_zone__g
 #    y = ALL_SFR__Tg[7]
 #    mask = ~(x.mask | y.mask)
 #    xm = x[mask]
 #    ym = y[mask]
 #    ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
 #  
 #    ax = axArr[0,1]
 #    x = ALL_dist_zone__g
 #    y = ALL_SFR__Tg[13] - ALL_SFR__Tg[7]
 #    mask = ~(x.mask | y.mask)
 #    xm = x[mask]
 #    ym = y[mask]
 #    ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
 #      
 #    ax = axArr[1,0]
 #    x = ALL_dist_zone__g
 #    y = ALL_SFR__Tg[17] - ALL_SFR__Tg[13]
 #    mask = ~(x.mask | y.mask)
 #    xm = x[mask]
 #    ym = y[mask]
 #    ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
 #      
 #    ax = axArr[1,1]
 #    x = ALL_dist_zone__g
 #    y = ALL_SFR__Tg[27] - ALL_SFR__Tg[17]
 #    mask = ~(x.mask | y.mask)
 #    xm = x[mask]
 #    ym = y[mask]
 #    ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE