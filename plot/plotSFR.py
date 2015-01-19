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
                     density_contour, plotRunningStatsAxis, \
                     plotTau, plotSFR, plotScatterColor
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

# zones
ALL_SFR__Tg = get_attrib_h5(h5, 'SFR__Tg')
#ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
#ALL_SFRSD__Tg = get_attrib_h5(h5, 'ALL_SFRSD__Tg')
ALL_SFRSD_kpc__Tg = get_attrib_h5(h5, 'SFRSD_kpc__Tg')
#ALL_SFRSD_Ha__g = get_attrib_h5(h5, 'ALL_SFRSD_Ha__g')
ALL_dist_zone__g = get_attrib_h5(h5, 'dist_zone__g')
#ALL_tau_V__Tg = get_attrib_h5(h5, 'ALL_tau_V__Tg')
ALL_tau_V_neb__g = get_attrib_h5(h5, 'tau_V_neb__g')
#ALL_L_int_Ha__g = get_attrib_h5(h5, 'ALL_L_int_Ha__g')

# Mcor and McorSD by zones for all galaxies
#ALL_Mcor__g = get_attrib_h5(h5, 'ALL_Mcor__g')
#ALL_McorSD__g = get_attrib_h5(h5, 'ALL_McorSD__g')

# galaxy wide quantities replicated by zones
#ALL_Mcor_GAL_zones__g = get_attrib_h5(h5, 'ALL_Mcor_GAL_zones__g')
#ALL_McorSD_GAL_zones__g = get_attrib_h5(h5, 'ALL_McorSD_GAL_zones__g')
#ALL_morfType_GAL_zones__g = get_attrib_h5(h5, 'ALL_morfType_GAL_zones__g')
#ALL_at_flux_GAL_zones__g = get_attrib_h5(h5, 'ALL_at_flux_GAL_zones__g')

# radius
#ALL_aSFRSD_kpc__Trg = get_attrib_h5(h5, 'ALL_aSFRSD_kpc__Trg')
#ALL_aSFRSD_Ha_kpc__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha_kpc__rg')
#ALL_tau_V__Trg = get_attrib_h5(h5, 'ALL_tau_V__Trg')
#ALL_tau_V_neb__rg = get_attrib_h5(h5, 'ALL_tau_V_neb__rg')

#correl_SFR__T = get_attrib_h5(h5, 'correl_SFR__T')
 
h5.close()

for iT,tSF in enumerate(tSF__T):
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # all zones 
#     # x = log10(SFR*(age)) 
#     # y = log10(SFRHa)
#     # z = (tauV*(age) - tauVHa)
#     x = np.ma.log10(ALL_SFR__Tg[iT])
#     y = np.ma.log10(ALL_SFR_Ha__g)
#     z = ALL_tau_V__Tg[iT] - ALL_tau_V_neb__g
#     mask = ~(x.mask | y.mask)
#     xm = x[mask]
#     ym = y[mask]
#     zm = z[mask]
#     xlabel = r'$\log\ SFR_\star\ [M_\odot\ yr^{-1}]$' 
#     ylabel = r'$\log\ SFR_{neb}\ [M_\odot\ yr^{-1}]$'
#     zlabel = r'$\tau_V^\star\ -\ \tau_V^{neb}$' 
#     fname = 'logSFR_logSFR_neb_dTau_age_%sMyr.png' % str(tSF / 1.e6)
#     #xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
#     #ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
#     plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, None, None, None, tSF, fname)
#     #plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
# 
#     # all zones 
#     # x = log10(SFRSD*(age)) 
#     # y = log10(SFRSDHa)
#     # z = (tauV*(age) - tauVHa)
#     x = np.ma.log10(ALL_SFRSD__Tg[iT])
#     y = np.ma.log10(ALL_SFRSD_Ha__g)
#     z = ALL_tau_V__Tg[iT] - ALL_tau_V_neb__g
#     mask = ~(x.mask | y.mask)
#     xm = x[mask]
#     ym = y[mask]
#     zm = z[mask]
#     xlabel = r'$\log\ SFRSD_\star\ [M_\odot\ yr^{-1}\ pc^{-2}]$' 
#     ylabel = r'$\log\ SFRSD_{neb}\ [M_\odot\ yr^{-1}\ pc^{-2}]$'
#     zlabel = r'$\tau_V^\star\ -\ \tau_V^{neb}$' 
#     fname = 'logSFRSD_logSFRSD_neb_dTau_age_%sMyr.png' % str(tSF / 1.e6)
#     #xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
#     #ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
#     plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, None, None, None, tSF, fname)
#     #plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    A_V_neb__g = np.log10(np.exp(1)) * ALL_tau_V_neb__g / 0.4
    SigmaGAS = 15. * A_V_neb__g 
    x = np.log10(SigmaGAS)
    y = np.log10(ALL_SFRSD_kpc__Tg[iT]) + 6.
    z = ALL_dist_zone__g
    mask = ~(x.mask | y.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\log\ \Sigma_{gas} [M_\odot pc^{-2}]$'
    ylabel = r'$\log\ \Sigma_{SFR}^\star\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'zone distance (HLR)'
    fname = 'logSigmaGas_logSFR_distZone_age_%sMyr.png' % str(tSF / 1.e6)
    #xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
    #ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, None, None, fname, None, tSF)
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # x = np.ma.log10(ALL_SFR__Tg[iT])
    # y = np.ma.log10(ALL_SFR_Ha__g)
    # xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
    # ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
    # fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
    # xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
    # ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
    # plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    # y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[:, :].flatten())
    # mask = ~(x.mask | y.mask)
    # xm = x[mask]
    # ym = y[mask]
    # xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
    # ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
    # xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    # ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    # fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr.png' % str(tSF / 1.e6)
    # plotSFR(xm,ym,xlabel,ylabel,xlim,ylim,tSF,fname)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
           
