#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import LogNorm
import sys
from plot_aux import H5SFRData, calcRunningStats

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

def plotSK(x, y, z, xlabel, ylabel, zlabel, xlim = False, ylim = False, draw_SK = False, age = False, fname = 'SK.png'):
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    sc = ax.scatter(x, y, c = z, cmap = 'spectral_r', marker = 'o', s = 5., edgecolor = 'none')
    nBox = 20
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    xPrc        = aux[8]
    yPrc        = aux[9]
    ax.plot(xMedian, yMedian, 'k', lw = 2)
    ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
    ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
    if draw_SK:
        ax.plot(np.asarray(ax.get_xlim()), SKzero + SKslope * np.asarray(ax.get_xlim()), 'b--', lw = 2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim != False:
        ax.set_xlim(xlim)
    if ylim != False:
        ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which = 'major')
    cb = f.colorbar(sc)
    cb.set_label(zlabel)
    if age != False:
        ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
    f.savefig(fname)
    plt.close(f)
    
try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)
    
H = H5SFRData(h5file)

tSF__T = H.tSF__T[0:20]

# zones
SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
SFRSD_Ha_kpc__g = H.get_data_h5('SFRSD_Ha__g') * 1e6
dist_zone__g = H.get_data_h5('dist_zone__g')
tau_V__Tg = H.get_data_h5('tau_V__Tg')
tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
tau_V_neb_err__g = H.get_data_h5('tau_V_neb_err__g')
L_int_Ha__g = H.get_data_h5('L_int_Ha__g')
Mcor__g = H.get_data_h5('Mcor__g')
McorSD__g = H.get_data_h5('McorSD__g')
alogZ_mass__g = H.get_data_h5('alogZ_mass__g')
logZ_neb_S06__g = H.get_data_h5('logZ_neb_S06__g')

# galaxy wide quantities replicated by zones
Mcor_GAL_zones__g = H.get_data_h5('Mcor_GAL_zones__g')
McorSD_GAL_zones__g = H.get_data_h5('McorSD_GAL_zones__g')
morfType_GAL_zones__g = H.get_data_h5('morfType_GAL_zones__g')
at_flux_GAL_zones__g = H.get_data_h5('at_flux_GAL_zones__g')

# radius
aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
tau_V__Trg = H.get_data_h5('tau_V__Trg')
tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')

SKzero = np.log10(1.6e-4)
SKslope = 1.4
logSigmaGas = (np.log10(SFRSD_Ha_kpc__g) - SKzero) / SKslope
c = 0 # np.log10(0.2)
logDGR = c + np.log10(tau_V_neb__g) - logSigmaGas
logO_H = logZ_neb_S06__g + np.log10(4.9e-4)

#################################################################################
#################################################################################
#################################################################################
 
for iT,tSF in enumerate(tSF__T):
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = logZ_neb_S06__g
    mask = x.mask | y.mask
    #mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{\star}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
    zlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]' 
    fname = 'tauV_SFRSD_logZneb_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = [-1.5, 0.5]
    ylim = [-3.0, 1]
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)
       
    ###
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = alogZ_mass__g
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{\star}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]'
    fname = 'tauV_SFRSD_alogZmass_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = [-1.5, 0.5]
    ylim = [-3.5, 0]
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)

    ###
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = logZ_neb_S06__g
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{\star}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
    fname = 'tauV_SFRSD_logZneb_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = [-1.5, 0.5]
    ylim = [-3.5, 0]
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)

    ###
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = np.ma.log10(McorSD__g)
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{\star}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
    fname = 'tauV_SFRSD_McorSD_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = [-1.5, 0.5]
    ylim = [-3.5, 0]
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)

    ### 
    x = np.ma.log10(tau_V_neb__g)
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = np.ma.log10(McorSD__g)
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{neb}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
    fname = 'tauVneb_SFRSD_McorSD_age_%sMyr.png' % str(tSF / 1.e6)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # xlim = [-1.5, 0.5]
    # ylim = [-3.5, 0]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    xlim = False
    ylim = False
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)
    
    ###
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = morfType_GAL_zones__g
    mask = x.mask | y.mask | (z <= 7) 
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{neb}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'Morph type'
    fname = 'tauV_SFRSD_MorphType_age_%sMyr.png' % str(tSF / 1.e6)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # xlim = [-1.5, 0.5]
    # ylim = [-3.5, 0]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    xlim = False
    ylim = False
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, False, tSF, fname)

    ###
    x = np.ma.log10(tau_V__Tg[iT])
    y = np.ma.log10(SFRSD__Tg[iT] * 1e6)
    z = np.ma.log10(Mcor_GAL_zones__g)
    mask = x.mask | y.mask 
    xm = x[~mask]
    ym = y[~mask]
    zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{neb}$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\log\ M_\star$ [$M_\odot$]'
    fname = 'tauV_SFRSD_McorGAL_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = [-1.5, 0.5]
    ylim = [-3.5, 0]
    plotSK(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, True, tSF, fname)
     
    ### 
    x = np.ma.log10(tau_V__Trg[iT].flatten())
    y = np.ma.log10(aSFRSD__Trg[iT].flatten() * 1e6)
    #z = np.ma.log10(McorSD__g)
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    #zm = z[~mask]
    xlabel = r'$\log\ \tau_V^{\star}(R)$'
    ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
    fname = 'tauV_aSFRSD_McorSD_age_%sMyr.png' % str(tSF / 1.e6)
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    sc = ax.scatter(xm, ym, c = 'k', cmap = 'spectral_r', marker = 'o', s = 3., edgecolor = 'none', alpha = 0.6)
    nBox = 20
    dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
    aux         = calcRunningStats(xm, ym, dxBox = dxBox, 
                                  xbinIni = xm.min(), xbinFin = xm.max(), 
                                  xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    xPrc        = aux[8]
    yPrc        = aux[9]
    ax.plot(xMedian, yMedian, 'k', lw = 2)
    ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
    ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
     
    #ax.plot(np.asarray(ax.get_xlim()), SKzero + SKslope * np.asarray(ax.get_xlim()), 'b--', lw = 2)
     
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    xlim = [-1.5, 0.5]
    ylim = [-3.5, 0]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which = 'major')
    #cb = f.colorbar(sc)
    #cb.set_label(zlabel)
    ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
    f.savefig(fname)
    plt.close(f)