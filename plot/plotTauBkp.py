#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import sys
from plot_aux import H5SFRData, plotRunningStatsAxis, \
                     plotScatterColor, calcRunningStats

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)
    
H = H5SFRData(h5file)

tSF__T = H.tSF__T[0:20]

# zones
SFR__Tg = H.get_data_h5('SFR__Tg')
SFR_Ha__g = H.get_data_h5('SFR_Ha__g')
SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
SFRSD_Ha_kpc__g = H.get_data_h5('SFRSD_Ha_kpc__g')
dist_zone__g = H.get_data_h5('dist_zone__g')
tau_V__Tg = H.get_data_h5('tau_V__Tg')
tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
L_int_Ha__g = H.get_data_h5('L_int_Ha__g')
F_obs_Ha__g = H.get_data_h5('F_obs_Ha__g')
Mcor__g = H.get_data_h5('Mcor__g')
McorSD__g = H.get_data_h5('McorSD__g')

# galaxy wide quantities replicated by zones
Mcor_GAL_zones__g = H.get_data_h5('Mcor_GAL_zones__g')
McorSD_GAL_zones__g = H.get_data_h5('McorSD_GAL_zones__g')
morfType_GAL_zones__g = H.get_data_h5('morfType_GAL_zones__g')
at_flux_GAL_zones__g = H.get_data_h5('at_flux_GAL_zones__g')

# radius
aSFRSD_kpc__Trg = H.get_data_h5('aSFRSD_kpc__Trg')
aSFRSD_Ha_kpc__rg = H.get_data_h5('aSFRSD_Ha_kpc__rg')
tau_V__Trg = H.get_data_h5('tau_V__Trg')
tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')

correl_SFR__T = H.get_data_h5('correl_SFR__T')

###########################################################################
###########################################################################
###########################################################################
iT = 4
x = np.ma.log10(SFR_Ha__g)
y1 = tau_V__Tg[iT]
y2 = tau_V_neb__g
mask = ~(x.mask | y1.mask | y2.mask)
xm = x[mask]
y1m = y1[mask]
y2m = y2[mask]    

xlabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'
ylabel = r'$\tau_V$'
fname = 'SFRNeb_tauV_tauVNeb_bestCorrAge.png'
f = plt.figure()
f.set_size_inches(10,10)
ax1 = plt.subplot2grid((3,1), (0, 0), rowspan = 2)
ax2 = plt.subplot2grid((3,1), (2, 0), sharex = ax1)
binsy = np.linspace(0., 4., 51)
binsx = np.linspace(xm.min(), xm.max(), 51)
ax1.scatter(xm, y1m, c = 'b', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^\star$')
#ax.hist2d(xm, y1m, bins = 50, cmap = 'Greys')
ax1.scatter(xm, y2m, c = 'r', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^{neb}$')
#ax.hist2d(xm, y2m, bins = 50, cmap = 'Blues')
#c = [ 'b', 'g', 'r', 'c', 'm', 'y' ]

for i in range(0, 10):
    tSF = tSF__T[i]
    yy1 = tau_V__Tg[i]
    yy2 = tau_V_neb__g
    mask = ~(x.mask | yy1.mask | yy2.mask)
    xm = x[mask]
    yym = yy2[mask] - yy1[mask]
    #plotRnningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '-', errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%c-' % c[i], errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax, xm, yy2m, r'$\tau_V^{neb}\ %.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%co-' % c[i], errorbar = False, nBox = 20)
    dxBox       = (xm.max() - xm.min()) / (20 - 1.)
    aux         = calcRunningStats(xm, yym, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]

    ax2.plot(xMedian, yMedian, '-', lw = 2, label = r'$%.2fMyr$' % (tSF / 1.e6))
    
ax1.set_ylim(0., 4.)
ax1.set_ylabel(ylabel)
ax1.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
ax1.legend(fontsize = 14, frameon = False, loc = 'upper right')

ax2.legend(bbox_to_anchor=(0.,3.), fontsize = 14, frameon = False, loc = 'upper left')
ax2.set_ylabel(r'$\tau_V^{neb} - \tau_V^\star$')
ax2.set_xlabel(xlabel)
#ax2.set_ylim([0, 4.49])

plt.setp(ax1.get_xticklabels(), visible = False)
plt.subplots_adjust(wspace=0, hspace=0)
f.savefig(fname)
plt.close(f)

###########################################################################
###########################################################################
###########################################################################
iT = 4
x = np.ma.log10(SFRSD_Ha_kpc__g)
y1 = tau_V__Tg[iT]
y2 = tau_V_neb__g
mask = ~(x.mask | y1.mask | y2.mask)
xm = x[mask]
y1m = y1[mask]
y2m = y2[mask]    

xlabel = r'$\log\ \Sigma_SFR^{neb}\ $[M${}_\odot\ $yr${}^{-1}\ $kpc${}^{-1}]$'
ylabel = r'$\tau_V$'
fname = 'SFRSDNeb_tauV_tauVNeb_bestCorrAge.png'
f = plt.figure()
f.set_size_inches(10,10)
ax1 = plt.subplot2grid((3,1), (0, 0), rowspan = 2)
ax2 = plt.subplot2grid((3,1), (2, 0), sharex = ax1)
binsy = np.linspace(0., 4., 51)
binsx = np.linspace(min(xm), max(xm), 51)
ax1.scatter(xm, y1m, c = 'b', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^\star$')
#ax.hist2d(xm, y1m, bins = 50, cmap = 'Greys')
ax1.scatter(xm, y2m, c = 'r', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^{neb}$')
#ax.hist2d(xm, y2m, bins = 50, cmap = 'Blues')
#c = [ 'b', 'g', 'r', 'c', 'm', 'y' ]

for i in range(0, 10):
    tSF = tSF__T[i]
    yy1 = tau_V__Tg[i]
    yy2 = tau_V_neb__g
    mask = ~(x.mask | yy1.mask | yy2.mask)
    xm = x[mask]
    yym = yy2[mask] - yy1[mask]
    #plotRunningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%c-' % c[i], errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '-', errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax, xm, yy2m, r'$\tau_V^{neb}\ %.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%co-' % c[i], errorbar = False, nBox = 20)
    dxBox       = (xm.max() - xm.min()) / (20 - 1.)
    aux         = calcRunningStats(xm, yym, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]

    ax2.plot(xMedian, yMedian, '-', lw = 2, label = r'$%.2fMyr$' % (tSF / 1.e6))
    
ax1.set_ylim(0, 4.)
ax1.set_ylabel(ylabel)
ax1.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
ax1.legend(fontsize = 14, frameon = False, loc = 'upper right')

ax2.legend(bbox_to_anchor=(0.,3.), fontsize = 14, frameon = False, loc = 'upper left')
ax2.set_ylabel(r'$\tau_V^{neb} - \tau_V^\star$')
ax2.set_xlabel(xlabel)
ax2.set_ylim([-0.5, 1.99])

plt.setp(ax1.get_xticklabels(), visible = False)
plt.subplots_adjust(wspace=0, hspace=0)
f.savefig(fname)
plt.close(f)

###########################################################################
###########################################################################
###########################################################################

iT = 4
x = np.ma.log10(F_obs_Ha__g)
y1 = tau_V__Tg[iT]
y2 = tau_V_neb__g
mask = ~(x.mask | y1.mask | y2.mask)
xm = x[mask]
y1m = y1[mask]
y2m = y2[mask]    

xlabel = r'$\log\ F_{obs}^{H\alpha}\ $[erg cm${}^{-2}$ s${}^{-1}]$'
ylabel = r'$\tau_V$'
fname = 'FobsHa_tauV_tauVNeb_bestCorrAge.png'
f = plt.figure()
f.set_size_inches(10,10)
ax1 = plt.subplot2grid((3,1), (0, 0), rowspan = 2)
ax2 = plt.subplot2grid((3,1), (2, 0), sharex = ax1)
binsy = np.linspace(0., 4., 51)
binsx = np.linspace(min(xm), max(xm), 51)
ax1.scatter(xm, y1m, c = 'b', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^\star$')
#ax.hist2d(xm, y1m, bins = 50, cmap = 'Greys')
ax1.scatter(xm, y2m, c = 'r', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^{neb}$')
#ax.hist2d(xm, y2m, bins = 50, cmap = 'Blues')
#c = [ 'b', 'g', 'r', 'c', 'm', 'y' ]

for i in range(0, 10):
    tSF = tSF__T[i]
    yy1 = tau_V__Tg[i]
    yy2 = tau_V_neb__g
    mask = ~(x.mask | yy1.mask | yy2.mask)
    xm = x[mask]
    yym = yy2[mask] - yy1[mask]
    #plotRunningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '-', errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax2, xm, yym, r'$%.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%c-' % c[i], errorbar = False, nBox = 20)
    #plotRunningStatsAxis(ax, xm, yy2m, r'$\tau_V^{neb}\ %.2fMyr$' % (tSF / 1.e6), plot_stats = 'median', color = '%co-' % c[i], errorbar = False, nBox = 20)
    nBox = 20
    dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
    aux         = calcRunningStats(xm, yym, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]

    ax2.plot(xMedian, yMedian, '-', lw = 2, label = r'$%.2fMyr$' % (tSF / 1.e6))
    
ax1.set_ylim(0, 4.)
ax1.set_ylabel(ylabel)
ax1.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
ax1.legend(fontsize = 14, frameon = False, loc = 'upper right')

ax2.legend(bbox_to_anchor=(0.,3.), fontsize = 14, frameon = False, loc = 'upper left')
ax2.set_ylabel(r'$\tau_V^{neb} - \tau_V^\star$')
ax2.set_xlabel(xlabel)
#ax2.set_ylim([0., 0.79])

plt.setp(ax1.get_xticklabels(), visible = False)
plt.subplots_adjust(wspace=0, hspace=0)
f.savefig(fname)
plt.close(f)

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# from matplotlib.colors import LogNorm
# iT = [7, 13, 17, 27]    
# f, axArr = plt.subplots(2,2)
# ax = axArr[0,0]
# x = dist_zone__g
# y = SFR__Tg[iT[0]]
# mask = ~(x.mask | y.mask)
# xm = x[mask]
# ym = y[mask]
# ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
# 
# ax = axArr[0,1]
# x = dist_zone__g
# y = SFR__Tg[iT[1]] - SFR__Tg[iT[0]]
# mask = ~(x.mask | y.mask)
# xm = x[mask]
# ym = y[mask]
# ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
#     
# ax = axArr[1,0]
# x = dist_zone__g
# y = SFR__Tg[iT[2]] - SFR__Tg[iT[1]]
# mask = ~(x.mask | y.mask)
# xm = x[mask]
# ym = y[mask]
# ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
#     
# ax = axArr[1,1]
# x = dist_zone__g
# y = SFR__Tg[iT[3]] - SFR__Tg[iT[2]]
# mask = ~(x.mask | y.mask)
# xm = x[mask]
# ym = y[mask]
# ax.hist2d(xm[(ym > -0.5) & (ym < 0.5)], ym[(ym > -0.5) & (ym < 0.5)], bins = 100, norm = LogNorm())
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
for iT,tSF in enumerate(tSF__T):
    x = tau_V__Tg[iT]
    y = tau_V_neb__g
    z = dist_zone__g
    mask = ~(x.mask | y.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\tau_V^\star$'
    ylabel = r'$\tau_V^{neb}$'
    zlabel = 'pixel distance (HLR)'
    fname = 'tauV_tauVNeb_pixDistHLR_age_%sMyr.png' % str(tSF / 1.e6)
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, [0, 2.5], [0, 2.5], fname, [0, 2.0], tSF)
    
    x = tau_V__Tg[iT]
    y = tau_V_neb__g
    z = np.ma.log10(L_int_Ha__g)
    mask = ~(x.mask | y.mask | z.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\tau_V^\star$'
    ylabel = r'$\tau_V^{neb}$'
    zlabel = r'$\log\ L_{H\alpha}\ [L_\odot]$'
    fname = 'tauV_tauVNeb_LHa_age_%sMyr.png' % str(tSF / 1.e6)
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, [0, 2.5], [0, 2.5], fname, [4.,6.], tSF)

    x = tau_V__Tg[iT]
    y = tau_V_neb__g
    z = np.ma.log10(F_obs_Ha__g)
    mask = ~(x.mask | y.mask | z.mask)
    xm = x[mask]
    ym = y[mask]
    zm = z[mask]
    xlabel = r'$\tau_V^\star$'
    ylabel = r'$\tau_V^{neb}$'
    zlabel = r'$\log\ F_{H\alpha}^{obs}\ [L_\odot]$'
    fname = 'tauV_tauVNeb_FobsHa_age_%sMyr.png' % str(tSF / 1.e6)
    plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, [0, 2.5], [0, 2.5], fname, [-15.5,-14.0], tSF)


    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # x = tau_V__Trg[iT, :, :].flatten()
    # y = np.ma.log10(aSFRSD_Ha_kpc__rg.flatten() / aSFRSD_kpc__Trg[iT, :, :].flatten())
    # xlabel = r'$\tau_V^\star(R)$'
    # ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    # fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    # plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
    #    
    # x = tau_V_neb__rg.flatten()
    # y = np.ma.log10(aSFRSD_Ha_kpc__rg.flatten() / aSFRSD_kpc__Trg[iT, :, :].flatten())
    # xlabel = r'$\tau_V^{neb}(R)$'
    # ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    # fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    # plotTau(x,y,xlabel,ylabel,None,None,tSF,fname)
    # f = plt.figure()
    # f.set_size_inches(10,10)
    # ax = f.gca()
    # sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', vmin = 4., vmax = 6.,  marker = 'o', s = 5., edgecolor = 'none')
    # binsx = np.linspace(min(xm), max(xm), 21)
    # binsy = np.linspace(min(ym), max(ym), 21)
    # density_contour(xm, ym, binsx, binsy, ax = ax, color = 'k')
    # ax.set_xlabel(xlabel)
    # ax.set_ylabel(ylabel)
    # ax.set_xlim(0, 2.5)
    # ax.set_ylim(0, 2.5)
    # #plotStatCorreAxis(ax, x, y, 0.03, 0.97, 16)
    # plotRunningStatsAxis(ax, x, y, 'k')    
    # cb = f.colorbar(sc)
    # cb.set_label(zlabel)
    # ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
    #    
    # for iT2 in range(0, 8):
    #     x = tau_V__Tg[iT2]
    #     y = tau_V_neb__g
    #     z = np.ma.log10(L_int_Ha__g)
    #     mask = ~(x.mask | y.mask | z.mask)
    #     xm = x[mask]
    #     ym = y[mask]
    #     zm = z[mask]
    #     plotRunningStatsAxis(ax, xm, ym, 'k')
    #        
    # f.savefig(fname)
    # plt.close(f)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    

