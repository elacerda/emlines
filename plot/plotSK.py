#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
#from matplotlib.colors import LogNorm
from CALIFAUtils.scripts import calc_running_stats
from CALIFAUtils.objects import H5SFRData
from CALIFAUtils.plots import plot_zbins

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# mpl.rcParams['font.size']       = 20
# mpl.rcParams['axes.labelsize']  = 20
# mpl.rcParams['axes.titlesize']  = 22
# mpl.rcParams['xtick.labelsize'] = 16
# mpl.rcParams['ytick.labelsize'] = 16 
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

def plotSK(x, y, z, xlabel, ylabel, zlabel, xlim = False, ylim = False, draw_SK = False, age = False, fname = 'SK.png'):
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    sc = ax.scatter(x, y, c = z, cmap = 'spectral_r', marker = 'o', s = 5., edgecolor = 'none')
    nBox = 20
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calc_running_stats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
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
    #if draw_SK:
    #    ax.plot(np.asarray(ax.get_xlim()), SKzero + SKslope * np.asarray(ax.get_xlim()), 'b--', lw = 2)
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

iT_values = [4,10,11,17]

#tSF__T = H.tSF__T[0:20]
tSF__T = np.asarray([H.tSF__T[i] for i in iT_values])

# zones
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
# SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
# SFRSD_Ha_kpc__g = H.get_data_h5('SFRSD_Ha__g') * 1e6
# dist_zone__g = H.get_data_h5('dist_zone__g')
# tau_V__Tg = H.get_data_h5('tau_V__Tg')
# tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
# tau_V_neb_err__g = H.get_data_h5('tau_V_neb_err__g')
# L_int_Ha__g = H.get_data_h5('L_int_Ha__g')
# Mcor__g = H.get_data_h5('Mcor__g')
# McorSD__g = H.get_data_h5('McorSD__g')
# #alogZ_mass__g = H.get_data_h5('alogZ_mass__g')
# logZ_neb_S06__g = H.get_data_h5('logZ_neb_S06__g')
# 
# # galaxy wide quantities replicated by zones
# Mcor_GAL_zones__g = H.reply_arr_by_zones(H.Mcor_GAL__g)
# McorSD_GAL_zones__g = H.reply_arr_by_zones(H.McorSD_GAL__g)
# morfType_GAL_zones__g = H.reply_arr_by_zones(H.morfType_GAL__g)
# #at_flux_GAL_zones__g = H.get_data_h5('at_flux_GAL_zones__g')
# 
# # radius
# aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
# aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
# tau_V__Trg = H.get_data_h5('tau_V__Trg')
# tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')
# 
# x_Y__Tg = H.get_data_h5('x_Y__Tg')
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#################################################################################
#################################################################################
#################################################################################
 
for iT,tSF in enumerate(tSF__T):
    logtau = {
        'logtauV' : dict(v = np.ma.log10(H.tau_V__Tg[iT]), label = r'$\log\ \tau_V^\star$'),
        'logtauVneb' : dict(v = np.ma.log10(H.tau_V_neb__g), label = r'$\log\ \tau_V^{neb}$'),
#        'logtauVGP' : dict(v = np.ma.log10(tau_V_GP__g), label = r'$\log\ \tau_V^{GP}$'),
    }
    SFRSD = {
        'logSFRSD' : dict(v = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'),
        'logSFRSDHa' : dict(v = np.ma.log10(H.SFRSD_Ha__g * 1e6), label = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$'),
    }
    z = {
        'OHIICHIM' : dict(v = H.O_HIICHIM__g, mask = True, label = r'12 + $\log\ O/H$ (HII-CHI-mistry, EPM, 2014)'),
        'logO3N2S06' : dict(v = H.logZ_neb_S06__g + np.log10(4.9e-4) + 12, mask = True, label = r'12 + $\log\ O/H$ (logO3N2, Stasinska, 2006)'),
        'logO3N2M13' : dict(v = H.O_O3N2_M13__g, mask = True, label = r'12 + $\log\ O/H$ (logO3N2, Marino, 2013)'),
        'alogZmass' : dict(v = H.alogZ_mass__Ug[-1], mask = True, label = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)),
        'xY' : dict(v = 100. * H.x_Y__Tg[iT], label = r'$x_Y$'),
        'logMcorSD' : dict(v = np.ma.log10(H.McorSD__g), mask = True, label = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'),
        'logMcorSDiT' : dict(v = np.ma.log10(H.McorSD__Tg[iT]), mask = True, label = r'$\log\ \mu_\star(t)$ [$M_\odot\ pc^{-2}$]'),
        'logMcorGAL' : dict(v = np.ma.log10(H.reply_arr_by_zones(H.Mcor_GAL__g)), label = r'$\log\ M^{gal}_\star$ [$M_\odot$]'),
        'morfType' : dict(v = H.reply_arr_by_zones(H.morfType_GAL__g), label = 'Morph. Type'),
    }
    for sfrsdname, sfrsdval in SFRSD.iteritems():
        for tauname, tauval in logtau.iteritems():
            for zname, zval in z.iteritems():
                plot_zbins(debug = True,
                           x = tauval['v'],
                           xlabel = tauval['label'],
                           y = sfrsdval['v'],
                           ylabel = sfrsdval['label'],
                           z = zval['v'],
                           zlabel = zval['label'],
                           zmask = zval.get('mask', False), 
                           xlim = [-1.5, 0.5],
                           ylim = [-3.0, 1],
                           #zlimprc = [ 2, 98 ],
                           zbins = 4,
                           zbins_rs_gaussian_smooth = True,
                           zbins_rs_gs_fwhm = 0.4,
                           kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                           running_stats = True,
                           rs_gaussian_smooth = True,
                           rs_gs_fwhm = 0.4,
                           kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                           rs_errorbar = False,
                           suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1.e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                           filename = '%s_%s_%s_%.2fMyr.png' % (tauname, sfrsdname, zname, tSF / 1e6),
                           x_major_locator = 0.5,
                           x_minor_locator = 0.125,
                           y_major_locator = 0.5,
                           y_minor_locator = 0.125,
                           kwargs_legend = dict(fontsize = 8),
                    )
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#       
#     ### 
#     x = np.ma.log10(tau_V__Trg[iT].flatten())
#     y = np.ma.log10(aSFRSD__Trg[iT].flatten() * 1e6)
#     #z = np.ma.log10(McorSD__g)
#     mask = x.mask | y.mask
#     xm = x[~mask]
#     ym = y[~mask]
#     #zm = z[~mask]
#     xlabel = r'$\log\ \tau_V^{\star}(R)$'
#     ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
#     fname = 'logtauV_logaSFRSD_age_%sMyr.png' % str(tSF / 1.e6)
#     f = plt.figure()
#     f.set_size_inches(10,8)
#     ax = f.gca()
#     sc = ax.scatter(xm, ym, c = 'k', marker = 'o', s = 3., edgecolor = 'none', alpha = 0.6)
#     nBox = 20
#     dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
#     aux         = calc_running_stats(xm, ym, dxBox = dxBox, 
#                                      xbinIni = xm.min(), xbinFin = xm.max(), 
#                                      xbinStep = dxBox)
#     xbinCenter  = aux[0]
#     xMedian     = aux[1]
#     xMean       = aux[2]
#     xStd        = aux[3]
#     yMedian     = aux[4]
#     yMean       = aux[5]
#     yStd        = aux[6]
#     nInBin      = aux[7]
#     xPrc        = aux[8]
#     yPrc        = aux[9]
#     ax.plot(xMedian, yMedian, 'k', lw = 2)
#     ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
#     ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
#     #ax.plot(np.asarray(ax.get_xlim()), SKzero + SKslope * np.asarray(ax.get_xlim()), 'b--', lw = 2)
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     xlim = [-1.5, 0.5]
#     ylim = [-3.0, 1]
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.xaxis.set_major_locator(MultipleLocator(0.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.yaxis.set_major_locator(MultipleLocator(0.5))
#     ax.yaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.grid(which = 'major')
#     #cb = f.colorbar(sc)
#     #cb.set_label(zlabel)
#     ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
#     f.savefig(fname)
#     plt.close(f)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
