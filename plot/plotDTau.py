#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys
from califa_scripts import H5SFRData
from plot_aux import plotScatterColor, calcRunningStats, \
                     plotOLSbisectorAxis, plot_text_ax

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

iT_values = [4, 10, 11, 17]

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    
    tSF__T = np.asarray([H.tSF__T[i] for i in iT_values])
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    
    # zones
    x_Y__Tg = H.get_data_h5('x_Y__Tg')
    tau_V__Tg = H.get_data_h5('tau_V__Tg')
    tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
    EW_Ha__g = H.get_data_h5('EW_Ha__g')
    dist_zone__g = H.get_data_h5('dist_zone__g')
    L_int_Ha__g = H.get_data_h5('L_int_Ha__g')
    SFR_Ha__g = H.get_data_h5('SFR_Ha__g')
    SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # SFR__Tg = H.get_data_h5('SFR__Tg')
    # SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
    # SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
    # Mcor__g = H.get_data_h5('Mcor__g')
    # McorSD__g = H.get_data_h5('McorSD__g')
    # logZ_neb_S06__g = H.get_data_h5('logZ_neb_S06__g')
    # 
    # # galaxy wide quantities replicated by zones
    # Mcor_GAL_zones__g = H.get_data_h5('Mcor_GAL_zones__g')
    # McorSD_GAL_zones__g = H.get_data_h5('McorSD_GAL_zones__g')
    # morfType_GAL_zones__g = H.get_data_h5('morfType_GAL_zones__g')
    # at_flux_GAL_zones__g = H.get_data_h5('at_flux_GAL_zones__g')
    # 
    # # radius
    # aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
    # aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
    # tau_V__Trg = H.get_data_h5('tau_V__Trg')
    # tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    ###########################################################################
    ###########################################################################
    ###########################################################################
    for iT, tSF in enumerate(tSF__T):
        x = 100. * x_Y__Tg[iT]
        y = tau_V_neb__g - tau_V__Tg[iT]
        #z = EW_Ha__g
        #mask = ~(x.mask | y.mask | z.mask | (EW_Ha__g < 1))
        mask = ~(x.mask | y.mask)
        xm = x[mask]
        ym = y[mask]
        #zm = z[mask]
        xlabel = r'$x_Y$'
        ylabel = r'$\tau_V^{neb} - \tau_V^\star$'
        #zlabel = r'$W_{H\alpha}$ [$\AA$]'
        #fname = 'xY_dtauV_EWHa_age_%sMyr.png' % str(tSF / 1.e6)
        fname = 'xY_dtauV_age_%sMyr.png' % str(tSF / 1.e6)
        f = plt.figure()
        f.set_size_inches(10,10)
        ax = f.gca()
        #sc = ax.scatter(xm, ym, c = zm, cmap = 'winter_r', marker = 'o', s = 5., edgecolor = 'none', vmax = 20, vmin = 3)
        sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 5., edgecolor = 'none')
        nBox = 20
        dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
        aux         = calcRunningStats(xm, ym, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
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
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        xlim = [xOkMin * 100., 80]
        ylim = [-1, 3]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(2.5))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.125))
        ax.grid(which = 'major')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # cb = f.colorbar(sc)
        # cb.set_label(zlabel)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        txt = r'%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
        f.suptitle(txt, fontsize = 14)
        f.savefig(fname)
        plt.close(f)
 
        x = 100. * x_Y__Tg[iT]
        y = tau_V_neb__g/tau_V__Tg[iT]
        z = EW_Ha__g
        #mask = ~(x.mask | y.mask | z.mask | (EW_Ha__g < 1))
        mask = ~(x.mask | y.mask)
        xm = x[mask]
        ym = y[mask]
        #zm = z[mask]
        xlabel = r'$X_Y$'
        ylabel = r'$\tau_V^{neb} / \tau_V^\star$'
        #zlabel = r'$W_{H\alpha}$ [$\AA$]'
        #fname = 'xY_tauVratio_EWHa_age_%sMyr.png' % str(tSF / 1.e6)
        fname = 'xY_tauVratio_age_%sMyr.png' % str(tSF / 1.e6)
        f = plt.figure()
        f.set_size_inches(10,10)
        ax = f.gca()
        #sc = ax.scatter(xm, ym, c = zm, cmap = 'winter_r', marker = 'o', s = 5., edgecolor = 'none', vmax = 20, vmin = 3)
        sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 5., edgecolor = 'none')
        nBox = 20
        dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
        aux         = calcRunningStats(xm, ym, dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
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
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        xlim = [xOkMin * 100., 80]
        ylim = [0, 5]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(2.5))
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.25))
        ax.grid(which = 'major')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # cb = f.colorbar(sc)
        # cb.set_label(zlabel)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        txt = r'%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
        f.suptitle(txt, fontsize = 14)
        f.savefig(fname)
        plt.close(f)
 
 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 #        x = tau_V__Tg[iT]
 #        y = tau_V_neb__g
 #        z = dist_zone__g
 #        mask = ~(x.mask | y.mask)
 #        xm = x[mask]
 #        ym = y[mask]
 #        zm = z[mask]
 #        xlabel = r'$\tau_V^\star$'
 #        ylabel = r'$\tau_V^{neb}$'
 #        zlabel = 'pixel distance (HLR)'
 #        fname = 'tauV_tauVNeb_pixDistHLR_age_%sMyr.png' % str(tSF / 1.e6)
 #        plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, 
 #                         [0, 2.5], [0, 2.5], fname, [0, 2.0], tSF)
 #         
 #        x = tau_V__Tg[iT]
 #        y = tau_V_neb__g
 #        z = np.ma.log10(L_int_Ha__g)
 #        mask = ~(x.mask | y.mask | z.mask)
 #        xm = x[mask]
 #        ym = y[mask]
 #        zm = z[mask]
 #        xlabel = r'$\tau_V^\star$'
 #        ylabel = r'$\tau_V^{neb}$'
 #        zlabel = r'$\log\ L_{H\alpha}\ [L_\odot]$'
 #        fname = 'tauV_tauVNeb_LHa_age_%sMyr.png' % str(tSF / 1.e6)
 #        plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, 
 #                         [0, 2.5], [0, 2.5], fname, [4.,6.], tSF)
 # 
 #        x = np.ma.log10(SFR_Ha__g)
 #        y1 = tau_V__Tg[iT]
 #        y2 = tau_V_neb__g
 #        mask = ~(x.mask | y1.mask | y2.mask)
 #        xm = x[mask]
 #        y1m = y1[mask]
 #        y2m = y2[mask]    
 #        xlabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'
 #        ylabel = r'$\tau_V$'
 #        fname = 'SFRNeb_tauV_tauVNeb_age_%sMyr.png' % str(tSF / 1.e6)
 #        f = plt.figure()
 #        f.set_size_inches(10,10)
 #        ax = f.gca()
 #        binsy = np.linspace(0., 4., 51)
 #        binsx = np.linspace(xm.min(), xm.max(), 51)
 #        ax.scatter(xm, y1m, c = 'b', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^\star$')
 #        #ax.hist2d(xm, y1m, bins = 50, cmap = 'Greys')
 #        ax.scatter(xm, y2m, c = 'r', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^{neb}$')
 #        #ax.hist2d(xm, y2m, bins = 50, cmap = 'Blues')
 #        #c = [ 'b', 'g', 'r', 'c', 'm', 'y' ]
 #        ax.set_ylim(0., 4.)
 #        ax.set_xlabel(xlabel)
 #        ax.set_ylabel(ylabel)
 #        ax.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
 #        ax.legend(fontsize = 14, frameon = False, loc = 'upper right')
 #        plt.setp(ax.get_xticklabels(), visible = False)
 #        plt.subplots_adjust(wspace=0, hspace=0)
 #        f.savefig(fname)
 #        plt.close(f)
 #         
 #        x = np.ma.log10(SFRSD_Ha__g * 1e6)
 #        y1 = tau_V__Tg[iT]
 #        y2 = tau_V_neb__g
 #        mask = ~(x.mask | y1.mask | y2.mask)
 #        xm = x[mask]
 #        y1m = y1[mask]
 #        y2m = y2[mask]    
 #        xlabel = r'$\log\ \Sigma_{SFR}^{neb}\ $[M${}_\odot\ $yr${}^{-1}\ $kpc${}^{-1}]$'
 #        ylabel = r'$\tau_V$'
 #        fname = 'SFRSDNeb_tauV_tauVNeb_age_%sMyr.png' % str(tSF / 1.e6)
 #        f = plt.figure()
 #        f.set_size_inches(10,10)
 #        ax = f.gca()
 #        binsy = np.linspace(0., 4., 51)
 #        binsx = np.linspace(min(xm), max(xm), 51)
 #        ax.scatter(xm, y1m, c = 'b', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^\star$')
 #        #ax.hist2d(xm, y1m, bins = 50, cmap = 'Greys')
 #        ax.scatter(xm, y2m, c = 'r', marker = 'o', s = 5., edgecolor = 'none', alpha = 0.4, label = r'$\tau_V^{neb}$')
 #        #ax.hist2d(xm, y2m, bins = 50, cmap = 'Blues')
 #        #c = [ 'b', 'g', 'r', 'c', 'm', 'y' ]
 #        ax.set_ylim(0, 4.)
 #        ax.set_xlabel(xlabel)
 #        ax.set_ylabel(ylabel)
 #        ax.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
 #        ax.legend(fontsize = 14, frameon = False, loc = 'upper right')
 #        plt.setp(ax.get_xticklabels(), visible = False)
 #        plt.subplots_adjust(wspace=0, hspace=0)
 #        f.savefig(fname)
 #        plt.close(f)
 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
