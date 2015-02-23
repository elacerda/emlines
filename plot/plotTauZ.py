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
from plot_aux import plotOLSbisectorAxis, plot_text_ax, \
                     calcRunningStats, gaussSmooth_YofX, \
                     OLS_bisector


mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

iT_values = [ 11, 17 ]

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    dtSF__T = H.tSF__T
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    #tauVtoAV = np.log10(np.exp(1)) / 0.4
    #zticks3 = [ -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2 ]
    
    alogZmass_edges_radius = [
        #-0.35755388,
        -0.35, 
        #-0.17530025,
        -0.2,
        #-0.06933334,
        -0.1,
        #0.01808038,
        0.,
    ]
    alogZmass_edges_zones = [
        #-0.3582655 , 
        -0.45, 
        #-0.12654235, 
        -0.2, 
        #0.01805482, 
        0., 
        #0.15001057,
        0.1,
    ]
    alogZmass_minmax_radius = [-0.75, 0.1]
    alogZmass_minmax_zones = [-0.95, 0.2]

    logZneb_edges_radius = [
        #-0.12385368,
        -0.12, 
        #-0.08062921,
        -0.08, 
        #-0.05185106,
        -0.05, 
        #-0.02014183,
        -0.02,
    ]
    logZneb_edges_zones = [
        #-0.10732491,
        -0.1, 
        #-0.0645055 ,
        -0.06,
        #-0.03323017,
        -0.03,
        #0.00041659,
        0.,
    ]
    
    logZneb_minmax_zones = [-0.2, 0.05]
    logZneb_minmax_radius = [-0.2, 0.05]        
    
    zticklabelcolors = [(1, 0, 0), (0, 1, 0), (0, 1, .5), (0, 0.5, 1), (0, 0, 1)]
    for i in iT_values:
        iT = i
        tSF = H.tSF__T[iT]
        suptitle_txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)    
        #######################################################################
        #######################################################################
        ##### Radius ##########################################################
        #######################################################################
        #######################################################################
        x = H.tau_V__Trg[iT]
        y = H.tau_V_neb__rg
        z = H.alogZ_mass__Urg[-1]
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        xlabel = r'$\tau_V^\star(R)$'
        ylabel = r'$\tau_V^{neb}(R)$'
        zlabel = r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)
        fname = 'tauVR_tauVNebR_alogZmassR_age_%sMyr.png' % str(tSF / 1.e6)
        xlim = [0, 1.5]
        ylim = [0, 2.5]
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', 
                        vmin = alogZmass_minmax_radius[0],
                        vmax = alogZmass_minmax_radius[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
        a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
        zticks_mask = [ 
            (zm <= alogZmass_edges_radius[0]), 
            (zm > alogZmass_edges_radius[0]) & (zm <= alogZmass_edges_radius[1]),
            (zm > alogZmass_edges_radius[1]) & (zm <= alogZmass_edges_radius[2]),
            (zm > alogZmass_edges_radius[2]) & (zm <= alogZmass_edges_radius[3]),
            (zm > alogZmass_edges_radius[3]),
        ]
        z_labels = [
            'Z <= %.2f' % alogZmass_edges_radius[0],
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[0],alogZmass_edges_radius[1]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[1],alogZmass_edges_radius[2]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[2],alogZmass_edges_radius[3]),
            'Z > %.2f' % alogZmass_edges_radius[3],
        ] 
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = '0.5', ls = '--', lw = 0.5)
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
        #####################
        ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
        ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
        for i, msk in enumerate(zticks_mask):
            y_pos = 0.07 + 0.05 * i
            X = xm[msk]
            Y = ym[msk]
            a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                         0.98, y_pos, 12, 
                                                         zticklabelcolors[i], 
                                                         rms = True, 
                                                         label = z_labels[i],
                                                         text = False)
        ax.legend(loc = 'upper left', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)
        f.savefig(fname)
        plt.close(f)

        x = H.tau_V__Trg[iT]
        y = H.tau_V_neb__rg
        z = H.logZ_neb_S06__rg
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        xlabel = r'$\tau_V^\star(R)$'
        ylabel = r'$\tau_V^{neb}(R)$'
        zlabel = r'$\log\ Z_{neb}(R)$ [$Z_\odot$]'
        fname = 'tauVR_tauVNebR_logZnebR_age_%sMyr.png' % str(tSF / 1.e6)
        xlim = [0, 1.5]
        ylim = [0, 2.5]
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r',
                        vmin = logZneb_minmax_radius[0],
                        vmax = logZneb_minmax_radius[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
        a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
        zticks_mask = [ 
            (zm <= logZneb_edges_radius[0]), 
            (zm > logZneb_edges_radius[0]) & (zm <= logZneb_edges_radius[1]),
            (zm > logZneb_edges_radius[1]) & (zm <= logZneb_edges_radius[2]),
            (zm > logZneb_edges_radius[2]) & (zm <= logZneb_edges_radius[3]),
            (zm > logZneb_edges_radius[3]),
        ] 
        z_labels = [
            'Z <= %.2f' % logZneb_edges_radius[0],
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[0],logZneb_edges_radius[1]),
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[1],logZneb_edges_radius[2]),
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[2],logZneb_edges_radius[3]),
            'Z > %.2f' % logZneb_edges_radius[3],
        ] 
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
        # cb = f.colorbar(sc, ticks=zticks)
        # cb.set_label(zlabel)
        # cb.ax.set_yticklabels(zticklabels)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = '0.5', ls = '--', lw = 0.5)
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
        #####################
        ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
        ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
        for i, msk in enumerate(zticks_mask):
            y_pos = 0.07 + 0.05 * i
            X = xm[msk]
            Y = ym[msk]
            a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                         0.98, y_pos, 12, 
                                                         zticklabelcolors[i], 
                                                         rms = True, 
                                                         label = z_labels[i],
                                                         text = False)
        ax.legend(loc = 'upper left', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)     #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        f.savefig(fname)
        plt.close(f)
        
        x = H.tau_V__Trg[iT]
        y = H.tau_V_neb__rg
        z = H.alogZ_mass__Urg[-1]
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        a, b, sigma_a, sigma_b = OLS_bisector(xm,ym)
        R = ym - (a * xm + b)
        zticks_mask = [ 
            (zm <= alogZmass_edges_radius[0]), 
            (zm > alogZmass_edges_radius[0]) & (zm <= alogZmass_edges_radius[1]),
            (zm > alogZmass_edges_radius[1]) & (zm <= alogZmass_edges_radius[2]),
            (zm > alogZmass_edges_radius[2]) & (zm <= alogZmass_edges_radius[3]),
            (zm > alogZmass_edges_radius[3]),
        ]
        z_labels = [
            'Z <= %.2f' % alogZmass_edges_radius[0],
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[0],alogZmass_edges_radius[1]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[1],alogZmass_edges_radius[2]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_radius[2],alogZmass_edges_radius[3]),
            'Z > %.2f' % alogZmass_edges_radius[3],
        ] 
        xlabel = r'$\tau_V^\star(R)$'
        ylabel = r'$\tau_V^{neb}(R)\ -\ y_{OLS}(R)$'
        zlabel = r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)
        fname = 'tauVR_ResTauVOLSR_alogZmassR_%sMyr.png' % str(tSF / 1.e6) 
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        xlim = [0, 1.5]
        ylim = [-2, 2]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if b > 0:
            txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ + %.2f' %  (a, b)
        else:
            txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ - %.2f' %  (a, b * -1.)
        plot_text_ax(ax, txt, 0.98, 0.02, 14, 'bottom', 'right', color = 'k')
        sc = ax.scatter(xm, R, c = zm, cmap = 'spectral_r',
                        vmin = alogZmass_minmax_radius[0],
                        vmax = alogZmass_minmax_radius[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.3)    
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        for i, msk in enumerate(zticks_mask):
            X = xm[msk]
            Y = R[msk]
            if len(X) > 3 and len(Y) > 3:
                nBox = 50
                dxBox       = (X.max() - X.min()) / (nBox - 1.)
                aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
                FWHM = 0.4
                xM = np.ma.masked_array(xMedian)
                yM = np.ma.masked_array(yMedian)
                mask = np.isnan(xM) | np.isnan(yM) 
                Xs, Ys = gaussSmooth_YofX(xM[~mask], yM[~mask], FWHM)
                ax.plot(Xs, Ys, c = zticklabelcolors[i], lw = 2, label = '%s' % z_labels[i])
        ax.legend(loc = 'upper right', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(1.00))
        ax.yaxis.set_minor_locator(MultipleLocator(0.20))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)
        #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        f.savefig(fname)
        plt.close(f)

        x = H.tau_V__Trg[iT]
        y = H.tau_V_neb__rg
        z = H.logZ_neb_S06__rg
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        a, b, sigma_a, sigma_b = OLS_bisector(xm,ym)
        R = ym - (a * xm + b)
        zticks_mask = [ 
            (zm <= logZneb_edges_radius[0]), 
            (zm > logZneb_edges_radius[0]) & (zm <= logZneb_edges_radius[1]),
            (zm > logZneb_edges_radius[1]) & (zm <= logZneb_edges_radius[2]),
            (zm > logZneb_edges_radius[2]) & (zm <= logZneb_edges_radius[3]),
            (zm > logZneb_edges_radius[3]),
        ]
        z_labels = [
            'Z <= %.2f' % logZneb_edges_radius[0],
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[0],logZneb_edges_radius[1]),
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[1],logZneb_edges_radius[2]),
            '%.2f < Z <= %.2f' % (logZneb_edges_radius[2],logZneb_edges_radius[3]),
            'Z > %.2f' % logZneb_edges_radius[3],
        ] 
        xlabel = r'$\tau_V^\star(R)$'
        ylabel = r'$\tau_V^{neb}(R)\ -\ y_{OLS}(R)$'
        zlabel = r'$\log\ Z_{neb}(R)$ [$Z_\odot$]'
        fname = 'tauVR_ResTauVOLSR_logZnebR_%sMyr.png' % str(tSF / 1.e6) 
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        xlim = [0, 1.5]
        ylim = [-2, 2]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if b > 0:
            txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ + %.2f' %  (a, b)
        else:
            txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ - %.2f' %  (a, b * -1.)
        plot_text_ax(ax, txt, 0.98, 0.02, 14, 'bottom', 'right', color = 'k')
        sc = ax.scatter(xm, R, c = zm, cmap = 'spectral_r',
                        vmin = alogZmass_minmax_radius[0],
                        vmax = alogZmass_minmax_radius[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.3)    
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        for i, msk in enumerate(zticks_mask):
            X = xm[msk]
            Y = R[msk]
            if len(X) > 3 and len(Y) > 3:
                nBox = 50
                dxBox       = (X.max() - X.min()) / (nBox - 1.)
                aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
                FWHM = 0.4
                xM = np.ma.masked_array(xMedian)
                yM = np.ma.masked_array(yMedian)
                mask = np.isnan(xM) | np.isnan(yM) 
                Xs, Ys = gaussSmooth_YofX(xM[~mask], yM[~mask], FWHM)
                ax.plot(Xs, Ys, c = zticklabelcolors[i], lw = 2, label = '%s' % z_labels[i])
        ax.legend(loc = 'upper right', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(1.00))
        ax.yaxis.set_minor_locator(MultipleLocator(0.20))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)
        #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        f.savefig(fname)
        plt.close(f)


        #######################################################################
        #######################################################################
        ##### Zones ###########################################################
        #######################################################################
        #######################################################################
        x = H.tau_V__Tg[iT]
        y = H.tau_V_neb__g
        z = H.alogZ_mass__Ug[-1]
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        xlabel = r'$\tau_V^\star$'
        ylabel = r'$\tau_V^{neb}$'
        zlabel = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)
        fname = 'tauV_tauVNeb_alogZmass_age_%sMyr.png' % str(tSF / 1.e6)
        xlim = [0, 1.5]
        ylim = [0, 2.5]
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r',
                        vmin = alogZmass_minmax_zones[0],
                        vmax = alogZmass_minmax_zones[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
        a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
        # cb = f.colorbar(sc, ticks=zticks)
        # cb.set_label(zlabel)
        # cb.ax.set_yticklabels(zticklabels)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        zticks_mask = [ 
            (zm <= alogZmass_edges_zones[0]), 
            (zm > alogZmass_edges_zones[0]) & (zm <= alogZmass_edges_zones[1]),
            (zm > alogZmass_edges_zones[1]) & (zm <= alogZmass_edges_zones[2]),
            (zm > alogZmass_edges_zones[2]) & (zm <= alogZmass_edges_zones[3]),
            (zm > alogZmass_edges_zones[3]),
        ]
        
        z_labels = [
            'Z <= %.2f' % alogZmass_edges_zones[0],
            '%.2f < Z <= %.2f' % (alogZmass_edges_zones[0],alogZmass_edges_zones[1]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_zones[1],alogZmass_edges_zones[2]),
            '%.2f < Z <= %.2f' % (alogZmass_edges_zones[2],alogZmass_edges_zones[3]),
            'Z > %.2f' % alogZmass_edges_zones[3],
        ] 
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        #cb.ax.set_yticklabels(zticklabels)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = '0.5', ls = '--', lw = 0.5)
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
        #####################
        ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
        ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
        for i, msk in enumerate(zticks_mask):
            y_pos = 0.07 + 0.05 * i
            X = xm[msk]
            Y = ym[msk]
            a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                         0.98, y_pos, 12, 
                                                         zticklabelcolors[i], 
                                                         rms = True, 
                                                         label = z_labels[i],
                                                         text = False)
        ax.legend(loc = 'upper left', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)        
        f.savefig(fname)
        plt.close(f)

        x = H.tau_V__Tg[iT]
        y = H.tau_V_neb__g
        z = H.logZ_neb_S06__g
        mask = x.mask | y.mask | z.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        zm = np.ma.masked_array(z, mask = mask)
        xlabel = r'$\tau_V^\star$'
        ylabel = r'$\tau_V^{neb}$'
        zlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
        fname = 'tauV_tauVNeb_logZneb_age_%sMyr.png' % str(tSF / 1.e6)
        xlim = [0, 1.5]
        ylim = [0, 2.5]
        f = plt.figure()
        f.set_size_inches(10, 8)
        ax = f.gca()
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r',
                        vmin = logZneb_minmax_zones[0],
                        vmax = logZneb_minmax_zones[1], 
                        marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
        a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
        zticks_mask = [ 
            (zm <= logZneb_edges_zones[0]), 
            (zm > logZneb_edges_zones[0]) & (zm <= logZneb_edges_zones[1]),
            (zm > logZneb_edges_zones[1]) & (zm <= logZneb_edges_zones[2]),
            (zm > logZneb_edges_zones[2]) & (zm <= logZneb_edges_zones[3]),
            (zm > logZneb_edges_zones[3]),
        ] 
        z_labels = [
            'Z <= %.2f' % logZneb_edges_zones[0],
            '%.2f < Z <= %.2f' % (logZneb_edges_zones[0],logZneb_edges_zones[1]),
            '%.2f < Z <= %.2f' % (logZneb_edges_zones[1],logZneb_edges_zones[2]),
            '%.2f < Z <= %.2f' % (logZneb_edges_zones[2],logZneb_edges_zones[3]),
            'Z > %.2f' % logZneb_edges_zones[3],
        ] 
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
        # cb = f.colorbar(sc, ticks=zticks)
        # cb.set_label(zlabel)
        # cb.ax.set_yticklabels(zticklabels)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        cb = f.colorbar(sc)
        cb.set_label(zlabel)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(xm, Y, c = '0.5', ls = '--', lw = 0.5)
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
        #####################
        ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
        ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
        for i, msk in enumerate(zticks_mask):
            y_pos = 0.07 + 0.05 * i
            X = xm[msk]
            Y = ym[msk]
            a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                         0.98, y_pos, 12, 
                                                         zticklabelcolors[i], 
                                                         rms = True, 
                                                         label = z_labels[i],
                                                         text = False)
        ax.legend(loc = 'upper left', fontsize = 14)
        ax.xaxis.set_major_locator(MultipleLocator(0.25))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.grid(which = 'major')
        f.suptitle(suptitle_txt, fontsize = 14)
        f.savefig(fname)
        plt.close(f)
        




exit(1)
















iT = 11
tSF = H.tSF__T[iT]
x = np.ma.log10(H.tau_V__Trg[iT].flatten())
y = np.ma.log10(H.aSFRSD__Trg[11].flatten() * 1e6)
z = H.alogZ_mass__Urg[-1].flatten()
mask = x.mask | y.mask | z.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
zm = np.ma.masked_array(z, mask = mask)
xlabel = r'$\log\ \tau_V^\star(R)$'
ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
zlabel = r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)
fname = 'logtauVR_logaSFRSDR_alogZmassR_age_%sMyr.png' % str(tSF / 1.e6)
xlim = [-1.5, 0.5]
ylim = [-3, 1]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', 
                vmin = alogZmass_minmax_radius[0],
                vmax = alogZmass_minmax_radius[1], 
                marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
zticks_mask = [ 
    (zm <= alogZmass_edges_radius[0]), 
    (zm > alogZmass_edges_radius[0]) & (zm <= alogZmass_edges_radius[1]),
    (zm > alogZmass_edges_radius[1]) & (zm <= alogZmass_edges_radius[2]),
    (zm > alogZmass_edges_radius[2]) & (zm <= alogZmass_edges_radius[3]),
    (zm > alogZmass_edges_radius[3]),
]

z_labels = [
    'Z <= %.2f' % alogZmass_edges_radius[0],
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[0],alogZmass_edges_radius[1]),
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[1],alogZmass_edges_radius[2]),
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[2],alogZmass_edges_radius[3]),
    'Z > %.2f' % alogZmass_edges_radius[3],
] 
cb = f.colorbar(sc)
cb.set_label(zlabel)
#####################
for i, msk in enumerate(zticks_mask):
    y_pos = 0.07 + 0.05 * i
    X = xm[msk]
    Y = ym[msk]
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                 0.98, y_pos, 12, 
                                                 zticklabelcolors[i], 
                                                 rms = True, 
                                                 label = z_labels[i],
                                                 text = False)
ax.legend(loc = 'upper left', fontsize = 14)
ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.grid(which = 'major')
f.suptitle(suptitle_txt, fontsize = 14)
f.savefig(fname)
plt.close(f)

x = np.ma.log10(H.tau_V__Trg[iT].flatten())
y = np.ma.log10(H.aSFRSD__Trg[11].flatten() * 1e6)
z = H.logZ_neb_S06__rg.flatten()
mask = x.mask | y.mask | z.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
zm = np.ma.masked_array(z, mask = mask)
xlabel = r'$\log\ \tau_V^\star(R)$'
ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
zlabel = r'$\log\ Z_{neb}(R)$ [$Z_\odot$]'
fname = 'logtauVR_logaSFRSDR_logZNebR_age_%sMyr.png' % str(tSF / 1.e6)
xlim = [-1.5, 0.5]
ylim = [-3, 1]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r',
                vmin = logZneb_minmax_radius[0],
                vmax = logZneb_minmax_radius[1], 
                marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
zticks_mask = [ 
    (zm <= logZneb_edges_radius[0]), 
    (zm > logZneb_edges_radius[0]) & (zm <= logZneb_edges_radius[1]),
    (zm > logZneb_edges_radius[1]) & (zm <= logZneb_edges_radius[2]),
    (zm > logZneb_edges_radius[2]) & (zm <= logZneb_edges_radius[3]),
    (zm > logZneb_edges_radius[3]),
] 
z_labels = [
    'Z <= %.2f' % logZneb_edges_radius[0],
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[0],logZneb_edges_radius[1]),
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[1],logZneb_edges_radius[2]),
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[2],logZneb_edges_radius[3]),
    'Z > %.2f' % logZneb_edges_radius[3],
] 
cb = f.colorbar(sc)
cb.set_label(zlabel)
#####################
for i, msk in enumerate(zticks_mask):
    y_pos = 0.07 + 0.05 * i
    X = xm[msk]
    Y = ym[msk]
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                 0.98, y_pos, 12, 
                                                 zticklabelcolors[i], 
                                                 rms = True, 
                                                 label = z_labels[i],
                                                 text = False)
ax.legend(loc = 'upper left', fontsize = 14)
ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.grid(which = 'major')
f.suptitle(suptitle_txt, fontsize = 14)
f.savefig(fname)
plt.close(f)

x = np.ma.log10(H.tau_V_neb__rg.flatten())
y = np.ma.log10(H.aSFRSD__Trg[11].flatten() * 1e6)
z = H.alogZ_mass__Urg[-1].flatten()
mask = x.mask | y.mask | z.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
zm = np.ma.masked_array(z, mask = mask)
xlabel = r'$\log\ \tau_V^{neb}(R)$'
ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
zlabel = r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9)
fname = 'logtauVNebR_logaSFRSDR_alogZmassR_age_%sMyr.png' % str(tSF / 1.e6)
xlim = [-1.3, 1]
ylim = [-3, 1]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', 
                vmin = alogZmass_minmax_radius[0],
                vmax = alogZmass_minmax_radius[1], 
                marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
zticks_mask = [ 
    (zm <= alogZmass_edges_radius[0]), 
    (zm > alogZmass_edges_radius[0]) & (zm <= alogZmass_edges_radius[1]),
    (zm > alogZmass_edges_radius[1]) & (zm <= alogZmass_edges_radius[2]),
    (zm > alogZmass_edges_radius[2]) & (zm <= alogZmass_edges_radius[3]),
    (zm > alogZmass_edges_radius[3]),
]

z_labels = [
    'Z <= %.2f' % alogZmass_edges_radius[0],
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[0],alogZmass_edges_radius[1]),
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[1],alogZmass_edges_radius[2]),
    '%.2f < Z <= %.2f' % (alogZmass_edges_radius[2],alogZmass_edges_radius[3]),
    'Z > %.2f' % alogZmass_edges_radius[3],
] 
cb = f.colorbar(sc)
cb.set_label(zlabel)
#####################
for i, msk in enumerate(zticks_mask):
    y_pos = 0.07 + 0.05 * i
    X = xm[msk]
    Y = ym[msk]
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                 0.98, y_pos, 12, 
                                                 zticklabelcolors[i], 
                                                 rms = True, 
                                                 label = z_labels[i],
                                                 text = False)
ax.legend(loc = 'upper left', fontsize = 14)
ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.grid(which = 'major')
f.suptitle(suptitle_txt, fontsize = 14)
f.savefig(fname)
plt.close(f)

x = np.ma.log10(H.tau_V_neb__rg.flatten())
y = np.ma.log10(H.aSFRSD__Trg[11].flatten() * 1e6)
z = H.logZ_neb_S06__rg.flatten()
mask = x.mask | y.mask | z.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
zm = np.ma.masked_array(z, mask = mask)
xlabel = r'$\log\ \tau_V^{neb}(R)$'
ylabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
zlabel = r'$\log\ Z_{neb}(R)$ [$Z_\odot$]'
fname = 'logtauVNebR_logaSFRSDR_logZNebR_age_%sMyr.png' % str(tSF / 1.e6)
xlim = [-1.3, 1]
ylim = [-3, 1]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r',
                vmin = logZneb_minmax_radius[0],
                vmax = logZneb_minmax_radius[1], 
                marker = 'o', s = 15., edgecolor = 'none', alpha = 0.6)
a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
zticks_mask = [ 
    (zm <= logZneb_edges_radius[0]), 
    (zm > logZneb_edges_radius[0]) & (zm <= logZneb_edges_radius[1]),
    (zm > logZneb_edges_radius[1]) & (zm <= logZneb_edges_radius[2]),
    (zm > logZneb_edges_radius[2]) & (zm <= logZneb_edges_radius[3]),
    (zm > logZneb_edges_radius[3]),
] 
z_labels = [
    'Z <= %.2f' % logZneb_edges_radius[0],
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[0],logZneb_edges_radius[1]),
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[1],logZneb_edges_radius[2]),
    '%.2f < Z <= %.2f' % (logZneb_edges_radius[2],logZneb_edges_radius[3]),
    'Z > %.2f' % logZneb_edges_radius[3],
] 
cb = f.colorbar(sc)
cb.set_label(zlabel)
#####################
for i, msk in enumerate(zticks_mask):
    y_pos = 0.07 + 0.05 * i
    X = xm[msk]
    Y = ym[msk]
    a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
                                                 0.98, y_pos, 12, 
                                                 zticklabelcolors[i], 
                                                 rms = True, 
                                                 label = z_labels[i],
                                                 text = False)
ax.legend(loc = 'upper left', fontsize = 14)
ax.xaxis.set_major_locator(MultipleLocator(0.25))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.grid(which = 'major')
f.suptitle(suptitle_txt, fontsize = 14)
#plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
f.savefig(fname)
plt.close(f)
