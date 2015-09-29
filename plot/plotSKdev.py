#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import CALIFAUtils as C
from CALIFAUtils.scripts import OLS_bisector, calc_running_stats
from CALIFAUtils.plots import plot_text_ax

debug = True

mpl.rcParams['font.size'] = 17
mpl.rcParams['axes.labelsize'] = 17
mpl.rcParams['axes.titlesize'] = 19
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
        
    H = C.H5SFRData(h5file)
    
    iT_values = [11, 17]
    
    #################################################################################
    #################################################################################
    #################################################################################
    
    updates = {
        'alogZmassR' : dict(lim = [-1.5, 0.5]),
        'atfluxR' : dict(lim = [7, 10]),
        'logO3N2S06R' : dict(lim = [-0.4, 0.1]),
        'logO3N2M13R' : dict(lim = [8, 9]),
        'logMcorSDR' : dict(lim = [2.5, 4.5]),
        'xYR' : dict(lim = [0, 50]),
        'logWHaWHbR' : dict(lim = [-0., 0.8]),
        'tauVdiffR' : dict(lim = [-0.75, 1.5]),
        'tauVRatioR' : dict(lim = [0., 6.]),
        'baR' : dict(lim = [0, 1.]),
    }
    
    xkeys = [ 'RadDist', 'morfTypeR', 'atfluxR', 'alogZmassR', 'logO3N2S06R', 'logMcorSDR', 'xYR', 'logWHaWHbR', 'tauVdiffR', 'tauVRatioR', 'baR' ]

    for iT in iT_values:
        tSF = H.tSF__T[iT]
        
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'logtauVR')
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'alogSFRSDR')
        xname = H.get_plot_dict(iT = iT, iU = -1)[xk[:-1]]['legendname']
        yname = H.get_plot_dict(iT = iT, iU = -1)[yk[1:-1]]['legendname']
        xm, ym = C.ma_mask_xyz(xv['v'], yv['v'])
        a, b, sig_a, sig_b = OLS_bisector(xm, ym)
        R = ym - (a * xm + b)
        Yrms = R.std()
        Yrms_str = r' : $y_{rms}$:%.2f' % Yrms
        if b > 0:
           txt_y = r'$y_{OLS}$ = %.2f %s + %.2f %s' % (yname, a, xname, b, Yrms_str)
        else:
           txt_y = r'$y_{OLS}\ \equiv$ %s = %.2f %s - %.2f %s' % (yname, a, xname, b * -1., Yrms_str)
        C.debug_var(debug, txt_y = txt_y)
        for xk in xkeys:
            fnamepref = 'SKdevOLS_%s' % xk
            if xk != 'RadDist':
                _, xv = H.get_plot_dict(iT = iT, iU = -1, key = xk)
                xv.update(updates.get(xk, {}))
            else:
                xv = dict(v = (np.ones_like(xv['v']).T * H.RbinCenter__r).T,
                          lim = [0, 2],
                          label = r'R [HLR]',
                          legendname = 'R',
                          minloc = 0.05,
                          majloc = 0.4
                )
            kwargs_figure = dict(figsize = (10, 8), dpi = 100)
            f = plt.figure(**kwargs_figure)
            ax = f.gca()
            plot_text_ax(ax, txt_y, 0.99, 0.01, 14, 'bottom', 'right', color = 'k')
            kw = C.plot_zbins(
                return_kwargs = True,
                f = f,
                ax = ax,
                debug = debug,
                x = xv['v'],
                y = R,
                xlim = xv['lim'],
                #ylim = yv['lim'],
                ylimprc = [0, 100],
                xlabel = xv['label'],
                ylabel = r'Residual %s' % yname,
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                spearmanr = True,
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
            )
            kwargs_suptitle = dict(fontsize = 12)
            suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            f.suptitle(suptitle, **kwargs_suptitle)
            filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.debug_var(debug, filename = filename)
            f.savefig(filename)
            
        #########################
        #########################
        
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'logtauVR')
        #x2k, x2v = H.get_plot_dict(iT = iT, iU = -1, key = 'logMcorSDR')
        #xv['v'] = xv['v'] - x2v['v']
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'alogSFRSDR')
        xname = H.get_plot_dict(iT = iT, iU = -1)[xk[:-1]]['legendname']
        yname = H.get_plot_dict(iT = iT, iU = -1)[yk[1:-1]]['legendname']
        xm, ym = C.ma_mask_xyz(xv['v'], yv['v'])
        nBox = len(xm.compressed()) / 20.
        dxBox = (xm.max() - xm.min()) / (nBox - 1.)
        kwargs_rs = dict(dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
        C.debug_var(debug, kwargs_rs = kwargs_rs)
        xbinCenter, xMedian, xMean, xStd, yMedian, yMean, yStd, nInBin, xPrc, yPrc = calc_running_stats(xm, ym, **kwargs_rs)
        m = (np.isnan(xMedian) | np.isnan(yMedian))
        p = np.polyfit(xMedian[~m], yMedian[~m], 3)
        R = ym - np.polyval(p, xm)
        Yrms = R.std()
        txt_y = r'$y_{rms}$:%.2f' % Yrms
        C.debug_var(debug, txt_y = txt_y)
        for xk in xkeys:
            fnamepref = 'SKdevMedian_%s' % xk
            if xk != 'RadDist':
                _, xv = H.get_plot_dict(iT = iT, iU = -1, key = xk)
                xv.update(updates.get(xk, {}))
            else:
                xv = dict(v = (np.ones_like(xv['v']).T * H.RbinCenter__r).T,
                          lim = [0, 2],
                          label = r'R [HLR]',
                          legendname = 'R',
                          minloc = 0.05,
                          majloc = 0.4
                )
            kwargs_figure = dict(figsize = (10, 8), dpi = 100)
            f = plt.figure(**kwargs_figure)
            ax = f.gca()
            plot_text_ax(ax, txt_y, 0.99, 0.01, 14, 'bottom', 'right', color = 'k')
            kw = C.plot_zbins(
                return_kwargs = True,
                f = f,
                ax = ax,
                debug = debug,
                x = xv['v'],
                y = R,
                xlim = xv['lim'],
                #ylim = yv['lim'],
                ylimprc = [0, 100],
                xlabel = xv['label'],
                ylabel = r'Residual %s' % yname,
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                spearmanr = True,
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
            )
            kwargs_suptitle = dict(fontsize = 12)
            suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            f.suptitle(suptitle, **kwargs_suptitle)
            filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.debug_var(debug, filename = filename)
            f.savefig(filename)

    #########################
    #########################

    xkeys = [ 'zoneDistHLR', 'logZoneArea', 'morfType', 'atflux', 'alogZmass', 'logO3N2S06', 'logMcorSD', 'xY', 'logWHaWHb', 'tauVdiff', 'tauVRatio', 'ba' ]

    for iT in iT_values:
        tSF = H.tSF__T[iT]
        
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'logtauV')
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'logSFRSD')
        xname = xv['legendname']
        yname = yv['legendname']
        xm, ym = C.ma_mask_xyz(xv['v'], yv['v'])
        a, b, sig_a, sig_b = OLS_bisector(xm, ym)
        R = ym - (a * xm + b)
        Yrms = R.std()
        Yrms_str = r' : $y_{rms}$:%.2f' % Yrms
        if b > 0:
           txt_y = r'$y_{OLS}$ = %.2f %s + %.2f %s' % (yname, a, xname, b, Yrms_str)
        else:
           txt_y = r'$y_{OLS}\ \equiv$ %s = %.2f %s - %.2f %s' % (yname, a, xname, b * -1., Yrms_str)
        C.debug_var(debug, txt_y = txt_y)
        for xk in xkeys:
            fnamepref = 'SKdevOLS_%s' % xk
            _, xv = H.get_plot_dict(iT = iT, iU = -1, key = xk)
            kwargs_figure = dict(figsize = (10, 8), dpi = 100)
            f = plt.figure(**kwargs_figure)
            ax = f.gca()
            plot_text_ax(ax, txt_y, 0.99, 0.01, 14, 'bottom', 'right', color = 'k')
            #xv.update(updates.get(xk, {}))
            kw = C.plot_zbins(
                return_kwargs = True,
                f = f,
                ax = ax,
                debug = debug,
                x = xv['v'],
                y = R,
                xlim = xv['lim'],
                #ylim = yv['lim'],
                ylimprc = [0, 100],
                xlabel = xv['label'],
                ylabel = r'Residual %s' % yname,
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                spearmanr = True,
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
            )
            kwargs_suptitle = dict(fontsize = 12)
            suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            f.suptitle(suptitle, **kwargs_suptitle)
            filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.debug_var(debug, filename = filename)
            f.savefig(filename)
            
        #########################
        #########################
        
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'logtauV')
        #x2k, x2v = H.get_plot_dict(iT = iT, iU = -1, key = 'logMcorSD')
        #xv['v'] = xv['v'] - x2v['v']
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'logSFRSD')
        xname = xv['legendname']
        yname = yv['legendname']
        xm, ym = C.ma_mask_xyz(xv['v'], yv['v'])
        nBox = len(xm.compressed()) / 20.
        dxBox = (xm.max() - xm.min()) / (nBox - 1.)
        kwargs_rs = dict(dxBox = dxBox, xbinIni = xm.min(), xbinFin = xm.max(), xbinStep = dxBox)
        C.debug_var(debug, kwargs_rs = kwargs_rs)
        xbinCenter, xMedian, xMean, xStd, yMedian, yMean, yStd, nInBin, xPrc, yPrc = calc_running_stats(xm, ym, **kwargs_rs)
        m = (np.isnan(xMedian) | np.isnan(yMedian))
        p = np.polyfit(xMedian[~m], yMedian[~m], 3)
        R = ym - np.polyval(p, xm)
        Yrms = R.std()
        txt_y = r'$y_{rms}$:%.2f' % Yrms
        C.debug_var(debug, txt_y = txt_y)
        for xk in xkeys:
            fnamepref = 'SKdevMedian_%s' % xk
            _, xv = H.get_plot_dict(iT = iT, iU = -1, key = xk)
            kwargs_figure = dict(figsize = (10, 8), dpi = 100)
            f = plt.figure(**kwargs_figure)
            ax = f.gca()
            plot_text_ax(ax, txt_y, 0.99, 0.01, 14, 'bottom', 'right', color = 'k')
            #xv.update(updates.get(xk, {}))
            kw = C.plot_zbins(
                return_kwargs = True,
                f = f,
                ax = ax,
                debug = debug,
                x = xv['v'],
                y = R,
                xlim = xv['lim'],
                #ylim = yv['lim'],
                ylimprc = [0, 100],
                xlabel = xv['label'],
                ylabel = r'Residual %s' % yname,
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                spearmanr = True,
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
            )
            kwargs_suptitle = dict(fontsize = 12)
            suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
            f.suptitle(suptitle, **kwargs_suptitle)
            filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.debug_var(debug, filename = filename)
            f.savefig(filename)
        