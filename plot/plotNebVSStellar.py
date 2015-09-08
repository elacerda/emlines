#!/usr/bin/python
#
# Lacerda@Granada - 20/Apr/2014
#
import sys
import numpy as np
import CALIFAUtils
import matplotlib as mpl
from matplotlib import pyplot as plt
from CALIFAUtils.plots import plot_zbins, plot_text_ax

debug = True
mask_radius = True
RNuc = 0.5

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
        
    H = CALIFAUtils.H5SFRData(h5file)
    
    #################################################################################
    #################################################################################
    #################################################################################
    
    iT = 11
    iU = -1
    tSF = H.tSF__T[iT]
    tZ = H.tZ__U[iU]
    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     f = plt.figure()
#     f.set_size_inches(10, 8)
#     f.set_dpi(100)
#     xk, xv = H.get_plot_dict(iT, iU, key = 'alogZmass')
#     yk, yv = H.get_plot_dict(iT, iU, key = 'logO3N2M13')
#     fnamepref = '%s_%s' % (xk, yk)
#     if mask_radius is True:
#         xv['v'][~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
#         xv['v'][~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
#         suptitle = r'NGals:%d  R > %.1fHLR  $x_Y$(min):%.0f%% ' % (H.N_gals, RNuc, H.xOkMin * 100.)
#         filename = '%s_maskradius_%.2fGyrs.png' % (fnamepref, tZ / 1e9)
#     else:
#         suptitle = r'NGals:%d  $x_Y$(min):%.0f%% ' % (H.N_gals, H.xOkMin * 100.)
#         filename = '%s_%.2fGyrs.png' % (fnamepref, tZ / 1e9)
#     plot_zbins(
#         f = f, 
#         ax = f.gca(),
#         debug = True,
#         x = xv['v'],
#         y = yv['v'],
#         xlabel = xv['label'],
#         ylabel = yv['label'],
#         xlim = xv['lim'],
#         ylim = yv['lim'],
#         running_stats = True,
#         rs_gaussian_smooth = True,
#         rs_percentiles = True,
#         rs_gs_fwhm = 8,
#         rs_frac_box = 80,
#         kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
#         kwargs_scatter = dict(c = 'b', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
#         kwargs_legend = dict(loc = 'best'),
#     )
#     f.suptitle(suptitle, fontsize = 12)
#     f.savefig(filename)
# 
#     f = plt.figure()
#     f.set_size_inches(10, 8)
#     f.set_dpi(100)
#     xk, xv = H.get_plot_dict(iT, iU, key = 'alogZmassR')
#     yk, yv = H.get_plot_dict(iT, iU, key = 'logO3N2M13R')
#     fnamepref = '%s_%s' % (xk, yk)
#     if mask_radius is True:
#         xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
#         xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
#         suptitle = r'NGals:%d  R > %.1fHLR  $x_Y$(min):%.0f%% ' % (H.N_gals, RNuc, H.xOkMin * 100.)
#         filename = '%s_maskradius_%.2fGyrs.png' % (fnamepref, tZ / 1e9)
#     else:
#         suptitle = r'NGals:%d  $x_Y$(min):%.0f%% ' % (H.N_gals, H.xOkMin * 100.)
#         filename = '%s_%.2fGyrs.png' % (fnamepref, tZ / 1e9)
#     plot_zbins(
#         f = f, 
#         ax = f.gca(),
#         debug = True,
#         x = xv['v'],
#         y = yv['v'],
#         xlabel = xv['label'],
#         ylabel = yv['label'],
#         xlim = xv['lim'],
#         ylim = yv['lim'],
#         z = H.Rtoplot(xv['v'].shape).flatten(),
#         zlabel = r'R (HLR)',
#         zmask = None,
#         zlim = [RNuc, 2],
#         running_stats = True,
#         rs_gaussian_smooth = True,
#         rs_percentiles = True,
#         rs_gs_fwhm = 8,
#         rs_frac_box = 80,
#         kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
#         kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
#         kwargs_legend = dict(loc = 'best'),
#         cb = True, 
#     )
#     f.suptitle(suptitle, fontsize = 12)
#     f.savefig(filename)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f, axArr = plt.subplots(1, 2)
    f.set_size_inches(15, 10)
    f.set_dpi(100)
    xk, xv = H.get_plot_dict(iT, iU, key = 'tauV')
    yk, yv = H.get_plot_dict(iT, iU, key = 'tauVNeb')
    fnamepref = '%s_%s' % (xk, yk)        
    if mask_radius is True:
        xv['v'][~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        xv['v'][~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $t_{SF}$:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        filename = '%s_maskradius_%.2fMyrs.png' % (fnamepref, tSF / 1e6)
    else:
        suptitle = r'NGals:%d  $t_{SF}$:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6)
    kw = plot_zbins(
        f = f, 
        ax = axArr[0],
        return_kwargs = True,
        debug = True,
        x = xv['v'],
        y = yv['v'],
        xlabel = xv['label'],
        ylabel = yv['label'],
        xlim = xv['lim'],
        ylim = yv['lim'],
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [0, 2],
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        ols = True,
        kwargs_ols = dict(c = 'r', pos_x = 0.99, pos_y = 0.02, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '-', lw = 2, label = 'OLS'),
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
        legend = False,
        cb = False,
    )
    ax_kw = kw['ax']
    ax_kw.plot(ax_kw.get_xlim(), (1/0.44) * np.asarray(ax_kw.get_xlim()), c = 'b', ls = '-.', lw = 3, label = r'$\tau_V^\star\ =\ 0.44\ \tau_V^{neb}$')
    ax_kw.legend(**dict(loc = 'upper left', fontsize = 14))

    xk, xv = H.get_plot_dict(iT, iU, key = 'xY')
    yk, yv = H.get_plot_dict(iT, iU, key = 'tauVRatio')
    fnamepref = '%s_%s' % (xk, yk)
    if mask_radius is True:
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
    plot_zbins(
        f = f, 
        ax = axArr[1],
        debug = True,
        x = xv['v'],
        y = yv['v'],
        xlabel = xv['label'],
        ylabel = yv['label'],
        xlim = xv['lim'],
        ylim = yv['lim'],
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        x_major_locator = xv['majloc'],
        x_minor_locator = xv['minloc'],
        y_major_locator = yv['majloc'],
        y_minor_locator = yv['minloc'],
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_percentiles = True,
        rs_gs_fwhm = 8,
        rs_frac_box = 80,
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
        kwargs_legend = dict(loc = 'best', fontsize = 14),
        cb = True,
    )
    f.suptitle(suptitle, fontsize = 14)
    f.savefig(filename)
         
    f, axArr = plt.subplots(1, 2)
    f.set_size_inches(15, 10)
    f.set_dpi(100)
    xk, xv = H.get_plot_dict(iT, iU, key = 'tauVR')
    yk, yv = H.get_plot_dict(iT, iU, key = 'tauVNebR')
    fnamepref = '%s_%s' % (xk, yk)        
    if mask_radius is True:
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $t_{SF}$:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        filename = '%s_maskradius_%.2fMyrs.png' % (fnamepref, tSF / 1e6)
    else:
        suptitle = r'NGals:%d  $t_{SF}$:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6)
    kw = plot_zbins(
        f = f, 
        ax = axArr[0],
        return_kwargs = True,
        debug = True,
        x = xv['v'],
        y = yv['v'],
        xlabel = xv['label'],
        ylabel = yv['label'],
        z = H.Rtoplot(shape = xv['v'].shape).flatten(),
        zlim = [0, 2],
        zlabel = r'R (HLR)',
        zmask = None, 
        xlim = xv['lim'],
        ylim = yv['lim'],
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        ols = True,
        kwargs_ols = dict(c = 'r', pos_x = 0.99, pos_y = 0.02, fs = 14, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '-', lw = 2, label = 'OLS'),
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
        legend = False,
        cb = False,
    )
    ax_kw = kw['ax']
    ax_kw.plot(ax_kw.get_xlim(), (1/0.44) * np.asarray(ax_kw.get_xlim()), c = 'b', ls = '-.', lw = 3, label = r'$\tau_V^\star\ =\ 0.44\ \tau_V^{neb}$')
    ax_kw.legend(**dict(loc = 'upper left', fontsize = 14))

    xk, xv = H.get_plot_dict(iT, iU, key = 'xYR')
    yk, yv = H.get_plot_dict(iT, iU, key = 'tauVRatioR')
    fnamepref = '%s_%s' % (xk, yk)
    if mask_radius is True:
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
        xv['v'][~(H.RbinCenter__r > RNuc)] = np.ma.masked
    plot_zbins(
        f = f, 
        ax = axArr[1],
        debug = True,
        x = xv['v'],
        y = yv['v'],
        xlabel = xv['label'],
        ylabel = yv['label'],
        xlim = xv['lim'],
        ylim = yv['lim'],
        z = H.Rtoplot(xv['v'].shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        x_major_locator = xv['majloc'],
        x_minor_locator = xv['minloc'],
        y_major_locator = yv['majloc'],
        y_minor_locator = yv['minloc'],
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_percentiles = True,
        rs_gs_fwhm = 8,
        rs_frac_box = 80,
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
        kwargs_legend = dict(loc = 'best'),
        cb = True, 
    )
    f.suptitle(suptitle, fontsize = 14)
    f.savefig(filename)