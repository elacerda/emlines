#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
import CALIFAUtils as C 

mask_radius = True
#mask_radius = False
RNuc = 0.5

def iTGen(tSF__T, iT_values = [ 11, 17 ]):
    for iT in iT_values:
        yield iT, tSF__T[iT]

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE CALIFAID' % (sys.argv[0])
    exit(1)

H = C.H5SFRData(h5file)
tSF__T = H.tSF__T
iT_values = [ 11 ]

for iT, tSF in iTGen(H.tSF__T, iT_values):
    
    tau_V_norm_GAL__g = (H.tau_V__Trg[iT][10, :] + H.tau_V__Trg[iT][9, :] / 2.)
    tau_V_norm__rg = H.tau_V__Trg[iT] / tau_V_norm_GAL__g
    SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][10, :] + H.aSFRSD__Trg[iT][9, :] / 2.)
    aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g

    x1 = np.ma.copy(np.ma.log10(H.tau_V__Tg[iT]))
    y1 = np.ma.copy(np.ma.log10(H.SFRSD__Tg[iT] * 1e6))
    x2 = np.ma.copy(np.ma.log10(H.tau_V__Trg[iT]))
    y2 = np.ma.copy(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6)) 
    x3 = np.ma.copy(np.ma.log10(tau_V_norm__rg))
    y3 = np.ma.copy(np.ma.log10(aSFRSD_norm__rg)) 

    x1label = r'$\log\ \tau_V^{\star}$'
    y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x2label = r'$\log\ \tau_V^{\star}(R)$'
    y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
    y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'

    if mask_radius is True:
        x1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        y1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        x2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        x3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_mask_%sMyr.png' % str(tSF / 1.e6)
    else:
        suptitle = r'NGals:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_%sMyr.png' % str(tSF / 1.e6)

    f, axArr = plt.subplots(1, 3)
    f.set_dpi(100)
    f.set_size_inches(15, 8)
    C.plot_zbins(
        f = f,
        ax = axArr[0],
        debug = True,
        x = x1,
        y = y1,
        xlabel = x1label,
        ylabel = y1label,
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[1],
        debug = True,
        x = x2.flatten(),
        y = y2.flatten(),
        z = H.Rtoplot(x2.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x2label,
        ylabel = y2label,
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[2],
        debug = True,
        x = x3.flatten(),
        y = y3.flatten(),
        z = H.Rtoplot(x3.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x3label,
        ylabel = y3label,
        xlim = [-1, 1],
        ylim = [-2, 2],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = True,
    )    
    f.suptitle(suptitle, fontsize = 12)
    plt.tight_layout()
    f.savefig(fname)
    plt.close(f)

    ############## Neb ##############
    ############## Neb ##############
    ############## Neb ##############

    tau_V_neb_norm_GAL__g = (H.tau_V_neb__rg[10, :] + H.tau_V_neb__rg[9, :] / 2.)
    tau_V_neb_norm__rg = H.tau_V_neb__rg / tau_V_neb_norm_GAL__g
    SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[10, :] + H.aSFRSD_Ha__rg[9, :] / 2.)
    aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g

    x1 = np.ma.copy(np.ma.log10(H.tau_V_neb__g))
    y1 = np.ma.copy(np.ma.log10(H.SFRSD_Ha__g * 1e6)) 
    x2 = np.ma.copy(np.ma.log10(H.tau_V_neb__rg))
    y2 = np.ma.copy(np.ma.log10(H.aSFRSD_Ha__rg * 1e6)) 
    x3 = np.ma.copy(np.ma.log10(tau_V_neb_norm__rg))
    y3 = np.ma.copy(np.ma.log10(aSFRSD_Ha_norm__rg)) 

    x1label = r'$\log\ \tau_V^{neb}$'
    y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
    x2label = r'$\log\ \tau_V^{neb}(R)$'
    y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
    y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'

    if mask_radius is True:
        x1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        y1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        x2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        x3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_neb_mask.png'
    else:
        suptitle = r'NGals:%d  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_neb.png'
    
    f, axArr = plt.subplots(1, 3)
    f.set_dpi(100)
    f.set_size_inches(15, 8)
    C.plot_zbins(
        f = f,
        ax = axArr[0],
        debug = True,
        x = x1,
        y = y1,
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x1label,
        ylabel = y1label,
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[1],
        debug = True,
        x = x2.flatten(),
        y = y2.flatten(),
        xlabel = x2label,
        ylabel = y2label,
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        z = H.Rtoplot(x2.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[2],
        debug = True,
        x = x3.flatten(),
        y = y3.flatten(),
        xlabel = x3label,
        ylabel = y3label,
        xlim = [-1, 1],
        ylim = [-2, 2],
        z = H.Rtoplot(x2.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = True,
    )    
    f.suptitle(suptitle, fontsize = 12)
    plt.tight_layout()
    f.savefig(fname)
    plt.close(f)

    ############## Mixed ##############
    ############## Mixed ##############
    ############## Mixed ############## 
     
    tau_V_neb_norm_GAL__g = (H.tau_V_neb__rg[10, :] + H.tau_V_neb__rg[9, :] / 2.)
    tau_V_neb_norm__rg = H.tau_V_neb__rg / tau_V_neb_norm_GAL__g
    SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][10, :] + H.aSFRSD__Trg[iT][9, :] / 2.)
    aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
 
    x1 = np.ma.copy(np.ma.log10(H.tau_V_neb__g))
    y1 = np.ma.copy(np.ma.log10(H.SFRSD__Tg[iT] * 1e6)) 
    x2 = np.ma.copy(np.ma.log10(H.tau_V_neb__rg))
    y2 = np.ma.copy(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6)) 
    x3 = np.ma.copy(np.ma.log10(tau_V_neb_norm__rg))
    y3 = np.ma.copy(np.ma.log10(aSFRSD_norm__rg)) 
 
    x1label = r'$\log\ \tau_V^{neb}$'
    y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x2label = r'$\log\ \tau_V^{neb}(R)$'
    y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
    y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
 
    if mask_radius is True:
        x1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        y1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        x2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        x3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_mix_mask_%sMyr.png' % str(tSF / 1.e6)
    else:
        suptitle = r'NGals:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_mix_%sMyr.png' % str(tSF / 1.e6)
 
    f, axArr = plt.subplots(1, 3)
    f.set_dpi(100)
    f.set_size_inches(15, 8)
    C.plot_zbins(
        f = f,
        ax = axArr[0],
        debug = True,
        x = x1,
        y = y1,
        xlabel = x1label,
        ylabel = y1label,
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[1],
        debug = True,
        x = x2.flatten(),
        y = y2.flatten(),
        z = H.Rtoplot(x2.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x2label,
        ylabel = y2label,
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[2],
        debug = True,
        x = x3.flatten(),
        y = y3.flatten(),
        z = H.Rtoplot(x3.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x3label,
        ylabel = y3label,
        xlim = [-1, 1],
        ylim = [-2, 2],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = True,
    )    
    f.suptitle(suptitle, fontsize = 12)
    plt.tight_layout()
    f.savefig(fname)
    plt.close(f)

    ############## Mixed ##############
    ############## Mixed ##############
    ############## Mixed ############## 
    tau_V_norm_GAL__g = (H.tau_V__Trg[iT][10, :] + H.tau_V__Trg[iT][9, :] / 2.)
    tau_V_norm__rg = H.tau_V__Trg[iT] / tau_V_norm_GAL__g
    SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[10, :] + H.aSFRSD_Ha__rg[9, :] / 2.)
    aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g
 
    x1 = np.ma.copy(np.ma.log10(H.tau_V__Tg[iT]))
    y1 = np.ma.copy(np.ma.log10(H.SFRSD_Ha__g * 1e6)) 
    x2 = np.ma.copy(np.ma.log10(H.tau_V__Trg[iT]))
    y2 = np.ma.copy(np.ma.log10(H.aSFRSD_Ha__rg * 1e6)) 
    x3 = np.ma.copy(np.ma.log10(tau_V_norm__rg))
    y3 = np.ma.copy(np.ma.log10(aSFRSD_Ha_norm__rg)) 

    x1label = r'$\log\ \tau_V^{\star}$'
    y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
    x2label = r'$\log\ \tau_V^{\star}(R)$'
    y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
    x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
    y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'
 
    if mask_radius is True:
        x1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        y1[~(H.zone_dist_HLR__g > RNuc)] = np.ma.masked
        x2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y2[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        x3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        y3[~(H.RbinCenter__r > RNuc)] = np.ma.masked
        suptitle = r'NGals:%d  R > %.1fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, RNuc, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_mix2_mask_%sMyr.png' % str(tSF / 1.e6)
    else:
        suptitle = r'NGals:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        fname = 'SKeSKradius_mix2_%sMyr.png' % str(tSF / 1.e6)
 
    f, axArr = plt.subplots(1, 3)
    f.set_dpi(100)
    f.set_size_inches(15, 8)
    C.plot_zbins(
        f = f,
        ax = axArr[0],
        debug = True,
        x = x1,
        y = y1,
        xlabel = x1label,
        ylabel = y1label,
        z = H.zone_dist_HLR__g,
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[1],
        debug = True,
        x = x2.flatten(),
        y = y2.flatten(),
        z = H.Rtoplot(x2.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x2label,
        ylabel = y2label,
        xlim = [-1.5, 0.5],
        ylim = [-3.5, 1],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = False,
    )    
    C.plot_zbins(
        f = f,
        ax = axArr[2],
        debug = True,
        x = x3.flatten(),
        y = y3.flatten(),
        z = H.Rtoplot(x3.shape).flatten(),
        zlabel = r'R (HLR)',
        zmask = None,
        zlim = [RNuc, 2],
        xlabel = x3label,
        ylabel = y3label,
        xlim = [-1, 1],
        ylim = [-2, 2],
        ols = True,
        running_stats = True,
        rs_gaussian_smooth = True,
        x_major_locator = 1.,
        x_minor_locator = 0.2,
        y_major_locator = 1.,
        y_minor_locator = 0.2,
        #rs_percentiles = True,
        rs_gs_fwhm = 8,
        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.5),
        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = ''),
        kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True),
        kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
        cb = True,
    )    
    f.suptitle(suptitle, fontsize = 12)
    plt.tight_layout()
    f.savefig(fname)
    plt.close(f)
