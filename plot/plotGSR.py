#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
from matplotlib.pyplot import cm
import sys
import CALIFAUtils as C
from CALIFAUtils.plots import plotScatterColorAxis, plot_text_ax, plot_zbins
import CALIFAUtils

debug = True
mask_radius = True
#RNuc = 0.5
#RNuc = 1.0

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
img_output_ext = 'png'


class read_kwargs(object):
    def __init__(self, default_val = None, **kwargs):
        self.default_val = None
        self.kwargs = kwargs  
        
    def __getattr__(self, attr):
        return self.kwargs.get(attr, self.default_val)


def f_plot(**kwargs):
    args = read_kwargs(**kwargs)
    C.debug_var(args.debug, kwargs = kwargs)
    H = args.H
    mask = args.x.mask | args.y.mask
    xm = np.ma.masked_array(args.x, mask = mask)
    ym = np.ma.masked_array(args.y, mask = mask)
    zm = np.ma.masked_array(args.z, mask = mask)
    f = plt.figure()
    f.set_size_inches(10, 8)
    ax = f.gca()
    plotScatterColorAxis(f, xm, ym, zm, args.xlabel, args.ylabel, args.zlabel, args.xlim, args.ylim, args.zlim, contour = args.contour, run_stats = args.run_stats, OLS = args.OLS)
    ax.xaxis.set_major_locator(MultipleLocator(args.x_major_locator))
    ax.xaxis.set_minor_locator(MultipleLocator(args.x_minor_locator))
    ax.yaxis.set_major_locator(MultipleLocator(args.y_major_locator))
    ax.yaxis.set_minor_locator(MultipleLocator(args.y_minor_locator))
    txt = r'DGR = $10^{%.2f}$' % (np.log10(args.DGR))
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    ax.grid(which = 'major')
    f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (args.tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax), fontsize=14)
    f.savefig(args.fname)
    plt.close(f)


if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE [age_index]' % (sys.argv[0])
        exit(1)
        
    try:
        iT = int(sys.argv[2])
    except IndexError:
        iT = 11
        print 'using default value for age index (%d)' % iT
        
    H = C.H5SFRData(h5file)
    iU = -1
    
    ########################
    #SK
    SK_zero = 1.6e-4
    SK_slope = 1.4
    pckpcconv = 1e6 ** (1./SK_slope)
    SK_SigmaGas__g = pckpcconv * (H.SFRSD__Tg[iT]/SK_zero) ** (1./SK_slope)
    SK_SigmaGas_Ha__g = pckpcconv * (H.SFRSD_Ha__g/SK_zero) ** (1./SK_slope)
    SK_SigmaGas__rg = pckpcconv * (H.aSFRSD__Trg[iT]/SK_zero) ** (1./SK_slope)
    SK_SigmaGas_Ha__rg = pckpcconv * (H.aSFRSD_Ha__rg/SK_zero) ** (1./SK_slope)
    SK_integrated_SigmaGas = pckpcconv * (H.integrated_SFRSD__Tg[iT]/SK_zero) ** (1./SK_slope)
    SK_integrated_SigmaGas_Ha = pckpcconv * (H.integrated_SFRSD__Tg[iT]/SK_zero) ** (1./SK_slope) 
    #SK_DGR
    dustdim = 0.2 # md / rhod
    SK_DGR__g = dustdim * H.tau_V__Tg[iT] / SK_SigmaGas__g
    SK_DGR_Ha__g = dustdim * H.tau_V_neb__g / SK_SigmaGas_Ha__g
    SK_DGR__rg = dustdim * H.tau_V__Trg[iT] / SK_SigmaGas__rg
    SK_DGR_Ha__rg = dustdim * H.tau_V_neb__rg / SK_SigmaGas_Ha__rg
    SK_GSR__g = SK_SigmaGas__g / H.McorSD__Tg[iT]
    SK_GSR_Ha__g = SK_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    SK_GSR__rg = SK_SigmaGas__rg / H.McorSD__Trg[iT]
    SK_GSR_Ha__rg = SK_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    #SK gas fraction
    SK_f_gas__g = 1. / (1. + (H.McorSD__Tg[iT] / SK_SigmaGas__g))
    SK_f_gas_Ha__g = 1. / (1. + (H.McorSD__Tg[iT] / SK_SigmaGas_Ha__g))
    SK_f_gas__rg = 1. / (1. + (H.McorSD__Trg[iT] / SK_SigmaGas__rg))
    SK_f_gas_Ha__rg = 1. / (1. + (H.McorSD__Trg[iT] / SK_SigmaGas_Ha__rg))
    ########################
    
    ########################
    # Remy-Ruyer
    # logDGR = c + np.log10(tau_V_neb__g) - logSigmaGas
    DGR = 10. ** (-2.21)
    #DGR__g = 10 ** (-1. * (2.21 + (8.69 - H.O_O3N2_M13__g)))
    #DGR = 10 ** (-1. * (2.21 + 2.02*(8.69 - H.O_O3N2_M13__rg)))
    k = dustdim / DGR
    SigmaGas__g = k * H.tau_V__Tg[iT]
    SigmaGas_Ha__g = k * H.tau_V_neb__g
    SigmaGas__rg = k * H.tau_V__Trg[iT]
    SigmaGas_Ha__rg = k * H.tau_V_neb__rg
    f_gas__g = 1. / (1. + (H.McorSD__Tg[iT] / SigmaGas__g))
    f_gas_Ha__g = 1. / (1. + (H.McorSD__Tg[iT] / SigmaGas_Ha__g))
    f_gas__rg = 1. / (1. + (H.McorSD__Trg[iT] / SigmaGas__rg))
    f_gas_Ha__rg = 1. / (1. + (H.McorSD__Trg[iT] / SigmaGas_Ha__rg))
    ########################
    # logO_H = logZ_neb_S06__g + np.log10(4.9e-4)

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plt.clf()
    # plt.scatter(H.Rtoplot().flatten(), np.ma.log10(SK_GSR__rg.flatten()), c = 'r', alpha = 0.1, s = 10, edgecolor = 'none')
    # plt.plot(H.RbinCenter__r, np.median(np.ma.log10(SK_GSR__rg), axis = 1), '*r', label = r'$\Sigma_{gas}$ from SK with $\Sigma_{SFR}^\star$', ms = 15)
    # plt.scatter(H.Rtoplot().flatten(), np.ma.log10(SK_GSR_Ha__rg.flatten()), c = 'b', alpha = 0.1, s = 10, edgecolor = 'none')
    # plt.plot(H.RbinCenter__r, np.ma.median(np.ma.log10(SK_GSR_Ha__rg), axis = 1), '*b', label = r'$\Sigma_{gas}$ from SK with $\Sigma_{SFR}^{H\alpha}$', ms = 15)
    # plt.scatter(H.Rtoplot().flatten(), np.ma.log10((f_gas__rg.flatten()/(1. - f_gas__rg.flatten()))), c = 'g', alpha = 0.1, s = 10, edgecolor = 'none')
    # plt.plot(H.RbinCenter__r, np.ma.median(np.ma.log10(f_gas__rg/(1. - f_gas__rg)), axis = 1), '*g', label = r'$\Sigma_{gas}$ from $\tau_V$', ms = 15)
    # plt.scatter(H.Rtoplot().flatten(), np.ma.log10((f_gas_Ha__rg.flatten()/(1. - f_gas_Ha__rg.flatten()))), c = 'c', alpha = 0.1, s = 10, edgecolor = 'none')
    # plt.plot(H.RbinCenter__r, np.ma.median(np.ma.log10(f_gas_Ha__rg/(1. - f_gas_Ha__rg)), axis = 1), '*c', label = r'$\Sigma_{gas}$ from $\tau_V^{neb}$', ms = 15)
    # plt.ylim(-4, 1)
    # plt.xlim(0, 2)
    # plt.ylabel(r'$\log$ gas-to-stars')
    # plt.xlabel(r'radius [HLR]')
    # plt.legend(fontsize=10)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plt.clf()
    # plt.scatter(H.Rtoplot().flatten(), np.ma.log10(SK_DGR__rg.flatten()), c = 'r', alpha = 0.1, s = 10, edgecolor = 'none')
    # plt.plot(H.RbinCenter__r, np.median(np.ma.log10(SK_DGR__rg), axis = 1), '*r', ms = 15)
    # plt.ylim(0, 2.5)
    # plt.xlim(0, 2)
    # plt.ylabel(r'$\log$ dust-to-gas')
    # plt.xlabel(r'radius [HLR]')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
        
    xaxis = {
        'R' : dict(
                v = H.Rtoplot(), 
                label = r'radius [HLR]',
                lim = [ 0, 2.0 ],
                majloc = 0.5,
                minloc = 0.1,
                limprc = [0, 100],
              ),
    }
    
    yaxis = {
        'logGasToStars' : dict(
                        v = np.ma.log10(SK_GSR__rg), 
                        label = r'$\log\ $gas-to-stars',
                        lim = [ -4, 1 ],
                        majloc = 1.,
                        minloc = 0.25,
                        limprc = [0, 100],
        ),
        'logDustToGas' : dict(
                        v = np.ma.log10(SK_DGR__rg), 
                        label = r'$\log$ dust-to-gas [from $\tau_V^\star$]',
                        lim = [ -3.5, -1.5 ],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
        ),
        'logDustHaToGas' : dict(
                        v = np.ma.log10(SK_DGR_Ha__rg), 
                        label = r'$\log$ dust-to-gas [from $\tau_V^{neb}$]',
                        lim = [ -3.5, -1.5 ],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
        ),
    }

    zk = 'McorGAL'
    zv = dict(
            v = np.ma.log10((H.Mcor_GAL__g[..., np.newaxis] * np.ones((20))).T),
            label = r'$\log\ M_\star^{GAL}$ [$M_\odot$]', 
    )
    
    for xk, xv in xaxis.iteritems():
        for yk, yv in yaxis.iteritems():
            if xk != yk:
                tSF = H.tSF__T[iT]
                f = plt.figure()
                f.set_dpi(100)
                f.set_size_inches(14, 8)
                ax = f.gca()
                i = 0
                filename = '%s_%s_%s_%.2fMyr.png' % (xk, yk, zk, tSF / 1e6)
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
                xm, ym, zm = CALIFAUtils.ma_mask_xyz(x = xv['v'], y = yv['v'], z = np.ma.log10((H.Mcor_GAL__g[..., np.newaxis] * np.ones((20))).T))
                zticks_mask = [
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 8.8), np.ma.less(zm, 9.6)),
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 9.6), np.ma.less(zm, 10.1)),
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 10.1), np.ma.less(zm, 10.6)),
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 10.6), np.ma.less(zm, 10.9)),
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 10.9), np.ma.less(zm, 11.2)),
                    np.ma.bitwise_and(np.ma.greater_equal(zm, 11.2), np.ma.less(zm, 11.6)),
                ]                    
                zbins_labels = ['9.6--8.8', '10.1--9.6', '10.6--10.1', '10.9--10.6', '11.2--10.9', '11.6--11.2']
                kw = plot_zbins(
                    return_kwargs = True,
                    f = f,
                    ax = ax,
                    debug = debug,
                    x = xv['v'],
                    y = yv['v'],
                    z = zv['v'],
                    zmask = False, 
                    xlabel = xv['label'],
                    ylabel = yv['label'],
                    zlabel = zv['label'],
                    xlim = xv['lim'],
                    ylim = yv['lim'],
                    zlim = [8.8, 11.6],
                    #zlim = zv['lim'],
                    #zmask = True,
                    #zticks = zticks,  
                    #zticklabels = zticklabels,
                    zbins = len(zticks_mask),
                    zbins_labels = zbins_labels,
                    zbins_rs_gaussian_smooth = True,
                    zbins_rs_gs_fwhm = 4,
                    zbins_rs_frac_box = 20,
                    zbins_mask = zticks_mask,
                    #zbins_colors = zbins_colors,
                    x_major_locator = xv['majloc'],
                    x_minor_locator = xv['minloc'],
                    y_major_locator = yv['majloc'],
                    y_minor_locator = yv['minloc'],
                    legend = True,
                    cmap = cm.RdYlBu,
                    #cmap = cm.hot,
                    cb = True,
                    #running_stats = True,
                    #rs_gaussian_smooth = True,
                    #rs_percentiles = True,
                    #rs_gs_fwhm = 8,
                    #rs_frac_box = 20,
                    #kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                    #rs_errorbar = False,
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8, label = ''),
                    kwargs_legend = dict(loc = 'lower right', fontsize = 12, frameon = False),
                    #kwargs_suptitle = dict(fontsize = 12),
                )
                ax = kw['ax']
                init_pos_x, final_pos_y, pos_y_step = 0.02, 0.02, 0.04
                fs = 14
                for j in xrange(len(zticks_mask)):
                    c = kw['zbins_colors'][j]
                    pos_y = final_pos_y + pos_y_step * (len(zticks_mask) - j - 1)
                    txt = str(len(zm[zticks_mask[j]].compressed())) 
                    kw_text = dict(pos_x = init_pos_x, pos_y = pos_y, fs = fs, va = 'bottom', ha = 'left', c = c)
                    plot_text_ax(ax, txt, **kw_text)
                f.suptitle(suptitle, fontsize = 12)
                f.savefig(filename)
                plt.close(f)


    zk = 'McorGAL'
    zv = dict(
            v = np.ma.log10((H.Mcor_GAL__g[..., np.newaxis] * np.ones((20))).T),
            label = r'$\log\ M_\star^{GAL}$ [$M_\odot$]',
            lim = [8.8, 11.6], 
    )
    
    Nbins = 4
    for xk, xv in xaxis.iteritems():
        for yk, yv in yaxis.iteritems():
            if xk != yk:
                tSF = H.tSF__T[iT]
                f = plt.figure()
                f.set_dpi(100)
                f.set_size_inches(14, 8)
                ax = f.gca()
                i = 0
                filename = '%s_%s_%s_%.2fMyr_4bins.png' % (xk, yk, zk, tSF / 1e6)
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
                kw = plot_zbins(
                    return_kwargs = True,
                    f = f,
                    ax = ax,
                    debug = debug,
                    x = xv['v'],
                    y = yv['v'],
                    z = zv['v'],
                    zmask = False, 
                    xlabel = xv['label'],
                    ylabel = yv['label'],
                    zlabel = zv['label'],
                    xlim = xv['lim'],
                    ylim = yv['lim'],
                    zlim = zv['lim'],
                    #zlim = zv['lim'],
                    #zmask = True,
                    #zticks = zticks,  
                    #zticklabels = zticklabels,
                    zbins = Nbins,
                    #zbins_labels = zbins_labels,
                    zbins_rs_gaussian_smooth = True,
                    zbins_rs_gs_fwhm = 4,
                    zbins_rs_frac_box = 20,
                    #zbins_mask = zticks_mask,
                    #zbins_colors = zbins_colors,
                    x_major_locator = xv['majloc'],
                    x_minor_locator = xv['minloc'],
                    y_major_locator = yv['majloc'],
                    y_minor_locator = yv['minloc'],
                    legend = True,
                    cmap = cm.RdYlBu,
                    #cmap = cm.hot,
                    cb = True,
                    #running_stats = True,
                    #rs_gaussian_smooth = True,
                    #rs_percentiles = True,
                    #rs_gs_fwhm = 8,
                    #rs_frac_box = 20,
                    #kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                    #rs_errorbar = False,
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8, label = ''),
                    kwargs_legend = dict(loc = 'lower right', fontsize = 12, frameon = False),
                    #kwargs_suptitle = dict(fontsize = 12),
                )
                ax = kw['ax']
                init_pos_x, final_pos_y, pos_y_step = 0.02, 0.02, 0.04
                fs = 14
                zticks_mask = kw['zbins_mask']
                zm = kw['zm']
                for j in xrange(len(zticks_mask)):
                    c = kw['zbins_colors'][j]
                    pos_y = final_pos_y + pos_y_step * (len(zticks_mask) - j - 1)
                    txt = str(len(zm[zticks_mask[j]].compressed())) 
                    kw_text = dict(pos_x = init_pos_x, pos_y = pos_y, fs = fs, va = 'bottom', ha = 'left', c = c)
                    plot_text_ax(ax, txt, **kw_text)
                f.suptitle(suptitle, fontsize = 12)
                f.savefig(filename)
                plt.close(f)
                
    zk, zv = H.get_plot_dict(iT, -1, key = 'logMcorSDR')
    Nbins = 6
    
    for xk, xv in xaxis.iteritems():
        for yk, yv in yaxis.iteritems():
            if xk != yk:
                tSF = H.tSF__T[iT]
                f = plt.figure()
                f.set_dpi(100)
                f.set_size_inches(14, 8)
                ax = f.gca()
                i = 0
                filename = '%s_%s_%s_%.2fMyr.png' % (xk, yk, zk, tSF / 1e6)
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
                kw = plot_zbins(
                    return_kwargs = True,
                    f = f,
                    ax = ax,
                    debug = debug,
                    x = xv['v'],
                    y = yv['v'],
                    z = zv['v'],
                    zmask = False, 
                    xlabel = xv['label'],
                    ylabel = yv['label'],
                    zlabel = zv['label'],
                    xlim = xv['lim'],
                    ylim = yv['lim'],
                    zlim = zv['lim'],
                    #zlim = zv['lim'],
                    #zmask = True,
                    #zticks = zticks,  
                    #zticklabels = zticklabels,
                    zbins = Nbins,
                    #zbins_labels = zbins_labels,
                    zbins_rs_gaussian_smooth = True,
                    zbins_rs_gs_fwhm = 4,
                    zbins_rs_frac_box = 20,
                    #zbins_mask = zticks_mask,
                    #zbins_colors = zbins_colors,
                    x_major_locator = xv['majloc'],
                    x_minor_locator = xv['minloc'],
                    y_major_locator = yv['majloc'],
                    y_minor_locator = yv['minloc'],
                    legend = True,
                    cmap = cm.RdYlBu,
                    #cmap = cm.hot,
                    cb = True,
                    #running_stats = True,
                    #rs_gaussian_smooth = True,
                    #rs_percentiles = True,
                    #rs_gs_fwhm = 8,
                    #rs_frac_box = 20,
                    #kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                    #rs_errorbar = False,
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8, label = ''),
                    kwargs_legend = dict(loc = 'lower right', fontsize = 12, frameon = False),
                    #kwargs_suptitle = dict(fontsize = 12),
                )
                ax = kw['ax']
                init_pos_x, final_pos_y, pos_y_step = 0.02, 0.02, 0.04
                fs = 14
                zticks_mask = kw['zbins_mask']
                zm = kw['zm']
                print zm.compressed().shape
                for j in xrange(len(zticks_mask)):
                    c = kw['zbins_colors'][j]
                    pos_y = final_pos_y + pos_y_step * (len(zticks_mask) - j - 1)
                    txt = str(len(zm[zticks_mask[j]].compressed())) 
                    kw_text = dict(pos_x = init_pos_x, pos_y = pos_y, fs = fs, va = 'bottom', ha = 'left', c = c)
                    plot_text_ax(ax, txt, **kw_text)
                f.suptitle(suptitle, fontsize = 12)
                f.savefig(filename)
                plt.close(f)                            