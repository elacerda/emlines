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
from CALIFAUtils.scripts import debug_var
#from CALIFAUtils.plots import add_subplot_axes

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
        
    H = H5SFRData(h5file)
    
    iT_values = [11, 17]
    
    #################################################################################
    #################################################################################
    #################################################################################
    
    for iT in iT_values:
        tSF = H.tSF__T[iT]
        xkeys = [ 'atflux', 'alogZmass', 'OHIICHIM', 'logO3N2S06', 'logO3N2M13', 'logMcorSD', 'xY', 'logWHaWHb' ]
        ykeys = [ 'tauVdiff', 'tauVRatio', 'logWHaWHb' ]
        for fnamepref, xv, yv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys):
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            if xk != yk:
                xv.update(dict(limprc=[0,100]))
                yv.update(dict(limprc=[0,100]))
                if yk == 'tauVRatio':
                    yv.update(dict(lim = [0, 10], limprc = None))
                else:
                    yv.update(dict(lim = None))
                plot_zbins(
                    debug = debug,
                    x = xv['v'],
                    xlabel = xv['label'],
                    y = yv['v'],
                    ylabel = yv['label'],
                    xlimprc = xv.get('limprc', None),
                    ylimprc = yv.get('limprc', None),
                    #xlim = xv['lim'],
                    ylim = yv['lim'],
                    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                    running_stats = True,
                    rs_percentiles = True,
                    rs_errorbar = False,
                    rs_gaussian_smooth = True,
                    rs_gs_fwhm = 0.4,
                    kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                    #x_major_locator = xv['majloc'],
                    #x_minor_locator = xv['minloc'],
                    #y_major_locator = yv['majloc'],
                    #y_minor_locator = yv['minloc'],
                    kwargs_suptitle = dict(fontsize = 12),
                    suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                    filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6),
                    kwargs_legend = dict(fontsize = 12),
                )
                plot_zbins(
                    debug = debug,
                    x = xv['v'],
                    xlabel = xv['label'],
                    y = yv['v'],
                    ylabel = yv['label'],
                    xlimprc = xv.get('limprc', None),
                    ylimprc = yv.get('limprc', None),
                    #xlim = xv['lim'],
                    #ylim = yv['lim'],
                    ylim = yv['lim'],
                    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                    running_stats = True,
                    rs_gaussian_smooth = True,
                    rs_percentiles = True,
                    rs_gs_fwhm = 0.4,
                    kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (xy) (run. stats)'),
                    #rs_errorbar = False,
                    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                    # rs_yx = True, 
                    # rs_yx_gaussian_smooth = True,
                    # rs_yx_percentiles = True,
                    # rs_yx_gs_fwhm = 0.4,
                    # kwargs_plot_rs_yx = dict(c = 'k', lw = 2, label = 'Median (yx) (run. stats)'),
                    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                    #rs_errorbar = False,
                    kwargs_suptitle = dict(fontsize = 12),
                    suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                    filename = '%s_bacolors_%.2fMyr.png' % (fnamepref, tSF / 1e6),
                    #x_major_locator = xv['majloc'],
                    #x_minor_locator = xv['minloc'],
                    #y_major_locator = yv['majloc'],
                    #y_minor_locator = yv['minloc'],
                    kwargs_legend = dict(fontsize = 12),
                    z = H.reply_arr_by_zones(H.ba_GAL__g),
                    #zbins = [0.39, 0.63, 1.],
                    zlabel = r'$\frac{b}{a}$',
                    zbins = 3,
                    zmask = False,
                    zlimprc = [0, 100],
                    zbins_rs_gaussian_smooth = True,
                    zbins_rs_gs_fwhm = 0.4,
                )
    ##########################################################################
    ##########################################################################
    ##########################################################################
                mask = []
                ba = H.reply_arr_by_zones(H.ba_GAL__g)
                balimsup = [0.39, 0.63, 1.]
                mask.append((ba <= 0.39))
                mask.append((np.bitwise_and(ba > 0.39, ba <= 0.63)))
                mask.append((np.bitwise_and(ba > 0.63, ba <= 1.)))
                f, axArr = plt.subplots(1, 3)
                f.set_size_inches(14, 7)
                i = 0
                axArr[0].set_xlabel(xv['label'])
                axArr[0].set_ylabel(yv['label'])
                #xlim = [ xv['v'].min(), xv['v'].max() ]
                #ylim = [ yv['v'].min(), yv['v'].max() ]
                xlim = xv['lim']
                ylim = yv['lim']
                for msk in mask:
                    ax = axArr[i]
                    plot_zbins(
                        f = f,
                        ax = ax,
                        debug = debug,
                        x = xv['v'][msk],
                        y = yv['v'][msk],
                        xlimprc = xv.get('limprc', None),
                        ylimprc = yv.get('limprc', None),
                        #xlim = xlim,
                        #ylim = ylim,
                        kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                        running_stats = True,
                        rs_gaussian_smooth = True,
                        rs_percentiles = True,
                        rs_gs_fwhm = 0.4,
                        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                        rs_errorbar = False,
                        #x_major_locator = xv['majloc'],
                        #x_minor_locator = xv['minloc'],
                        #y_major_locator = yv['majloc'],
                        #y_minor_locator = yv['minloc'],
                        kwargs_legend = dict(fontsize = 12),
                    )
                    i += 1
                    
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)                            
                f.suptitle(suptitle, fontsize = 12)
                filename = '%s_ba_%.2fMyr.png' % (fnamepref, tSF / 1e6)
                debug_var(debug, filename = filename)
                f.savefig(filename)
                plt.close(f)