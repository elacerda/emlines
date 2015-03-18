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

mpl.rcParams['font.size']       = 17
mpl.rcParams['axes.labelsize']  = 17
mpl.rcParams['axes.titlesize']  = 19
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

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
        print iT, tSF/1e6
        xaxis = {
            'atflux' : dict(
                        v = H.at_flux__g, 
                        label = r'$\langle \log\ t \rangle_L$ [yr]',
                        lim = [7, 10],
                        majloc = 0.6,
                        minloc = 0.12,
                        limprc = [1, 99],
                       ),
            'alogZmass' : dict(
                            v = H.alogZ_mass__Ug[-1], 
                            label = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9),
                            lim = [-2.2, 0.25],
                            majloc = 0.35,
                            minloc = 0.07,
                            limprc = [1, 99],
                          ),
            'logMcorSD' : dict(
                            v = np.ma.log10(H.McorSD__g), 
                            label = r'$\log\ \mu_\star$ [$M_\odot \ pc^{-2}$]',
                            lim = [1, 4.6],
                            majloc = 1.,
                            minloc = 0.2,
                            limprc = [1, 99],                            
                          ),
            'xY' : dict(
                    v = 100. * H.x_Y__Tg[iT], 
                    label = '$x_Y$ [%]',
                    lim = [0, 60],
                    majloc = 12.,
                    minloc = 2.,
                    limprc = [1, 99],                    
                   ),
            'logWHaWHb' : dict(
                            v = np.ma.log10(H.EW_Ha__g/H.EW_Hb__g), 
                            label = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$',
                            lim = [0.2, 0.8],
                            majloc = 0.12,
                            minloc = 0.024,
                            limprc = [1, 99],
                          ),
        }
        yaxis = {
             'tauVdiff' : dict(
                             v = H.tau_V_neb__g - H.tau_V__Tg[iT], 
                             label = r'$\tau_V^{neb}\ -\ \tau_V^\star$',
                             lim = [-1.2, 2.6],
                             majloc = 0.75,
                             minloc = 0.15,
                             limprc = [1, 99],
                          ),
            'tauVRatio' : dict(
                            v = H.tau_V_neb__g/H.tau_V__Tg[iT], 
                            label = r'$\frac{\tau_V^{neb}}{\tau_V^\star}$',
                            lim = [0, 15],
                            majloc = 3.,
                            minloc = 0.6,
                            limprc = [1, 99],
                          ),
             'logWHaWHb' : dict(
                             v = np.ma.log10(H.EW_Ha__g/H.EW_Hb__g), 
                             label = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$',
                             lim = [0.2, 0.8],
                             majloc = 0.12,
                             minloc = 0.024,
                             limprc = [1, 99],
                           ),                 
        }
        for xk, xv in xaxis.iteritems():
            for yk, yv in yaxis.iteritems():
                if xk != yk:
                    plot_zbins(
                        debug = debug,
                        x = xv['v'],
                        xlabel = xv['label'],
                        y = yv['v'],
                        ylabel = yv['label'],
                        xlimprc = xv['limprc'],
                        ylimprc = yv['limprc'],
                        #xlim = xv['lim'],
                        #ylim = yv['lim'],
                        kwargs_figure = dict(figsize=(10,8), dpi = 100),
                        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                        #kwargs_scatter = dict(c = '0.5', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                        running_stats = True,
                        rs_gaussian_smooth = True,
                        rs_percentiles = True,
                        rs_gs_fwhm = 0.4,
                        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                        rs_errorbar = False,
                        kwargs_suptitle = dict(fontsize = 12),
                        suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                        filename = '%s_%s_%.2fMyr.png' % (xk, yk, tSF / 1e6),
                        #x_major_locator = xv['majloc'],
                        #x_minor_locator = xv['minloc'],
                        #y_major_locator = yv['majloc'],
                        #y_minor_locator = yv['minloc'],
                        kwargs_legend = dict(fontsize = 12),
                        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                        # z = H.reply_arr_by_zones(H.ba_GAL__g),
                        # #zbins = [0.39, 0.63, 1.],
                        # zbins = 3,
                        # zmask = False,
                        # zlimprc = [1, 99],
                        # zbins_rs_gaussian_smooth = True,
                        # zbins_rs_gs_fwhm = 0.4,
                        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                    )
                    plot_zbins(
                        debug = debug,
                        x = xv['v'],
                        xlabel = xv['label'],
                        y = yv['v'],
                        ylabel = yv['label'],
                        xlimprc = xv['limprc'],
                        ylimprc = yv['limprc'],
                        #xlim = xv['lim'],
                        #ylim = yv['lim'],
                        kwargs_figure = dict(figsize=(10,8), dpi = 100),
                        kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                        #kwargs_scatter = dict(c = '0.5', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                        running_stats = True,
                        rs_gaussian_smooth = True,
                        rs_percentiles = True,
                        rs_gs_fwhm = 0.4,
                        kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                        rs_errorbar = False,
                        kwargs_suptitle = dict(fontsize = 12),
                        suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                        filename = '%s_%s_bacolors_%.2fMyr.png' % (xk, yk, tSF / 1e6),
                        #x_major_locator = xv['majloc'],
                        #x_minor_locator = xv['minloc'],
                        #y_major_locator = yv['majloc'],
                        #y_minor_locator = yv['minloc'],
                        kwargs_legend = dict(fontsize = 12),
                        z = H.reply_arr_by_zones(H.ba_GAL__g),
                        #zbins = [0.39, 0.63, 1.],
                        zbins = 3,
                        zmask = False,
                        zlimprc = [1, 99],
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
                    f.set_size_inches(14,7)
                    i = 0
                    axArr[0].set_xlabel(xv['label'])
                    axArr[0].set_ylabel(yv['label'])
                    xlim = [ xv['v'].min(), xv['v'].max() ]
                    ylim = [ yv['v'].min(), yv['v'].max() ]
                    for msk in mask:
                        ax = axArr[i]
                        plot_zbins(
                            f = f,
                            ax = ax,
                            debug = debug,
                            x = xv['v'][msk],
                            y = yv['v'][msk],
                            xlim = xlim,
                            ylim = ylim,
                            kwargs_figure = dict(figsize=(10,8), dpi = 100),
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
                    filename = '%s_%s_ba_%.2fMyr.png' % (xk, yk, tSF / 1e6)
                    debug_var(debug, filename = filename)
                    f.savefig(filename)
                    plt.close(f)
