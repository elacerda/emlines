#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import CALIFAUtils
import matplotlib as mpl
from matplotlib import pyplot as plt
from CALIFAUtils.plots import plot_zbins, plot_text_ax

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
        
    H = CALIFAUtils.H5SFRData(h5file)
    
    iT_values = [11, 17]
    
    #################################################################################
    #################################################################################
    #################################################################################
    for iT in iT_values:
        tSF = H.tSF__T[iT]

        f = plt.figure()
        f.set_size_inches(10, 8)
        f.set_dpi(100)
        ax = f.gca()
        #ax.set_ylabel(r'$\tau_V (R)$')
        ax.set_ylabel(r'$\tau_V$')
         
        #xkeys = [ 'alogSFRSDR' ]
        #ykeys = [ 'tauVR', 'tauVNebR' ]

        xkeys = [ 'logSFRSD' ]
        ykeys = ['tauV', 'tauVNeb']
        ykeys = ['logtauV', 'logtauVNeb']
         
        color_scatter = { ykeys[-1] : 'b', ykeys[0] : 'r' }
        color_rs = { ykeys[-1] : '#0088ff', ykeys[0] : '#ff3333' }
        pos_ols_xy = { ykeys[-1] : [ 0.99, 0.02 ], ykeys[0] : [ 0.99, 0.07 ] }
         
        for fnamepref, xv, yv in H.plot_xyz_keys_iter(xkeys = xkeys, ykeys = ykeys, iT = 11, iU = -1):
            print xv
            ax.set_xlabel(xv['label'])
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            c_scatter = color_scatter[yk]
            c_rs = color_rs[yk]
            pos_x = pos_ols_xy[yk][0]
            pos_y = pos_ols_xy[yk][1]
            plot_zbins(
                f = f,
                ax = ax,
                debug = True,
                x = xv['v'],
                y = yv['v'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                ols = True,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(c = c_scatter, marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
                kwargs_ols = dict(c = c_rs, pos_x = pos_x, pos_y = pos_y, fs = 12, rms = True, text = True),
                kwargs_ols_plot = dict(c = c_rs, ls = '--', lw = 2, label = yv['label']),
                #running_stats = True,
                #rs_gaussian_smooth = True,
                #rs_percentiles = True,
                #rs_gs_fwhm = 0.4,
                #kwargs_plot_rs = dict(c = c_rs, lw = 2, label = 'Median (run. stats)'),
                #rs_errorbar = False,
                #x_major_locator = xv['majloc'],
                #x_minor_locator = xv['minloc'],
                #y_major_locator = yv['majloc'],
                #y_minor_locator = yv['minloc'],
            )
     
        suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)                            
        f.suptitle(suptitle, fontsize = 12)
        filename = '%s_%.2fMyr.png' % ('logSFRSD_tauV', tSF / 1e6)
        #filename = '%s_%.2fMyr.png' % ('logSFRSDR_tauVR', tSF / 1e6)
        CALIFAUtils.debug_var(True, filename = filename)
        f.savefig(filename)
        plt.close(f)
     
        #################################################################################
         
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xkeys = [ 'logtauVR', 'logtauVNebR' ]
        # ykeys = [ 'alogSFRSDR', 'alogSFRSDHaR' ]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xkeys = [ 'logtauV', 'logtauVNeb' ]
        ykeys = [ 'logSFRSD', 'logSFRSDHa' ]
        for fnamepref, xv, yv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys):
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            plot_zbins(
                debug = True,
                x = xv['v'],
                y = yv['v'],
                xlabel = xv['label'],
                ylabel = yv['label'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(c = 'b', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
                ols = True,
                kwargs_ols = dict(c = 'k', pos_x = 0.99, pos_y = 0.02, fs = 12, rms = True, text = True),
                kwargs_ols_plot = dict(c = 'k', ls = '-.', lw = 1.2, label = 'OLS'),
                running_stats = True,
                rs_gaussian_smooth = True,
                rs_percentiles = True,
                rs_gs_fwhm = 0.4,
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                kwargs_suptitle = dict(fontsize = 12),
                kwargs_legend = dict(loc = 'best'),
                filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6),
            )

        #################################################################################

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xkeys = [ 'logtauVR', 'logtauVNebR' ]
        # ykeys = [ 'alogSFRSDR', 'alogSFRSDHaR' ]
        # zkeys = [ 'logO3N2M13R', 'logO3N2S06R', 'alogZmassR' ]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xkeys = [ 'logtauV', 'logtauVNeb' ]
        ykeys = [ 'logSFRSD', 'logSFRSDHa' ]
        zkeys = [ 'logO3N2M13', 'logO3N2S06', 'alogZmass' ]
        for fnamepref, xv, yv, zv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys, zkeys = zkeys):
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            zk = fnamepref.split('_')[2]
            plot_zbins(
                debug = True,
                x = xv['v'],
                y = yv['v'],
                z = zv['v'],
                zbins = zv.get('bins', 4),
                zbins_mask = zv.get('ticks_mask', None),
                ticklabels = zv.get('ticklabels', None),
                zticks = zv.get('ticks', None),
                zmask = zv.get('mask', None),
                zbins_rs_gaussian_smooth = True,
                zbins_rs_gs_fwhm = 0.4,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                #kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (xy) (run. stats)'),
                kwargs_legend = dict(loc = 'best'),
                kwargs_suptitle = dict(fontsize = 14),
                xlabel = xv['label'],
                ylabel = yv['label'],
                zlabel = zv['label'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                #zlim = zv['lim'],
                zlimprc = zv.get('limprc', [2, 98]),
                x_major_locator = xv['majloc'],
                x_minor_locator = xv['minloc'],
                y_major_locator = yv['majloc'],
                y_minor_locator = yv['minloc'],
                suptitle = suptitle,
                filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6),
            )

        #################################################################################

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xkeys = [ 'logtauVR', 'logtauVNebR' ]
        # ykeys = [ 'alogSFRSDR', 'alogSFRSDHaR' ]
        # zkeys = [ 'morfTypeR' ]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xkeys = [ 'logtauV', 'logtauVNeb' ]
        ykeys = [ 'logSFRSD', 'logSFRSDHa' ]
        zkeys = [ 'morfType' ]
        for fnamepref, xv, yv, zv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys, zkeys = zkeys):
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            zk = fnamepref.split('_')[2]
            mask = xv['v'].mask | yv['v'].mask
            zm = np.ma.masked_array(zv['v'], mask = mask)
            ticks_mask = [(zm > 8.9) & (zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.) & (zm > 9) & (zm <= 11.5)]
            zv.update(
                dict(
                   ticks_mask = ticks_mask,
                   ticks = [9., 9.5, 10, 10.5, 11., 11.5],
                   ticklabels = ['Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd'],
                   bins = len(ticks_mask),
                   bins_mask = ticks_mask,
               )
            )
            plot_zbins(
                debug = True,
                x = xv['v'],
                y = yv['v'],
                z = zv['v'],
                zbins = zv.get('bins', 4),
                zbins_mask = zv.get('ticks_mask', None),
                zticklabels = zv.get('ticklabels', None),
                zticks = zv.get('ticks', None),
                zmask = zv.get('mask', None),
                zbins_rs_gaussian_smooth = True,
                zbins_rs_gs_fwhm = 0.4,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_legend = dict(loc = 'best'),
                kwargs_suptitle = dict(fontsize = 14),
                xlabel = xv['label'],
                ylabel = yv['label'],
                zlabel = zv['label'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                zlim = zv['lim'],
                #zlimprc = [0, 100], 
                x_major_locator = xv['majloc'],
                x_minor_locator = xv['minloc'],
                y_major_locator = yv['majloc'],
                y_minor_locator = yv['minloc'],
                suptitle = suptitle,
                filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6),
            )

        #################################################################################

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xkeys = [ 'logtauVR', 'logtauVNebR' ]
        # ykeys = [ 'alogSFRSDR', 'alogSFRSDHaR' ]
        # zkeys = [ 'baR' ]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xkeys = [ 'logtauV', 'logtauVNeb' ]
        ykeys = [ 'logSFRSD', 'logSFRSDHa' ]
        zkeys = [ 'ba' ]
        for fnamepref, xv, yv, zv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys, zkeys = zkeys):
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            zk = fnamepref.split('_')[2]
            zv.update(dict(limprc = [0, 100], bins = 3))
            plot_zbins(
                debug = True,
                x = xv['v'],
                y = yv['v'],
                z = zv['v'],
                zbins = zv['bins'],
                zmask = zv.get('mask', None),
                zbins_rs_gaussian_smooth = True,
                zbins_rs_gs_fwhm = 0.4,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_legend = dict(loc = 'best'),
                kwargs_suptitle = dict(fontsize = 14),
                xlabel = xv['label'],
                ylabel = yv['label'],
                zlabel = zv['label'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                zlim = zv['lim'],
                zlimprc = zv['limprc'],
                x_major_locator = xv['majloc'],
                x_minor_locator = xv['minloc'],
                y_major_locator = yv['majloc'],
                y_minor_locator = yv['minloc'],
                suptitle = suptitle,
                filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6),
            )

        #################################################################################

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVR')
        # yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVNebR')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauV')
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVNeb')
        f = plt.figure()
        f.set_size_inches(10, 8)
        f.set_dpi(100)
        ax = f.gca()
        #ax.set_ylabel(r'$\tau_V$')
        ax.set_ylabel(r'$\tau_V (R)$')
        plot_zbins(
            f = f,
            ax = ax,
            debug = True,
            x = xv['v'],
            y = yv['v'],
            xlabel = xv['label'],
            ylabel = yv['label'],
            xlim = xv['lim'],
            ylim = yv['lim'],
            kwargs_figure = dict(figsize = (10, 8), dpi = 100),
            kwargs_scatter = dict(c = 'b', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
            ols = True,
            kwargs_ols = dict(c = 'k', pos_x = 0.99, pos_y = 0.02, fs = 12, rms = True, text = True),
            kwargs_ols_plot = dict(c = 'k', ls = '-.', lw = 1.2, label = 'OLS'),
            running_stats = True,
            rs_gaussian_smooth = True,
            rs_percentiles = True,
            rs_gs_fwhm = 0.4,
            kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
            kwargs_legend = dict(loc = 'best'),
        )     
        x = xv['v']
        y = yv['v']
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        #####################
        # y holding x0=0
        #####################
        A = ym.sum() / xm.sum()
        Y = A * xm
        Yrms = (ym - Y).std()
        ax.plot(ax.get_xlim(), A * np.asarray(ax.get_xlim()), c = '0.5', ls = '-.', lw = 1.5, label = r'OLS${}_{x0=0}$')
        txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
        plot_text_ax(ax, txt, 0.99, 0.07, 12, 'bottom', 'right', color = '0.5')
        #####################
        suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)                            
        f.suptitle(suptitle, fontsize = 12)
        filename = '%s_%.2fMyr.png' % ('tauVR_tauVNebR_yholdingx00', tSF / 1e6)
        #filename = '%s_%.2fMyr.png' % ('tauV_tauVNeb_yholdingx00', tSF / 1e6)
        CALIFAUtils.debug_var(True, filename = filename)
        f.savefig(filename)
        plt.close(f)
      
        #################################################################################
 
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVR')
        # yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVNebR')
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xk, xv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauV')
        yk, yv = H.get_plot_dict(iT = iT, iU = -1, key = 'tauVNeb')

        Rmin = 0.2
        mask = np.bitwise_or(xv['v'] >= 1, H.zone_dist_HLR__g < Rmin)
        #mask = (xv['v'][1:, :] >= 1)
        plot_zbins(
            debug = True,
            x = np.ma.masked_array(xv['v'], mask = (xv['v']).mask | mask),
            #x = np.ma.masked_array(xv['v'][1:, :], mask = (xv['v'][1:,:]).mask | mask),
            y = yv['v'],
            #y = yv['v'][1:, :],
            xlabel = xv['label'],
            ylabel = yv['label'],
            xlim = [0, 1.],
            #ylim = yv['lim'],
            ylim = [0, 2.5],
            kwargs_figure = dict(figsize = (10, 8), dpi = 100),
            kwargs_scatter = dict(c = 'b', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.3, label = ''),
            ols = True,
            kwargs_ols = dict(c = 'k', pos_x = 0.99, pos_y = 0.02, fs = 12, rms = True, text = True),
            kwargs_ols_plot = dict(c = 'k', ls = '-.', lw = 2, label = 'OLS'),
            running_stats = True,
            rs_gaussian_smooth = True,
            rs_percentiles = True,
            rs_gs_fwhm = 0.4,
            kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
            suptitle = r'NGals:%d  R >= %.1f HLR  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, Rmin, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
            kwargs_suptitle = dict(fontsize = 12),
            kwargs_legend = dict(loc = 'best'),
            filename = '%s_%s_mask_%.2fMyrs.png' % (xk, yk, tSF / 1e6),
        )

        #################################################################################

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xkeys = [ 'tauVR', 'tauVNebR' ]
        # ykeys = [ 'alogSFRSDR', 'alogSFRSDHaR' ]
        # zkeys = [ 'baR' ]
        # for fnamepref, xv, yv, zv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys, zkeys = zkeys):
        #     xk = fnamepref.split('_')[0]
        #     yk = fnamepref.split('_')[1]
        #     zk = fnamepref.split('_')[2]
        #     zv.update(dict(limprc = [0, 100], bins = 3))
        #     z = np.asarray(zv['v'], dtype = xv['v'].dtype)
        #     x = np.ma.log10(xv['v'] * z)
        #     plot_zbins(
        #         debug = True,
        #         x = x,
        #         y = yv['v'],
        #         z = z,
        #         zbins = zv['bins'],
        #         zmask = zv.get('mask', None),
        #         zbins_rs_gaussian_smooth = True,
        #         zbins_rs_gs_fwhm = 0.4,
        #         kwargs_figure = dict(figsize = (10, 8), dpi = 100),
        #         kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
        #         kwargs_legend = dict(loc = 'best'),
        #         kwargs_suptitle = dict(fontsize = 14),
        #         xlabel = r'$\log\ (\tau_V ^\star\ \frac{b}{a})$',
        #         ylabel = yv['label'],
        #         zlabel = zv['label'],
        #         #xlim = xv['lim'],
        #         #ylim = yv['lim'],
        #         zlim = zv['lim'],
        #         zlimprc = zv['limprc'],
        #         #x_major_locator = xv['majloc'],
        #         #x_minor_locator = xv['minloc'],
        #         #y_major_locator = yv['majloc'],
        #         #y_minor_locator = yv['minloc'],
        #         suptitle = suptitle,
        #         filename = '%s_%.2fMyrs.png' % (fnamepref, tSF / 1e6),
        #     )
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


sys.exit('bai')

dict_aux = dict(
    running_stats = True,
    rs_percentiles = True,
    rs_errorbar = False,
    rs_gaussian_smooth = True,
    rs_gs_fwhm = 0.4,
    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
    kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
    kwargs_suptitle = dict(fontsize = 12),
    kwargs_legend = dict(fontsize = 12),
)
