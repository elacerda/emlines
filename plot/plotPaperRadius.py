#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import CALIFAUtils as C

debug = True
mask_radius = True

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
    
    for iT in iT_values:
        tSF = H.tSF__T[iT]
        xkeys = [ 'atfluxR', 'alogZmassR', 'logO3N2S06R', 'logO3N2M13R', 'logMcorSDR', 'xYR', 'logWHaWHbR' ]
        ykeys = [ 'tauVdiffR', 'tauVRatioR', 'logWHaWHbR' ]
        for fnamepref, xv, yv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys):
            if mask_radius is True:
                m = (H.RbinCenter__r < 0.5) 
                xv['v'][~m] = np.ma.masked
                yv['v'][~m] = np.ma.masked
                filename = '%s_maskradius_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            else:  
                filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            if xk == yk:
                continue
            xv.update(dict(limprc = [0, 100]))
            yv.update(dict(limprc = [0, 100]))
            xv.update(updates[xk])
            yv.update(updates[yk])
            C.plot_zbins(
                debug = debug,
                x = xv['v'],
                y = yv['v'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                xlabel = xv['label'],
                ylabel = yv['label'],
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_suptitle = dict(fontsize = 12),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                filename = filename,
            )
             
            mask = []
            ba = H.reply_arr_by_radius(H.ba_GAL__g)
            balimsup = [0.39, 0.63, 1.]
            titles = [ 
                r'$\frac{b}{a}\ \leq\ 0.39$',
                r'$0.39\ <\ \frac{b}{a}\ \leq\ 0.63$',
                r'$0.63\ <\ \frac{b}{a}\ \leq\ 1.00$',
            ]
            mask.append((ba <= 0.39))
            mask.append((np.bitwise_and(ba > 0.39, ba <= 0.63)))
            mask.append((np.bitwise_and(ba > 0.63, ba <= 1.)))
            if mask_radius is True:
                filename = '%s_bacolors_maskradius_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            else:
                filename = '%s_bacolors_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.plot_zbins(
                debug = debug,
                x = xv['v'],
                xlabel = xv['label'],
                y = yv['v'],
                ylabel = yv['label'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                running_stats = True,
                rs_gaussian_smooth = True,
                rs_percentiles = True,
                rs_gs_fwhm = 4.,
                kwargs_plot_rs = dict(c = 'k', ls = '--', lw = 2, label = 'Median (xy) (run. stats)'),
                kwargs_suptitle = dict(fontsize = 12),
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
                z = H.reply_arr_by_radius(H.ba_GAL__g),
                zlabel = r'$\frac{b}{a}$',
                zname = r'$\frac{b}{a}$',
                zbins = len(titles),
                zbins_labels = titles,
                zbins_mask = mask,
                zmask = False,
                zlimprc = [0, 100],
                zbins_rs_gaussian_smooth = True,
                zbins_rs_gs_fwhm = 4.,
                filename = filename,
            )

            f, axArr = plt.subplots(1, 3)
            f.set_size_inches(14, 6)
            i = 0
            axArr[0].set_xlabel(xv['label'])
            axArr[0].set_ylabel(yv['label'])
            for msk in mask:
                ax = axArr[i]
                kwargs = dict(
                    f = f,
                    ax = ax,
                    debug = debug,
                    x = xv['v'][msk],
                    y = yv['v'][msk],
                    xlim = xv['lim'],
                    ylim = yv['lim'],
                    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                    running_stats = True,
                    rs_gaussian_smooth = True,
                    rs_percentiles = True,
                    rs_gs_fwhm = 4.,
                    kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                    rs_errorbar = False,
                    kwargs_legend = dict(fontsize = 12, loc = 'best'),
                    kwargs_title = dict(fontsize = 10),
                    title = titles[i],
                )
                if i > 0:
                    kwargs.update({'legend' : False})
                C.plot_zbins(**kwargs)
                i += 1
                 
            suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)                            
            f.suptitle(suptitle, fontsize = 12)
            if mask_radius is True:
                filename = '%s_ba_maskradius_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            else:
                filename = '%s_ba_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.debug_var(debug, filename = filename)
            f.savefig(filename)
            plt.close(f)
        ##########################################################################
        ##########################################################################
        ##########################################################################
    
        ##########################################################################
        ##########################################################################
        ##########################################################################
        xkeys = [ 'baR' ]
        ykeys = [ 'tauVdiffR', 'tauVRatioR' , 'logWHaWHbR' ]
        for fnamepref, xv, yv in H.plot_xyz_keys_iter(iT = iT, iU = -1, xkeys = xkeys, ykeys = ykeys):
            if mask_radius is True:
                m = (H.RbinCenter__r < 0.5) 
                xv['v'][~m] = np.ma.masked
                yv['v'][~m] = np.ma.masked
                filename = '%s_maskradius_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            else:  
                filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            xk = fnamepref.split('_')[0]
            yk = fnamepref.split('_')[1]
            xv.update(dict(limprc = [0, 100]))
            yv.update(dict(limprc = [0, 100]))
            xv.update(updates[xk])
            yv.update(updates[yk])
            if mask_radius is True:
                filename = '%s_maskradius_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            else:
                filename = '%s_%.2fMyr.png' % (fnamepref, tSF / 1e6)
            C.plot_zbins(
                debug = debug,
                x = xv['v'],
                y = yv['v'],
                xlim = xv['lim'],
                ylim = yv['lim'],
                xlabel = xv['label'],
                ylabel = yv['label'],
                running_stats = True,
                rs_percentiles = True,
                rs_errorbar = False,
                rs_gaussian_smooth = True,
                rs_gs_fwhm = 4.,
                kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
                kwargs_suptitle = dict(fontsize = 12),
                kwargs_legend = dict(fontsize = 12, loc = 'best'),
                suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                filename = filename,
            )
            zkeys = [ 'atfluxR', 'alogZmassR', 'logO3N2S06R', 'logMcorSDR', 'xYR', 'logWHaWHbR' ]
            for zk in zkeys:
                if zk == yk:
                    continue
                _, zv = H.get_plot_dict(iT = iT, iU = -1, key = zk)
                zv.update(updates[zk])
                zname = H.get_plot_dict(iT = iT, iU = -1)[zk[:-1]]['legendname']
                mask = []
                zlimsup = np.percentile(zv['v'], [ 33, 67, 100 ]).tolist()
                C.debug_var(debug, zlimsup = zlimsup)
                titles = [ 
                    r'%s$\ \leq\ %.2f$' % (zname, zlimsup[0]),
                    r'$%.2f\ <\ $%s$\ \leq\ %.2f$' % (zlimsup[0], zname, zlimsup[1]),
                    r'$%.2f\ <\ $%s$\ \leq\ %.2f$' % (zlimsup[1], zname, zlimsup[2]),
                ]
                C.debug_var(debug, legends = titles)
                mask.append((zv['v'] <= zlimsup[0]))
                mask.append((np.bitwise_and(zv['v'] > zlimsup[0], zv['v'] <= zlimsup[1])))
                mask.append((np.bitwise_and(zv['v'] > zlimsup[1], zv['v'] <= zlimsup[2])))
                if mask_radius is True:
                    filename = '%s_%s_maskradius_%.2fMyr.png' % (fnamepref, zk, tSF / 1e6)
                else:
                    filename = '%s_%s_%.2fMyr.png' % (fnamepref, zk, tSF / 1e6)
                C.plot_zbins(
                    debug = debug,
                    x = xv['v'],
                    xlabel = xv['label'],
                    y = yv['v'],
                    ylabel = yv['label'],
                    xlim = xv['lim'],
                    ylim = yv['lim'],
                    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
                    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
                    running_stats = True,
                    rs_gaussian_smooth = True,
                    rs_percentiles = True,
                    rs_gs_fwhm = 4.,
                    kwargs_plot_rs = dict(c = 'k', ls = '--', lw = 2, label = 'Median (xy) (run. stats)'),
                    kwargs_suptitle = dict(fontsize = 12),
                    suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                    kwargs_legend = dict(fontsize = 12, loc = 'best'),
                    zname = zname,
                    z = zv['v'],
                    zlabel = zv['label'],
                    zbins = len(titles),
                    zbins_labels = titles,
                    zbins_mask = mask,
                    zmask = False,
                    zlim = zv['lim'],
                    zbins_rs_gaussian_smooth = True,
                    zbins_rs_gs_fwhm = 4.,
                    filename = '%s_%s_%.2fMyr.png' % (fnamepref, zk, tSF / 1e6),
                )