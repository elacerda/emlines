#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import sys
from CALIFAUtils.plots import plotScatterColorAxis, plot_text_ax, plot_zbins
from CALIFAUtils.objects import H5SFRData, read_kwargs

debug = True

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
img_output_ext = 'png'
    

def f_plot(**kwargs):
    args = read_kwargs(**kwargs)
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
        
    H = H5SFRData(h5file)
    iU = -1

    #################################################################################
    # SKzero = np.log10(1.6e-4)
    # SKslope = 1.4
    # logSigmaGas = (np.log10(SFRSD_Ha__g * 1e6) - SKzero) / SKslope
    # c = np.log10(0.2)
    # logDGR = c + np.log10(tau_V_neb__g) - logSigmaGas
    # logO_H = logZ_neb_S06__g + np.log10(4.9e-4)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    #DGR(SigmaGas - KS)
    SKzero = 1.6e-4
    SKslope = 1.4
    dustdim = 0.2
    DGR__g = dustdim * H.tau_V__Tg[iT] / (H.SFRSD__Tg[iT]/SKzero) ** (1./SKslope)
    DGR_Ha__g = dustdim * H.tau_V_neb__g / (H.SFRSD_Ha__g/SKzero) ** (1./SKslope)
    DGR__rg = dustdim * H.tau_V_neb__rg / (H.aSFRSD_Ha__rg/1.6e-4) ** (1./SKslope)
    ########################
    
    DGR = 10. ** (-2.21)
    k = dustdim / DGR
    k__g = dustdim / DGR__g
    k__rg = dustdim / DGR__rg
    #SigmaGas__g = k * H.tau_V_neb__g
    SigmaGas__g = k * H.tau_V__Tg[iT]
    #SigmaGas__g = k__g * H.tau_V_neb__g
    #SigmaGas__g = k__g * H.tau_V__Tg[iT]
    #SigmaGas__rg = k * H.tau_V_neb__rg
    SigmaGas__rg = k * H.tau_V__Trg[iT]
    #SigmaGas__rg = k__rg * H.tau_V_neb__rg
    #SigmaGas__rg = k__rg * H.tau_V__Trg[iT]
    f_gas__g = 1. / (1. + (H.McorSD__Tg[iT] / SigmaGas__g))
    f_gas__rg = 1. / (1. + (H.McorSD__Trg[iT] / SigmaGas__rg))
    RbinCenter__rg = ((np.ones_like(f_gas__rg).T * H.RbinCenter__r).T).flatten()

    #################################################################################
    #################################################################################
    #################################################################################
    xaxis = {
        'alogZmass' : dict(
                        v = H.alogZ_mass__Ug[-1] + np.log10(4.9e-4) + 12, 
                        label = r'12 + $\log\ O/H$ ($\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr))' % (H.tZ__U[-1] / 1e9),
                        #lim = [7., 9.5],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                      ),
        'OHIICHIM' : dict(
                        v = H.O_HIICHIM__g, 
                        label = r'12 + $\log\ O/H$ (HII-CHI-mistry, EPM, 2014)',
                        #lim = [7., 9.5],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                     ),
        'logO3N2S06' : dict(
                        v = H.logZ_neb_S06__g + np.log10(4.9e-4) + 12, 
                        label = r'12 + $\log\ O/H$ (logO3N2, Stasinska, 2006)',
                        #lim = [7., 9.5],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                       ),
        'logO3N2M13' : dict(
                        v = H.O_O3N2_M13__g, 
                        label = r'12 + $\log\ O/H$ (logO3N2, Marino, 2013)',
                        #lim = [7., 9.5],
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                       ),
        'logfgas' : dict(
                        v = np.ma.log10(f_gas__g),
                        label = r'$\log\ f_{gas}$',
                        #lim =,
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                    ),
        'logtauV' : dict(
                        v = np.ma.log10(H.tau_V__Tg[iT]), 
                        label = r'$\log\ \tau_V^\star$',
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                    ),
        'logtauVneb' : dict(
                        v = np.ma.log10(H.tau_V_neb__g), 
                        label = r'$\log\ \tau_V^{neb}$',
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                       ),
    }
    yaxis = {
        'logSFRSD' : dict(
                        v = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), 
                        label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$',
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                     ),
        'logSFRSDHa' : dict(
                        v = np.ma.log10(H.SFRSD_Ha__g * 1e6), 
                        label = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$',
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                       ),
        'logDGR' : dict(
                    v = np.ma.log10(DGR__g),
                    label = r'$\log$ DGR',
                    majloc = 0.5,
                    minloc = 0.1,
                    limprc = [0, 100],
                   ),
        'logDGRHa' : dict(
                        v = np.ma.log10(DGR_Ha__g),
                        label = r'$\log$ DGR',
                        majloc = 0.5,
                        minloc = 0.1,
                        limprc = [0, 100],
                   ),
    }
    for xk, xv in xaxis.iteritems():
        for yk, yv in yaxis.iteritems():
            if xk != yk:
                tSF = H.tSF__T[iT]
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
                )
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plotArgs = [
    #     dict(x = H.logZ_neb_S06__g,
    #          y = np.ma.log10(DGR__g),
    #          z = np.ma.log10(H.McorSD__g), 
    #          xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]',
    #          ylabel = r'$\log$ DGR',
    #          zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]',
    #          fname_pref = 'logZneb_logDGR_McorSD',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None,
    #          x_major_locator = 0.1,
    #          x_minor_locator = 0.02,
    #          y_major_locator = 0.5,
    #          y_minor_locator = 0.1,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     dict(x = H.alogZ_mass__Ug[iU],
    #          y = np.ma.log10(DGR__g),
    #          z = H.dist_zone__g,  
    #          xlabel = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[iU] / 1e9),
    #          ylabel = r'$\log$ DGR',
    #          zlabel = r'zone distance [HLR]',
    #          fname_pref = 'alogZmass_logDGR_zoneDistance',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None,
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.5,
    #          y_minor_locator = 0.1,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     dict(x = H.logZ_neb_S06__g, 
    #          y = np.ma.log10(DGR__g),
    #          z = np.ma.log10(H.reply_arr_by_zones(H.McorSD_GAL__g)),
    #          xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]',
    #          ylabel = r'$\log$ DGR',
    #          zlabel = r'$\log\ M_\star$ [$M_\odot\ pc^{-2}$]',
    #          fname_pref = 'logZneb_logDGR_McorSDGAL', 
    #          xlim = None,
    #          ylim = None,
    #          zlim = None,
    #          x_major_locator = 0.1,
    #          x_minor_locator = 0.02,
    #          y_major_locator = 0.5,
    #          y_minor_locator = 0.1,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #     # dict(
    #     #      x = np.ma.log10(f_gas__rg.flatten()), 
    #     #      y = H.alogZ_mass__Urg[-1].flatten(), 
    #     #      z = RbinCenter__rg, 
    #     #      xlabel = r'$\log\ f_{gas}(R)$', 
    #     #      ylabel = r'$\langle \log\ Z_\star \rangle_M(R)$ [$Z_\odot$]', 
    #     #      zlabel = r'R [HLR]', 
    #     #      fname_pref = 'logfgas_alogZmass_radius',
    #     #      #xlim = [-7.5, -3.5],
    #     #      xlim = None,
    #     #      ylim = [-1.5, 0.3], 
    #     #      zlim = None, 
    #     #      x_major_locator = 0.5,
    #     #      x_minor_locator = 0.1,
    #     #      y_major_locator = 0.2,
    #     #      y_minor_locator = 0.04,
    #     #      contour = False, 
    #     #      run_stats = True, 
    #     #      OLS = False
    #     #      ),
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #     dict(x = np.ma.log10(f_gas__rg.flatten()),
    #          y = H.logZ_neb_S06__rg.flatten(),
    #          z = RbinCenter__rg,
    #          xlabel = r'$\log\ f_{gas}(R)$',
    #          ylabel = r'$\langle \log\ Z_{neb} \rangle(R)$ [$Z_\odot$]',
    #          zlabel = r'R [HLR]',
    #          fname_pref = 'logfgas_logZneb_radius',
    #          #xlim = [-7.5, -3.5],
    #          xlim = None,
    #          ylim = None,
    #          #ylim = [-0.5, 0.2],
    #          zlim = None,
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.1,
    #          y_minor_locator = 0.02,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     dict(x = np.ma.log10(f_gas__g),
    #          y = H.logZ_neb_S06__g,
    #          z = H.dist_zone__g,
    #          xlabel = r'$\log\ f_{gas}$',
    #          ylabel = r'$\langle \log\ Z_{neb} \rangle$ [$Z_\odot$]',
    #          zlabel = r'zone distance [HLR]',
    #          fname_pref = 'logfgas_logZneb_zoneDistance',
    #          #xlim = [-7.5, -3.5],
    #          xlim = None,
    #          ylim = None,
    #          #ylim = [-0.5, 0.2],
    #          zlim = None,
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.1,
    #          y_minor_locator = 0.02,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #     # dict(x = np.ma.log10(f_gas__rg.flatten()),
    #     #      y = np.ma.log10(H.McorSD__Trg[iT].flatten()),
    #     #      z = RbinCenter__rg,
    #     #      xlabel = r'$\log\ f_{gas}(R)$',
    #     #      ylabel = r'$\log\ \langle \mu_\star \rangle(R)$ [$M_\odot \ pc^{-2}$]',
    #     #      zlabel = r'R [HLR]',
    #     #      fname_pref = 'logfgas_logMcorSD_radius',
    #     #      #xlim = [-7.5, -3.5],
    #     #      xlim = None,
    #     #      ylim = [1, 4.7],
    #     #      zlim = None,
    #     #      x_major_locator = 0.5,
    #     #      x_minor_locator = 0.1,
    #     #      y_major_locator = 0.5,
    #     #      y_minor_locator = 0.1,
    #     #      contour = False, 
    #     #      run_stats = True, 
    #     #      OLS = False,
    #     #      ),
    #     # dict(x = np.ma.log10(f_gas__rg.flatten()),
    #     #      y = np.ma.log10(SigmaGas__rg.flatten()),
    #     #      z = RbinCenter__rg,
    #     #      xlabel = r'$\log\ f_{gas}(R)$',
    #     #      ylabel = r'$\log\ \langle \Sigma_{gas} \rangle(R)$ [$M_\odot \ pc^{-2}$]',
    #     #      zlabel = r'R [HLR]',
    #     #      fname_pref = 'logfgas_logSigmaGas_radius',
    #     #      #xlim = [-7.5, -3.5],
    #     #      xlim = None,
    #     #      ylim = None,
    #     #      zlim = None,
    #     #      x_major_locator = 0.5,
    #     #      x_minor_locator = 0.1,
    #     #      y_major_locator = 0.25,
    #     #      y_minor_locator = 0.05,
    #     #      contour = False, 
    #     #      run_stats = True, 
    #     #      OLS = False,
    #     #      ),
    #     # dict(x = np.ma.log10(f_gas__rg.flatten()),
    #     #      y = np.ma.log10(H.McorSD__Trg[iT].flatten() / H.aSFRSD__Trg[iT].flatten()),
    #     #      z = RbinCenter__rg,
    #     #      xlabel = r'$\log\ f_{gas}(R)$',
    #     #      ylabel = r'$\log\ \langle \frac{\mu_\star}{\Sigma_{SFR}} \rangle$ [yr]',
    #     #      zlabel = r'R [HLR]',
    #     #      fname_pref = 'logfgas_McorSD_SFRSD_radius',
    #     #      #xlim = [-7.5, -3.5],
    #     #      xlim = None,
    #     #      ylim = None,
    #     #      zlim = None,
    #     #      x_major_locator = 0.5,
    #     #      x_minor_locator = 0.1,
    #     #      y_major_locator = 0.5,
    #     #      y_minor_locator = 0.1,
    #     #      contour = False, 
    #     #      run_stats = True, 
    #     #      OLS = False,
    #     #      ),
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #     dict(x = np.ma.log10(f_gas__g),
    #          y = np.ma.masked_where(H.O_HIICHIM__g == 0.,H.O_HIICHIM__g, copy=True),
    #          z = H.dist_zone__g,
    #          xlabel = r'$\log\ f_{gas}$', 
    #          ylabel = r'O_HIICHIM', 
    #          zlabel = r'zone distance [HLR]', 
    #          fname_pref = 'logfgas_O_HIICHIM_zoneDistance',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None, 
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.2,
    #          y_minor_locator = 0.04,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False
    #          ),
    #     dict(x = np.ma.masked_where(H.O_HIICHIM__g == 0.,H.O_HIICHIM__g, copy=True),
    #          y = np.ma.log10(DGR__g),
    #          z = H.dist_zone__g,  
    #          xlabel = 'O_HIICHIM', 
    #          ylabel = r'$\log$ DGR',
    #          zlabel = r'zone distance [HLR]',
    #          fname_pref = 'OHIICHIM_logDGR_zoneDistance',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None,
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.5,
    #          y_minor_locator = 0.1,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     dict(x = np.ma.masked_where(H.O_O3N2_M13__g == 0.,H.O_O3N2_M13__g, copy=True),
    #          y = np.ma.log10(DGR__g),
    #          z = H.dist_zone__g,  
    #          xlabel = 'O3N2_M13', 
    #          ylabel = r'$\log$ DGR',
    #          zlabel = r'zone distance [HLR]',
    #          fname_pref = 'O3N2M13_logDGR_zoneDistance',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None,
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.5,
    #          y_minor_locator = 0.1,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False,
    #          ),
    #     dict(x = np.ma.log10(f_gas__g),
    #          y = np.ma.masked_where(H.O_O3N2_M13__g == 0.,H.O_O3N2_M13__g, copy=True),
    #          z = H.dist_zone__g,
    #          xlabel = r'$\log\ f_{gas}$', 
    #          ylabel = r'O3N2_M13', 
    #          zlabel = r'zone distance [HLR]', 
    #          fname_pref = 'logfgas_O3N2M13_zoneDistance',
    #          xlim = None,
    #          ylim = None,
    #          zlim = None, 
    #          x_major_locator = 0.5,
    #          x_minor_locator = 0.1,
    #          y_major_locator = 0.2,
    #          y_minor_locator = 0.04,
    #          contour = False, 
    #          run_stats = True, 
    #          OLS = False
    #          ),                
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #     # dict(x = np.ma.log10(SigmaGas__rg.flatten()), 
    #     #      y = np.ma.log10(H.aSFRSD_Ha__rg.flatten() * 1e6),
    #     #      #y = np.ma.log10(H.aSFRSD__Trg[iT].flatten() * 1e6),
    #     #      z = RbinCenter__rg.flatten(),
    #     #      xlabel = r'$\log\ \Sigma_{gas}(R)\ [M_\odot yr^{-1} pc^{-2}]$',
    #     #      ylabel = r'$\log\ \overline{\Sigma_{SFR}}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$',
    #     #      zlabel = r'R [HLR]',
    #     #      fname = 'logSigmaGas_logSigmaSFR_radius_%dgals',
    #     #      xlim = [-4.5, -1.5],
    #     #      ylim = None,
    #     #      zlim = None,
    #     #      x_major_locator = 0.5,
    #     #      x_minor_locator = 0.1,
    #     #      y_major_locator = 0.5,
    #     #      y_minor_locator = 0.1,
    #     #      contour = False, 
    #     #      run_stats = True, 
    #     #      OLS = True,
    #     #      ),
    #     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # ]
    #     
    # for D in plotArgs:
    #     tSF = H.tSF__T[iT]
    #     fname = '%s_%.2fMyr_%dgals.%s' % (D['fname_pref'], tSF/1e6, H.N_gals, img_output_ext)
    #     D.update(tSF = H.tSF__T[iT], H = H, DGR = DGR, fname = fname)
    #     f_plot(**D)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            