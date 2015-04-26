#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from CALIFAUtils.plots import plot_zbins
from CALIFAUtils.objects import H5SFRData
from matplotlib.ticker import MultipleLocator

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

iT_values = [4, 10, 11, 17]

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = H5SFRData(h5file)
    tSF__T = np.asarray([H.tSF__T[i] for i in iT_values])
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    
    ###########################################################################
    ###########################################################################
    ###########################################################################
    for iT,tSF in enumerate(tSF__T):        
        x = H.tau_V__Tg[iT]
        y = H.tau_V_neb__g
        z =  H.reply_arr_by_zones(H.morfType_GAL__g)
        mask = x.mask | y.mask
        zm = np.ma.masked_array(z, mask = mask)
        zticks_mask = [(zm > 8.9) & (zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
        zticks = [9., 9.5, 10, 10.5, 11., 11.5]
        zticklabels = ['Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd']

        plot_zbins(
                   zname ='morph',
                   debug = True, 
                   x = H.tau_V__Tg[iT],
                   xlim = [0, 1.5],
                   xlabel = r'$\tau_V^\star$',
                   y = H.tau_V_neb__g,
                   ylim = [0, 2.5],
                   ylabel = r'$\tau_V^{neb}$',
                   z = H.reply_arr_by_zones(H.morfType_GAL__g),
                   zlabel = 'morph. type',
                   zbins_rs_gaussian_smooth = True,
                   zbins_rs_gs_fwhm = 0.4,
                   zbins = len(zticks_mask),
                   zbins_mask = zticks_mask,
                   zticks = zticks,  
                   zticklabels = zticklabels,
                   zlim = [9, 11.5],
                   kwargs_scatter = dict(marker = 'o', s = 6, edgecolor = 'none', alpha = 0.4, label = ''),
                   suptitle = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1.e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax),
                   filename = 'tauV_tauVneb_morphType_%.2fMyr.png' % (tSF / 1e6),
                   ols = True,
                   running_stats = True,
                   kwargs_plot_rs = dict(c = 'k', lw = 2),
                   rs_errorbar = False,
                   kwargs_ols_plot = dict(ls = '--', lw = 0.7, label = ''),
                   x_major_locator = 0.25,
                   x_minor_locator = 0.05,
                   y_major_locator = 0.25,
                   y_minor_locator = 0.05,
                   kwargs_legend = dict(fontsize = 8),
        )
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     ###########################################################################
#     ###########################################################################
#     ###########################################################################
#     #tauVtoAV = np.log10(np.exp(1)) / 0.4
#     zticks = [9., 9.5, 10, 10.5, 11., 11.5]
# 
#     zticklabels = ['Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd']
#     zticklabelcolors = [(1, 0, 0), (1, .5, 0), (0, 1, 0), (0, 1, .5), (0, 0.5, 1), (0, 0, 1)]
#     
#     zticklabels2 = ['Sa + Sab', 'Sb', 'Sbc', 'Sc + Scd']
#     zticklabelcolors2 = [(1, 0, 0), (0.5, 1, 0), (0, .5, 1), (0, 0, 1)]
#     
#     zticks3 = [ -0.8, -0.6, -0.4, -0.2, -0.1, 0, 0.1, 0.2 ]
#
#     for iT,tSF in enumerate(tSF__T):        
#         x = H.tau_V__Trg[iT]
#         y = H.tau_V_neb__rg
#         z = H.reply_arr_by_radius(H.morfType_GAL__g)
#         mask = x.mask | y.mask
#         xm = np.ma.masked_array(x, mask = mask)
#         ym = np.ma.masked_array(y, mask = mask)
#         zm = np.ma.masked_array(z, mask = mask)
#         xlabel = r'$\tau_V^\star(R)$'
#         ylabel = r'$\tau_V^{neb}(R)$'
#         zlabel = r'morf. type'
#         fname = 'tauVR_tauVNebR_age_%sMyr.png' % str(tSF / 1.e6)
#         xlim = [0, 1.5]
#         ylim = [0, 2.5]
#         f = plt.figure()
#         f.set_size_inches(10, 8)
#         ax = f.gca()
#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(ylabel)
#         ax.set_xlim(xlim)
#         ax.set_ylim(ylim)
#         sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', marker = 'o', s = 15., edgecolor = 'none', vmax = 11.5, vmin = 9., alpha = 0.6)
#         a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
#         zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
#         cb = f.colorbar(sc, ticks=zticks)
#         cb.set_label(zlabel)
#         cb.ax.set_yticklabels(zticklabels)
#         #####################
#         # y holding x0=0
#         #####################
#         A = ym.sum() / xm.sum()
#         Y = A * xm
#         Yrms = (ym - Y).std()
#         ax.plot(xm, Y, c = '0.5', ls = '--', lw = 0.5)
#         txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
#         plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
#         #####################
#         ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
#         ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
#         for i, msk in enumerate(zticks_mask):
#             print zticklabels2[i]
#             y_pos = 0.07 + 0.05 * i
#             X = xm[msk]
#             Y = ym[msk]
#             a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
#                                                          0.98, y_pos, 12, 
#                                                          zticklabelcolors2[i], 
#                                                          rms = True, 
#                                                          label = zticklabels2[i],
#                                                          text = False)
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         # for i, morfType in enumerate(zticks):
#         #     X = xm[zm == morfType]
#         #     Y = ym[zm == morfType]
#         #     if len(X) > 3 and len(Y) > 3:
#         #         nBox = 50
#         #         dxBox       = (X.max() - X.min()) / (nBox - 1.)
#         #         aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
#         #         xbinCenter  = aux[0]
#         #         xMedian     = aux[1]
#         #         xMean       = aux[2]
#         #         xStd        = aux[3]
#         #         yMedian     = aux[4]
#         #         yMean       = aux[5]
#         #         yStd        = aux[6]
#         #         nInBin      = aux[7]
#         #         xPrc        = aux[8]
#         #         yPrc        = aux[9]
#         #         FWHM = 0.4
#         #         xM = np.ma.masked_array(xMedian)
#         #         yM = np.ma.masked_array(yMedian)
#         #         mask = np.isnan(xM) | np.isnan(yM) 
#         #         Xs, Ys = gaussSmooth_YofX(xM[~mask], yM[~mask], FWHM)
#         #         ax.plot(Xs, Ys, c = zticklabelcolors[i], lw = 2, label = '%s' % zticklabels[i])
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         ax.legend(loc = 'upper left', fontsize = 14)
#         ax.xaxis.set_major_locator(MultipleLocator(0.25))
#         ax.xaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.yaxis.set_major_locator(MultipleLocator(0.25))
#         ax.yaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.grid(which = 'major')
#         txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
#         f.suptitle(txt, fontsize = 14)
#         #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#         f.savefig(fname)
#         plt.close(f)
# 
#         x = H.tau_V__Tg[iT]
#         y = H.tau_V_neb__g
#         z = H.reply_arr_by_zones(H.morfType_GAL__g)
#         mask = x.mask | y.mask
#         xm = np.ma.masked_array(x, mask = mask)
#         ym = np.ma.masked_array(y, mask = mask)
#         zm = np.ma.masked_array(z, mask = mask)
#         zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
#         xlabel = r'$\tau_V^\star$'
#         ylabel = r'$\tau_V^{neb}$'
#         zlabel = r'morf. type'
#         fname = 'tauV_tauVNeb_age_%sMyr.png' % str(tSF / 1.e6)
#         xlim = [0, 1.5]
#         ylim = [0, 2.5]
#         f = plt.figure()
#         f.set_size_inches(10, 8)
#         ax = f.gca()
#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(ylabel)
#         ax.set_xlim(xlim)
#         ax.set_ylim(ylim)
#         sc = ax.scatter(xm, ym, c = zm, cmap = 'spectral_r', marker = 'o', s = 15., edgecolor = 'none', vmax = 11.5, vmin = 9., alpha = 0.3)
#         a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 16, 'k', rms = True)
#         cb = f.colorbar(sc, ticks=zticks)
#         cb.set_label(zlabel)
#         cb.ax.set_yticklabels(zticklabels)
#         #####################
#         # y holding x0=0
#         #####################
#         A = ym.sum() / xm.sum()
#         Y = A * xm
#         Yrms = (ym - Y).std()
#         ax.plot(ax.get_xlim(), A * np.asarray(ax.get_xlim()), c = '0.5', ls = '-', lw = 1.5)
#         txt = r'y$(x_{0}=0)$ = %.2fx $y_{rms}$:%.2f' % (A, Yrms)
#         plot_text_ax(ax, txt, 0.98, 0.08, 16, 'bottom', 'right', color = '0.5')
#         #####################
#         ax.plot(ax.get_xlim(), (1/0.44) * np.asarray(ax.get_xlim()), ls = '--', lw = 2., c = '0.3', label = r'$\tau_V^\star = 0.44\tau_V^{neb}$')
#         ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".7", lw = 2., label = r'$\tau_V^\star = \tau_V^{neb}$')
#         for i, msk in enumerate(zticks_mask):
#             print zticklabels2[i]
#             y_pos = 0.13 + 0.05 * i
#             X = xm[msk]
#             Y = ym[msk]
#             a, b, sigma_a, sigma_b = plotOLSbisectorAxis(ax, X, Y, 
#                                                          0.98, y_pos, 12, 
#                                                          zticklabelcolors2[i], 
#                                                          rms = True, 
#                                                          label = zticklabels2[i],
#                                                          text = False)
#         ax.legend(loc = 'upper left', fontsize = 14)
#         ax.xaxis.set_major_locator(MultipleLocator(0.25))
#         ax.xaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.yaxis.set_major_locator(MultipleLocator(0.25))
#         ax.yaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.grid(which = 'major')
#         txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
#         f.suptitle(txt, fontsize = 14)
#         #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#         f.savefig(fname)
#         plt.close(f)
# 
#         x = H.tau_V__Tg[iT]
#         y = H.tau_V_neb__g
#         mask = x.mask | y.mask
#         xm = np.ma.masked_array(x, mask = mask)
#         ym = np.ma.masked_array(y, mask = mask)
#         a, b, sigma_a, sigma_b = OLS_bisector(xm,ym)
#         R = ym - (a * xm + b)
#         z = H.reply_arr_by_zones(H.morfType_GAL__g)
#         #Z = H.alogZ_mass__Ug[-1]
#         mask = R.mask
#         #mask = R.mask | Z.mask
#         xm = np.ma.masked_array(x, mask = mask)
#         ym = np.ma.masked_array(y, mask = mask)
#         zm = np.ma.masked_array(z, mask = mask)
#         #Zm = np.ma.masked_array(Z, mask = mask)
#         zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
#         xlabel = r'$\tau_V^\star$'
#         ylabel = r'$\tau_V^{neb}\ -\ y_{OLS}$'
#         zlabel = r'morf. type'
#         #Zlabel = r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]'
#         fname = 'tauV_ResTauVOLS_morf_%sMyr.png' % str(tSF / 1.e6)
#         #fname = 'tauV_ResTauVOLS_alogZmass_morf_%sMyr.png' % str(tSF / 1.e6) 
#         f = plt.figure()
#         f.set_size_inches(10, 8)
#         ax = f.gca()
#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(ylabel)
#         xlim = [0, 1.5]
#         ylim = [-3, 3]
#         ax.set_xlim(xlim)
#         ax.set_ylim(ylim)
#         if b > 0:
#             txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ + %.2f' %  (a, b)
#         else:
#             txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ - %.2f' %  (a, b * -1.)
#         plot_text_ax(ax, txt, 0.98, 0.02, 14, 'bottom', 'right', color = 'k')
#         sc = ax.scatter(xm, R, c = zm, cmap = 'spectral_r', marker = 'o', s = 15., edgecolor = 'none', vmax = 11.5, vmin = 9., alpha = 0.3)
#         #sc = ax.scatter(xm, R, c = Zm, cmap = 'spectral_r', marker = 'o', s = 15., edgecolor = 'none', alpha = 0.3)    
#         cb = f.colorbar(sc, ticks=zticks)
#         cb.set_label(zlabel)
#         cb.ax.set_yticklabels(zticklabels)
#         #cb = f.colorbar(sc)
#         #cb.set_label(Zlabel)
#         for i, msk in enumerate(zticks_mask):
#             X = xm[msk]
#             Y = R[msk]
#             if len(X) > 3 and len(Y) > 3:
#                 nBox = 50
#                 dxBox       = (X.max() - X.min()) / (nBox - 1.)
#                 aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
#                 xbinCenter  = aux[0]
#                 xMedian     = aux[1]
#                 xMean       = aux[2]
#                 xStd        = aux[3]
#                 yMedian     = aux[4]
#                 yMean       = aux[5]
#                 yStd        = aux[6]
#                 nInBin      = aux[7]
#                 xPrc        = aux[8]
#                 yPrc        = aux[9]
#                 FWHM = 0.4
#                 xM = np.ma.masked_array(xMedian)
#                 yM = np.ma.masked_array(yMedian)
#                 mask = np.isnan(xM) | np.isnan(yM) 
#                 Xs, Ys = gaussSmooth_YofX(xM[~mask], yM[~mask], FWHM)
#                 ax.plot(Xs, Ys, c = zticklabelcolors2[i], lw = 2, label = '%s' % zticklabels2[i])
#         ax.legend(loc = 'upper right', fontsize = 14)
#         ax.xaxis.set_major_locator(MultipleLocator(0.25))
#         ax.xaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.yaxis.set_major_locator(MultipleLocator(1.00))
#         ax.yaxis.set_minor_locator(MultipleLocator(0.20))
#         ax.grid(which = 'major')
#         txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
#         f.suptitle(txt, fontsize = 14)
#         #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#         f.savefig(fname)
#         plt.close(f)
# 
#         x = H.tau_V__Trg[iT]
#         y = H.tau_V_neb__rg
#         mask = x.mask | y.mask
#         xm = np.ma.masked_array(x, mask = mask)
#         ym = np.ma.masked_array(y, mask = mask)
#         a, b, sigma_a, sigma_b = OLS_bisector(xm,ym)
#         R = ym - (a * xm + b)
#         z = H.reply_arr_by_radius(H.morfType_GAL__g)
#         zm = np.ma.masked_array(z, mask = R.mask)
#         zticks_mask = [(zm <= 9.5), (zm == 10), (zm == 10.5), (zm >= 11.)]
#         xlabel = r'$\tau_V^\star$'
#         ylabel = r'$\tau_V^{neb}\ -\ y_{OLS}$'
#         zlabel = r'morf. type'
#         fname = 'tauVR_ResTauVOLSR_morf_%sMyr.png' % str(tSF / 1.e6) 
#         f = plt.figure()
#         f.set_size_inches(10, 8)
#         ax = f.gca()
#         ax.set_xlabel(xlabel)
#         ax.set_ylabel(ylabel)
#         xlim = [0, 1.5]
#         ylim = [-3, 3]
#         ax.set_xlim(xlim)
#         ax.set_ylim(ylim)
#         if b > 0:
#             txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ + %.2f' %  (a, b)
#         else:
#             txt = r'$y_{OLS}$ = %.2f$\tau_V^\star$ - %.2f' %  (a, b * -1.)
#         plot_text_ax(ax, txt, 0.98, 0.02, 14, 'bottom', 'right', color = 'k')
#         sc = ax.scatter(xm, R, c = zm, cmap = 'spectral_r', marker = 'o', s = 15., edgecolor = 'none', vmax = 11.5, vmin = 9., alpha = 0.3)    
#         cb = f.colorbar(sc, ticks=zticks)
#         cb.set_label(zlabel)
#         cb.ax.set_yticklabels(zticklabels)
#         for i, msk in enumerate(zticks_mask):
#             X = xm[msk]
#             Y = R[msk]
#             if len(X) > 3 and len(Y) > 3:
#                 nBox = 50
#                 dxBox       = (X.max() - X.min()) / (nBox - 1.)
#                 aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
#                 xbinCenter  = aux[0]
#                 xMedian     = aux[1]
#                 xMean       = aux[2]
#                 xStd        = aux[3]
#                 yMedian     = aux[4]
#                 yMean       = aux[5]
#                 yStd        = aux[6]
#                 nInBin      = aux[7]
#                 xPrc        = aux[8]
#                 yPrc        = aux[9]
#                 FWHM = 0.4
#                 xM = np.ma.masked_array(xMedian)
#                 yM = np.ma.masked_array(yMedian)
#                 mask = np.isnan(xM) | np.isnan(yM) 
#                 Xs, Ys = gaussSmooth_YofX(xM[~mask], yM[~mask], FWHM)
#                 ax.plot(Xs, Ys, c = zticklabelcolors2[i], lw = 2, label = '%s' % zticklabels2[i])
#         ax.legend(loc = 'upper right', fontsize = 14)
#         ax.xaxis.set_major_locator(MultipleLocator(0.25))
#         ax.xaxis.set_minor_locator(MultipleLocator(0.05))
#         ax.yaxis.set_major_locator(MultipleLocator(1.00))
#         ax.yaxis.set_minor_locator(MultipleLocator(0.20))
#         ax.grid(which = 'major')
#         txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF__T[iT] / 1.e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax)
#         f.suptitle(txt, fontsize = 14)
#         #plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#         f.savefig(fname)
#         plt.close(f)
# 
#         
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
