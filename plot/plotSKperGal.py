#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import sys
import time
import numpy as np
import argparse as ap
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
from CALIFAUtils.objects import H5SFRData
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.scripts import read_one_cube
from CALIFAUtils.scripts import DrawHLRCircle
from CALIFAUtils.scripts import get_morfologia
from CALIFAUtils.globals import califa_work_dir
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.scripts import calc_running_stats
from CALIFAUtils.scripts import DrawHLRCircleInSDSSImage

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'califaID' : 'K0073',
        'itSF' : 11,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--califaID', '-g',
                        metavar = 'FILE',
                        type = str,
                        default = default['califaID'])
    parser.add_argument('--itSF', '-T',
                        help = 'age index',
                        metavar = '',
                        type = int,
                        default = default['itSF'])

    return parser.parse_args()

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

imgDir = califa_work_dir + 'images/'

Zsun = 0.019

if __name__ == '__main__':
    t_init_gal = time.clock()

    args = parser_args()
    debug = args.debug
    h5fname = args.hdf5
    galName = args.califaID
    iT = args.itSF
    
    H = H5SFRData(args.hdf5)
    tSF__T = H.tSF__T
    ageMyr = tSF__T[iT] / 1e6
    
    if debug:
        print 'califaID: ', galName
        print 'h5fname: ', h5fname
        print 'iTSF: (%.2fMyr)' % (iT, tSF__T[iT] / 1e6)
    
    if (len(np.where(H.califaIDs == galName)[0]) == 0):
        exit('<<< plot: %s: no data.' % galName)
    
    # global
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax
    RbinCenter__r = H.RbinCenter__r
    # ALL gal
    ##stellar
    x_Y__g = H.x_Y__Tg[iT]
    aSFRSD__rg = H.aSFRSD__Trg[iT]
    tau_V__rg = H.tau_V__Trg[iT]
    ##nebular
    aSFRSD_Ha__rg = H.aSFRSD_Ha__rg
    tau_V_neb__rg = H.tau_V_neb__rg
    # one gal
    ##stellar
    tau_V__z = getattr(H, '%s_tau_V__Tg' % galName)[iT]
    atau_V__r = getattr(H, '%s_tau_V__Trg' % galName)[iT]
    SFRSD__z = getattr(H, '%s_SFRSD__Tg' % galName)[iT]
    aSFRSD__r = getattr(H, '%s_aSFRSD__Trg' % galName)[iT]
    x_Y__z = getattr(H, '%s_x_Y__Tg' % galName)[iT]
    ##nebular
    EW_Ha__z = getattr(H, '%s_EW_Ha__g' % galName)
    EW_Hb__z = getattr(H, '%s_EW_Hb__g' % galName)
    tau_V_neb__z = getattr(H, '%s_tau_V_neb__g' % galName)
    atau_V_neb__r = getattr(H, '%s_tau_V_neb__rg' % galName)
    tau_V_neb_err__z = getattr(H, '%s_tau_V_neb_err__g' % galName)
    SFRSD_Ha__z = getattr(H, '%s_SFRSD_Ha__g' % galName)
    aSFRSD_Ha__r = getattr(H, '%s_aSFRSD_Ha__rg' % galName)
    
    galaxyImgFile = imgDir + galName + '.jpg'
    K = read_one_cube(galName, EL = True)
    
    # Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
    
    tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
    
    #stellar
    tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
    SFRSD__yx = K.zoneToYX(SFRSD__z, extensive = False)
    x_Y__yx = K.zoneToYX(x_Y__z, extensive = False)
    #nebular
    tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
    SFRSD_Ha__yx = K.zoneToYX(SFRSD_Ha__z, extensive = False)
    EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
    EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
    tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z.data, extensive = False)
    F_obs_Ha__z = K.EL.flux[K.EL.lines.index('6563'), :]
    F_obs_Ha__yx = K.zoneToYX(F_obs_Ha__z, extensive = True)
    #mixed
    deltaTau__z = tau_V_neb__z - tau_V__z 
    deltaTau__yx = K.zoneToYX(deltaTau__z, extensive = False)
    
    t_calc = time.clock()
    print 'calc: elapsed time: %.2f' % (t_calc - t_init_gal)
    
    NRows = 3
    NCols = 4
    
    f, axArr = plt.subplots(NRows, NCols)
    f.set_size_inches((NCols * 5.34, NRows * 5.))
    
    for ax in f.axes:
        ax.set_axis_off()
    
    age = tSF__T[iT]

    ax = axArr[0, 0]
    ax.set_axis_on()
    galimg = plt.imread(galaxyImgFile)[::-1, :, :]
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    ax.imshow(galimg, origin = 'lower')
    DrawHLRCircleInSDSSImage(ax, K.HLR_pix, pa, ba)
    
    ax = axArr[0, 1]
    ax.set_axis_on()
    xlabel = r'EW(H$\alpha$) [$\AA$]'
    im = ax.imshow(EW_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 20, vmin = 3)
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    sigma_dev = 1
    ax = axArr[0, 2]
    ax.set_axis_on()
    xlabel = r'$F_{obs}(H\alpha)$ [erg cm${}^{-2}$ s${}^{-1}$]'
    vmax = F_obs_Ha__yx.mean() + sigma_dev * F_obs_Ha__yx.std()
    vmin = 0
    im = ax.imshow(F_obs_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[0, 3]
    ax.set_axis_on()
    xlabel = r'$\epsilon\tau_V^{neb}$'
    im = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 1)
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[1, 0]
    ax.set_axis_on()
    x = np.ma.log10(tau_V__rg.flatten())
    y = np.ma.log10(aSFRSD__rg.flatten() * 1e6)
    mask = x.mask | y.mask  
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlabel = r'$\log\ \tau_V^{\star}(R)$'
    ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
    xlim = [np.log10(tauVOkMin), 0.5]
    ylim = [-3.5, 1]
    sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calc_running_stats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
    xbinCenter = aux[0]
    xMedian = aux[1]
    xMean = aux[2]
    xStd = aux[3]
    yMedian = aux[4]
    yMean = aux[5]
    yStd = aux[6]
    nInBin = aux[7]
    xPrc = aux[8]
    yPrc = aux[9]
    ax.plot(xMedian, yMedian, 'k', lw = 2)
    ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
    ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
    txt = '%.2f Myr' % (age / 1e6)
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14)
    ##########################
    x = np.ma.log10(atau_V__r)
    y = np.ma.log10(aSFRSD__r * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
    # vmax = xr.mean() + 2. * xr.std()
    # vmin = xr.mean() - 2. * xr.std()
    # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'jet_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
    cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
    cb.set_label(r'R [HLR]')
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which = 'major')
    
    sigma_dev = 3.

    ax = axArr[1, 1]
    ax.set_axis_on()
    x = np.ma.log10(tau_V__yx)
    y = np.ma.log10(SFRSD__yx * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylim = [-3.5, 1]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # vmax = y.mean() + sigma_dev * y.std()
    # vmin = y.mean() - sigma_dev * y.std()
    # im = ax.imshow(ym, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin, vmax = vmax)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    txt = 'not masked: %d' % (~mask).sum()
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)
    
    ax = axArr[1, 2]
    ax.set_axis_on()
    x = np.ma.log10(tau_V__yx)
    y = np.ma.log10(SFRSD__yx * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlabel = r'$\log\ \tau_V^\star$' 
    ylim = [-3.5, 1]
    vmin = np.log10(tauVOkMin)
    im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
    #im = ax.imshow(x, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[1, 3]
    ax.set_axis_on()
    xlabel = r'$x_Y [\%]$'
    vmin = xOkMin * 100.
    im = ax.imshow(100. * x_Y__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(xlabel, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[2, 0]
    ax.set_axis_on()
    x = np.ma.log10(tau_V_neb__rg.flatten())
    y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlabel = r'$\log\ \tau_V^{neb}(R)$'
    ylabel = r'$\log\ \langle \Sigma_{SFR}^{neb}(R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
    xlim = [np.log10(tauVNebOkMin), 0.5]
    ylim = [-3.5, 1]
    sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calc_running_stats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
    xbinCenter = aux[0]
    xMedian = aux[1]
    xMean = aux[2]
    xStd = aux[3]
    yMedian = aux[4]
    yMean = aux[5]
    yStd = aux[6]
    nInBin = aux[7]
    xPrc = aux[8]
    yPrc = aux[9]
    ax.plot(xMedian, yMedian, 'k', lw = 2)
    ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
    ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
    txt = '%.2f Myr' % (age / 1e6)
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14)
    ##########################
    x = np.ma.log10(atau_V_neb__r)
    y = np.ma.log10(aSFRSD_Ha__r * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
    # vmax = xr.mean() + 2. * xr.std()
    # vmin = xr.mean() - 2. * xr.std()
    # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'jet_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
    cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
    cb.set_label(r'R [HLR]')
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax.grid(which = 'major')

    ax = axArr[2, 1]
    ax.set_axis_on()
    x = np.ma.log10(tau_V_neb__yx)
    y = np.ma.log10(SFRSD_Ha__yx * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    label = r'$\log\ \Sigma_{SFR}^{neb} [M_\odot yr^{-1} kpc^{-2}]$' 
    ylim = [-3.5, 1]
    im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    txt = 'not masked: %d' % (~mask).sum()
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    ax.set_title(label, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)
    
    ax = axArr[2, 2]
    ax.set_axis_on()
    x = np.ma.log10(tau_V_neb__yx)
    y = np.ma.log10(SFRSD_Ha__yx * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    label = r'$\log\ \tau_V^{neb}$' 
    vmin = np.log10(tauVNebOkMin)
    im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(label, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[2, 3]
    ax.set_axis_on()
    label = r'$\delta\ \tau_V$'
    im = ax.imshow(deltaTau__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    mask = deltaTau__yx.mask
    txt = 'not masked: %d' % (~mask).sum()
    plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    ax.set_title(label, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # x = np.ma.log10(tau_V__Tz[iT])
    # y = np.ma.log10(SFRSD__Tz[iT] * 1e6)
    # mask = x.mask | y.mask
    # xm = np.ma.masked_array(x, mask = mask)
    # ym = np.ma.masked_array(y, mask = mask)
    # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
    # xr__yx = K.zoneToYX(xr, extensive = False)
    # xlabel = r'xr'
    # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # # vmax = xr__yx.mean() + 2. * xr__yx.std()
    # # vmin = xr__yx.mean() - 2. * xr__yx.std()
    # # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
    # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    # DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    # txt = 'Nz: %d' % K.N_zone 
    # plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
    # txt = 'masked: %d' % mask.sum() 
    # plot_text_ax(ax, txt, 0.02, 0.92, 14, 'top', 'left')
    # ax.set_title(xlabel, y=-0.15)
    # ax.grid()
    # f.colorbar(ax = ax, mappable = im, use_gridspec = True)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    f.suptitle(r'%s - morph:%s  b/a:%.2f  age:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (galName, tipos, ba, ageMyr, xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize = 24)
    f.subplots_adjust(left = 0.07, bottom = 0.1, right = 0.99, wspace = 0.1, top = 0.9)
    f.savefig('%s_mosaic.png' % galName)
    plt.close(f)
    t_plot = time.clock()
    print 'plot: elapsed time: %.2f' % (t_plot - t_calc)
    print 'total: elapsed time: %.2f' % (t_plot - t_init_gal)
