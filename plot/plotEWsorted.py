#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import numpy as np
from pycasso import fitsQ3DataCube
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
from get_morfologia import get_morfologia
from scipy import stats as st
from lines import *
import os
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
import time
import sys
import argparse as ap
from plot_aux import H5SFRData, plotOLSbisectorAxis, \
                     plot_text_ax, calcRunningStats, \
                     DrawHLRCircleInSDSSImage, DrawHLRCircle

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

default = {
    'debug' : False,
    'hdf5' : None,
    'califaID' : 'K0073',
    'itSF' : 11,
}

def parser_args():
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
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

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
#galaxiesListFile    = CALIFAWorkDir + 'listAll.txt'
baseCode = 'Bgsd6e'
#versionSuffix       = 'px1_q043.d14a'
versionSuffix = 'v20_q043.d14a'
#superFitsDir        = '/Volumes/backupzeira/CALIFA/q043/px1/'
superFitsDir = CALIFAWorkDir + 'gal_fits/' + versionSuffix + '/'

#emLinesFitsDir      = CALIFAWorkDir + 'superfits/' + versionSuffix + '/'
emLinesFitsDir = CALIFAWorkDir + 'rgb-gas/' + versionSuffix + '/'
imgDir = CALIFAWorkDir + 'images/'

Zsun = 0.019

if __name__ == '__main__':
    t_init_gal = time.clock()
    args = parser_args()
    
    H = H5SFRData(args.hdf5)
    tSF__T = H.tSF__T
    
    debug = args.debug
    h5fname = args.hdf5
    galName = args.califaID
    iT = args.itSF
    
    if debug:
        print 'califaID: ', galName
        print 'h5fname: ', h5fname
        print 'iTSF: (%.2fMyr)' % (iT, tSF__T[iT] / 1e6)
    
    if (len(np.where(H.califaIDs == galName)[0]) == 0):
        exit('<<< plot: %s: no data.' % galName)
    
    # global
    xOkMin              = H.xOkMin
    tauVOkMin           = H.tauVOkMin
    tauVNebOkMin        = H.tauVNebOkMin
    tauVNebErrMax       = H.tauVNebErrMax
    
    # ALL gal
    dist_zone_HLR__g    = H.get_data_h5('dist_zone__g')
    ##stellar
    x_young__g          = H.get_data_h5('x_young__g')
    tau_V__g            = H.get_data_h5('tau_V__Tg')[iT]
    ##nebular
    EW_Hb__g            = H.get_data_h5('EW_Hb__g')
    EW_Ha__g            = H.get_data_h5('EW_Ha__g')
    tau_V_neb__g        = H.get_data_h5('tau_V_neb__g')
    
    # one gal
    dist_zone_HLR__z    = H.get_prop_gal('dist_zone__g', galName)
    ##stellar
    tau_V__z            = H.get_prop_gal('tau_V__Tg', galName)[iT]
    atau_V__r           = H.get_prop_gal('tau_V__Trg', galName)[iT]
    SFR__z              = H.get_prop_gal('SFR__Tg', galName)[iT]
    SFRSD__z            = H.get_prop_gal('SFRSD__Tg', galName)[iT]
    aSFRSD__r           = H.get_prop_gal('aSFRSD__Trg', galName)[iT]
    x_young__z          = H.get_prop_gal('x_young__g', galName)
    ##nebular
    EW_Ha__z            = H.get_prop_gal('EW_Ha__g', galName)
    EW_Hb__z            = H.get_prop_gal('EW_Hb__g', galName)
    tau_V_neb__z        = H.get_prop_gal('tau_V_neb__g', galName)
    atau_V_neb__r       = H.get_prop_gal('tau_V_neb__rg', galName)
    tau_V_neb_err__z    = H.get_prop_gal('tau_V_neb_err__g', galName)
    SFR_Ha__z           = H.get_prop_gal('SFR_Ha__g', galName)
    SFRSD_Ha__z         = H.get_prop_gal('SFRSD_Ha__g', galName)
    aSFRSD_Ha__r        = H.get_prop_gal('aSFRSD_Ha__rg', galName)
    
    #mixed
    delta_tau__z        = tau_V_neb__z - tau_V__z
    delta_tau__g        = tau_V_neb__g - tau_V__g
    EW_HaHb__z          = EW_Ha__z / EW_Hb__z
    EW_HaHb__g          = EW_Ha__g / EW_Hb__g
    
    mask_nuc__g         = dist_zone_HLR__g <= 0.5
    mask_nuc__z         = dist_zone_HLR__z <= 0.5  
    mask_disc__g        = (dist_zone_HLR__g > 0.5) & (dist_zone_HLR__g <= 2)
    mask_disc__z        = (dist_zone_HLR__z > 0.5) & (dist_zone_HLR__z <= 2)
    mask_off__g         = dist_zone_HLR__g > 2.
    mask_off__z         = dist_zone_HLR__z > 2.

    CALIFASuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
    CALIFAFitsFile = superFitsDir + galName + CALIFASuffix
    galaxyImgFile = imgDir + galName + '.jpg'
     
    K = fitsQ3DataCube(CALIFAFitsFile)
    #Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
     
    tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
     
    emLinesSuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
    emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
    K.loadEmLinesDataCube(emLinesFitsFile)
    
    #nebular 
    EW_HaHb__yx         = K.zoneToYX(EW_HaHb__z, extensive = False)
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # #stellar
    # tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
    # SFRSD__yx = K.zoneToYX(SFR__z, extensive = True)
    # x_young__yx = K.zoneToYX(x_young__z, extensive = False)
    # #nebular
    # tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
    # SFRSD_Ha__yx = K.zoneToYX(SFR_Ha__z, extensive = True)
    # EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
    # EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
    # tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z.data, extensive = False)
    # #mixed
    # deltaTau__z = tau_V_neb__z - tau_V__z 
    # deltaTau__yx = K.zoneToYX(deltaTau__z, extensive = False)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    t_calc = time.clock()
    print 'calc: elapsed time: %.2f' % (t_calc - t_init_gal)
    
    NRows = 2
    NCols = 3
    
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
    xlabel = r'$\tau_V^{neb}$'
    ylabel = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$'
    x = tau_V_neb__g
    y = np.ma.log10(EW_HaHb__g)
    mask = x.mask | y.mask | (EW_Ha__g < 1.) | (EW_Hb__g < 1.)
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlim = [tauVNebOkMin, 4.0]
    ylim = [0.25, 1.2]
    sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14, color = 'k')
    ##########################
    x = tau_V_neb__z
    y = np.ma.log10(EW_HaHb__z)
    # nuc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_nuc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'b', marker = 'o', s = 30, edgecolor = 'black', label = 'nucleus')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.08, 14, color = 'b')
    # disc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_disc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'r', marker = 'o', s = 30, edgecolor = 'black', label = 'disc')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.14, 14, color = 'r')
    # > 2HLR
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_off__z) 
    # xm = np.ma.masked_array(x, mask = mask)
    # ym = np.ma.masked_array(y, mask = mask)
    # sc = ax.scatter(xm, ym, c = 'g', marker = 'o', s = 30, edgecolor = 'black', label = '> 2HLR')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.125))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which = 'major')
    ax.legend(fontsize = 14, frameon = False, loc = 'upper left')

    ax = axArr[0, 2]
    ax.set_axis_on()
    xlabel = r'$\tau_V^{neb}\ -\ \tau_V$'
    ylabel = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$'
    x = delta_tau__g
    y = np.ma.log10(EW_HaHb__g)
    mask = x.mask | y.mask | (EW_Ha__g < 1.) | (EW_Hb__g < 1.)
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlim = [-1.5, 2.5]
    ylim = [0.25, 1.2]
    sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14, color = 'k')
    ##########################
    x = delta_tau__z
    y = np.ma.log10(EW_HaHb__z)
    # nuc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_nuc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'b', marker = 'o', s = 30, edgecolor = 'black', label = 'nucleus')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.08, 14, color = 'b')
    # disc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_disc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'r', marker = 'o', s = 30, edgecolor = 'black', label = 'disc')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.14, 14, color = 'r')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # # > 2HLR
    # mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_off__z) 
    # xm = np.ma.masked_array(x, mask = mask)
    # ym = np.ma.masked_array(y, mask = mask)
    # sc = ax.scatter(xm, ym, c = 'g', marker = 'o', s = 30, edgecolor = 'black', label = '> 2HLR')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.125))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which = 'major')
    #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')

    ax = axArr[1, 0]
    ax.set_axis_on()
    label = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$'
    lim = [-0.25, 1.2]
    vmin = lim[0]
    vmax = lim[1]
    im = ax.imshow(np.ma.log10(EW_HaHb__yx), origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin, vmax = vmax)
    #im = ax.imshow(np.ma.log10(EW_HaHb__yx), origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
    DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
    ax.set_title(label, y = -0.15)
    ax.grid()
    f.colorbar(ax = ax, mappable = im, use_gridspec = True)

    ax = axArr[1, 1]
    ax.set_axis_on()
    xlabel = r'$\tau_V$'
    ylabel = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$'
    x = tau_V__g
    y = np.ma.log10(EW_HaHb__g)
    mask = x.mask | y.mask | (EW_Ha__g < 1.) | (EW_Hb__g < 1.)
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlim = [tauVOkMin, 2.5]
    ylim = [-0.75, 1.2]
    sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14, color = 'k')
    ##########################
    x = tau_V__z
    y = np.ma.log10(EW_HaHb__z)
    # nuc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_nuc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'b', marker = 'o', s = 30, edgecolor = 'black', label = 'nucleus')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.08, 14, color = 'b')
    # disc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_disc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'r', marker = 'o', s = 30, edgecolor = 'black', label = 'disc')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.14, 14, color = 'r')
    # > 2HLR
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_off__z) 
    # xm = np.ma.masked_array(x, mask = mask)
    # ym = np.ma.masked_array(y, mask = mask)
    # sc = ax.scatter(xm, ym, c = 'g', marker = 'o', s = 30, edgecolor = 'black', label = '> 2HLR')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which = 'major')
    #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')

    ax = axArr[1, 2]
    ax.set_axis_on()
    xlabel = r'$x_Y [\%]$'
    ylabel = r'$\log\ \frac{W_{H\alpha}}{W_{H\beta}}$'
    x = x_young__g * 100.
    y = np.ma.log10(EW_HaHb__g)
    mask = x.mask | y.mask | (EW_Ha__g < 1.) | (EW_Hb__g < 1.)
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    xlim = [0,100.]
    ylim = [-0.75, 1.2]
    sc = ax.scatter(xm, ym, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
    nBox = 20
    dxBox = (xm.max() - xm.min()) / (nBox - 1.)
    X = x[~mask]
    Y = y[~mask]
    aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
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
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14, color = 'k')
    ##########################
    x = x_young__z * 100.
    y = np.ma.log10(EW_HaHb__z)
    # nuc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_nuc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'b', marker = 'o', s = 30, edgecolor = 'black', label = 'nucleus')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.08, 14, color = 'b')
    # disc
    mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_disc__z) 
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    sc = ax.scatter(xm, ym, c = 'r', marker = 'o', s = 30, edgecolor = 'black', label = 'disc')
    a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.14, 14, color = 'r')
    # > 2HLR
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask = x.mask | y.mask | (EW_Ha__z < 1.) | (EW_Hb__z < 1.) | ~(mask_off__z) 
    # xm = np.ma.masked_array(x, mask = mask)
    # ym = np.ma.masked_array(y, mask = mask)
    # sc = ax.scatter(xm, ym, c = 'g', marker = 'o', s = 30, edgecolor = 'black', label = '> 2HLR')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ##########################
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.grid(which = 'major')
    #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')

    f.suptitle(r'%s - morph:%s  tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (galName, tipos, (age / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax))
    f.subplots_adjust(left = 0.07, bottom = 0.1, right = 0.95, wspace = 0.3, top = 0.9)
    #f.subplots_adjust(left = 0.07, bottom = 0.1, right = 0.95, wspace = 0.15, top = 0.9)
    f.savefig('%s_EW_mosaic.png' % galName)
    plt.close(f)
    t_plot = time.clock()
    print 'plot: elapsed time: %.2f' % (t_plot - t_calc)
    print 'total: elapsed time: %.2f' % (t_plot - t_init_gal)

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     ax = axArr[0, 2]
#     ax.set_axis_on()
#     xlabel = r'EW(H$\beta$) [$\AA$]'
#     im = ax.imshow(EW_Hb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[0, 3]
#     ax.set_axis_on()
#     xlabel = r'$\epsilon\tau_V^{neb}$'
#     im = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 1)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[1, 0]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V__rg.flatten())
#     y = np.ma.log10(aSFRSD__rg.flatten() * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     xlabel = r'$\log\ \tau_V^{\star}(R)$'
#     ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     xlim = [np.log10(tauVOkMin), 0.5]
#     ylim = [-3.5, 1]
#     sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
#     nBox = 20
#     dxBox = (xm.max() - xm.min()) / (nBox - 1.)
#     X = x[~mask]
#     Y = y[~mask]
#     aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
#     xbinCenter = aux[0]
#     xMedian = aux[1]
#     xMean = aux[2]
#     xStd = aux[3]
#     yMedian = aux[4]
#     yMean = aux[5]
#     yStd = aux[6]
#     nInBin = aux[7]
#     xPrc = aux[8]
#     yPrc = aux[9]
#     ax.plot(xMedian, yMedian, 'k', lw = 2)
#     ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
#     ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
#     txt = '%.2f Myr' % (age / 1e6)
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14)
#     ##########################
#     x = np.ma.log10(atau_V__r)
#     y = np.ma.log10(aSFRSD__r * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # vmax = xr.mean() + 2. * xr.std()
#     # vmin = xr.mean() - 2. * xr.std()
#     # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
#     cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
#     cb.set_label(r'R [HLR]')
#     ##########################
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.xaxis.set_major_locator(MultipleLocator(0.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.yaxis.set_major_locator(MultipleLocator(0.5))
#     ax.yaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.grid(which = 'major')
#     
#     sigma_dev = 3.
# 
#     ax = axArr[1, 1]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V__yx)
#     y = np.ma.log10(SFRSD__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     ylim = [-3.5, 1]
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # vmax = y.mean() + sigma_dev * y.std()
#     # vmin = y.mean() - sigma_dev * y.std()
#     # im = ax.imshow(ym, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin, vmax = vmax)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     txt = 'not masked: %d' % (~mask).sum()
#     print (~mask).sum()
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     
#     ax = axArr[1, 2]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V__yx)
#     y = np.ma.log10(SFRSD__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     xlabel = r'$\log\ \tau_V^\star$' 
#     ylim = [-3.5, 1]
#     vmin = np.log10(tauVOkMin)
#     im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
#     #im = ax.imshow(x, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[1, 3]
#     ax.set_axis_on()
#     xlabel = r'$x_Y [\%]$'
#     x_young__yx = K.zoneToYX(x_young__z, extensive = False)
#     im = ax.imshow(100. * x_young__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(xlabel, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[2, 0]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__rg.flatten())
#     y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     xlabel = r'$\log\ \tau_V^{neb}(R)$'
#     ylabel = r'$\log\ \langle \Sigma_{SFR}^{neb}(R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
#     xlim = [np.log10(tauVNebOkMin), 0.5]
#     ylim = [-3.5, 1]
#     sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
#     nBox = 20
#     dxBox = (xm.max() - xm.min()) / (nBox - 1.)
#     X = x[~mask]
#     Y = y[~mask]
#     aux = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
#     xbinCenter = aux[0]
#     xMedian = aux[1]
#     xMean = aux[2]
#     xStd = aux[3]
#     yMedian = aux[4]
#     yMean = aux[5]
#     yStd = aux[6]
#     nInBin = aux[7]
#     xPrc = aux[8]
#     yPrc = aux[9]
#     ax.plot(xMedian, yMedian, 'k', lw = 2)
#     ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
#     ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
#     txt = '%.2f Myr' % (age / 1e6)
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     a_ols, b_ols = plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.02, 14)
#     ##########################
#     x = np.ma.log10(atau_V_neb__r)
#     y = np.ma.log10(aSFRSD_Ha__r * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # vmax = xr.mean() + 2. * xr.std()
#     # vmin = xr.mean() - 2. * xr.std()
#     # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
#     cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
#     cb.set_label(r'R [HLR]')
#     ##########################
#     ax.set_xlabel(xlabel)
#     ax.set_ylabel(ylabel)
#     ax.set_xlim(xlim)
#     ax.set_ylim(ylim)
#     ax.xaxis.set_major_locator(MultipleLocator(0.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.yaxis.set_major_locator(MultipleLocator(0.5))
#     ax.yaxis.set_minor_locator(MultipleLocator(0.125))
#     ax.grid(which = 'major')
# 
#     ax = axArr[2, 1]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__yx)
#     y = np.ma.log10(SFRSD_Ha__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     label = r'$\log\ \Sigma_{SFR}^{neb} [M_\odot yr^{-1} kpc^{-2}]$' 
#     ylim = [-3.5, 1]
#     im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     txt = 'not masked: %d' % (~mask).sum()
#     print (~mask).sum() 
#     plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     
#     ax = axArr[2, 2]
#     ax.set_axis_on()
#     x = np.ma.log10(tau_V_neb__yx)
#     y = np.ma.log10(SFRSD_Ha__yx * 1e6)
#     mask = x.mask | y.mask
#     xm = np.ma.masked_array(x, mask = mask)
#     ym = np.ma.masked_array(y, mask = mask)
#     label = r'$\log\ \tau_V^{neb}$' 
#     vmin = np.log10(tauVNebOkMin)
#     im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     ax = axArr[2, 3]
#     ax.set_axis_on()
#     label = r'$\delta\ \tau_V$' 
#     im = ax.imshow(deltaTau__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     ax.set_title(label, y = -0.15)
#     ax.grid()
#     f.colorbar(ax = ax, mappable = im, use_gridspec = True)
# 
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # x = np.ma.log10(tau_V__Tz[iT])
#     # y = np.ma.log10(SFRSD__Tz[iT] * 1e6)
#     # mask = x.mask | y.mask
#     # xm = np.ma.masked_array(x, mask = mask)
#     # ym = np.ma.masked_array(y, mask = mask)
#     # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
#     # xr__yx = K.zoneToYX(xr, extensive = False)
#     # xlabel = r'xr'
#     # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # # vmax = xr__yx.mean() + 2. * xr__yx.std()
#     # # vmin = xr__yx.mean() - 2. * xr__yx.std()
#     # # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
#     # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
#     # DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
#     # txt = 'Nz: %d' % K.N_zone 
#     # plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
#     # txt = 'masked: %d' % mask.sum() 
#     # plot_text_ax(ax, txt, 0.02, 0.92, 14, 'top', 'left')
#     # ax.set_title(xlabel, y=-0.15)
#     # ax.grid()
#     # f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#     #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

