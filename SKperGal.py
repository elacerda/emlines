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
    'underS06' : False,
    'hdf5' : None,
    'plot' : False,
    'weiradprof' : True,
    'minpopx' : 0.05,
    'mintauv' : 0.05,
    'mintauvneb' : 0.05,
    'maxtauvneberr' : 999.,
    'rbinini' : 0.0,
    'rbinfin' : 2.0,
    'rbinstep' : 0.1,
    'califaID' : 'K0277'
}

def parser_args():
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--underS06', 
                        action = 'store_true',
                        default = default['underS06'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default['weiradprof'])
    parser.add_argument('--plot', '-P',
                        action = 'store_true',
                        default = default['plot'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--califaID', '-g',
                        metavar = 'FILE',
                        type = str,
                        default = default['califaID'])
    parser.add_argument('--minpopx', 
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default['minpopx'])
    parser.add_argument('--mintauv', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['mintauv'])
    parser.add_argument('--mintauvneb', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['mintauvneb'])
    parser.add_argument('--maxtauvneberr', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['maxtauvneberr'])
    parser.add_argument('--rbinini', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['rbinini'])
    parser.add_argument('--rbinfin', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['rbinfin'])
    parser.add_argument('--rbinstep', 
                        metavar = 'FRAC',
                        type = float,
                        default = default['rbinstep'])

    return parser.parse_args()

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

args = parser_args()

debug = args.debug
BPTLowS06 = args.underS06
h5fname = args.hdf5
weiRadProf = args.weiradprof
xOkMin = args.minpopx
tauVOkMin = args.mintauv
tauVNebOkMin = args.mintauvneb
tauVNebErrMax = args.maxtauvneberr
RbinIni = args.rbinini
RbinFin = args.rbinfin
RbinStep = args.rbinstep
galName = args.califaID

if debug:
    print 'califaID: ', galName
    print 'BPTLowS06: ', BPTLowS06
    print 'h5fname: ', h5fname
    print 'weiRadProf: ', weiRadProf
    print 'xOkMin: ', xOkMin
    print 'tauVOkMin: ', tauVOkMin
    print 'tauVNebOkMin: ', tauVNebOkMin
    print 'tauVNebErrMax: ', tauVNebErrMax
    print 'RbinIni: ', RbinIni
    print 'RbinFin: ', RbinFin
    print 'RbinStep: ', RbinStep

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
qCCM = {
    '4861' : 1.16427,
    '5007' : 1.12022,
    '6563' : 0.81775,
    '6583' : 0.81466,
}

Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins = len(RbinCenter__r)

RColor = [ 'r', 'y', 'b', 'k' ]
RRange = [  .5, 1., 1.5, 2.  ]

# SFR-time-scale array (index __T)
base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
tSF__T = base.ageBase
N_T = base.nAges

if __name__ == '__main__':
    t_init_gal = time.clock()
    
    CALIFASuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
    CALIFAFitsFile = superFitsDir + galName + CALIFASuffix
    emLinesSuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
    emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
    galaxyImgFile = imgDir + galName + '.jpg'
    
    # both files
    if not (os.path.isfile(CALIFAFitsFile) and os.path.isfile(emLinesFitsFile)):
        exit('<<< %s: miss files' % galName)
    
    K = fitsQ3DataCube(CALIFAFitsFile)
    tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
    
    # read FITSFILE containing galaxy emission lines measured by R.G.B.
    # read_rgb_fits returns False if emLinesFitsFile does not exists.
    #read = read_rgb_fits(emLinesFitsFile, read_lines)
    K.loadEmLinesDataCube(emLinesFitsFile)
    
    # Problem in FITS file
    if K.EL.flux[0, :].sum() == 0.:
        exit('<<< %s: EmLines FITS problem' % galName)
    
    # Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
    
    ##########################
    ###### MASK EmLines ######
    ##########################        
    # minimum value of f_lz / err_f_lz
    minSNR = 3.
    
    i_Hb = K.EL.lines.index('4861')
    i_O3 = K.EL.lines.index('5007')
    i_Ha = K.EL.lines.index('6563')
    i_N2 = K.EL.lines.index('6583')
    
    Ha = K.EL.flux[i_Ha, :]
    eHa = K.EL.eflux[i_Ha, :]
    Hb = K.EL.flux[i_Hb, :]
    eHb = K.EL.eflux[i_Hb, :]
    O3 = K.EL.flux[i_O3, :]
    eO3 = K.EL.eflux[i_O3, :]
    N2 = K.EL.flux[i_N2, :]
    eN2 = K.EL.eflux[i_N2, :]
    
    HbOk = np.array((Hb / eHb) >= minSNR, dtype = np.bool)
    O3Ok = np.array((O3 / eO3) >= minSNR, dtype = np.bool)
    HaOk = np.array((Ha / eHa) >= minSNR, dtype = np.bool)
    N2Ok = np.array((N2 / eN2) >= minSNR, dtype = np.bool)
    
    maskLinesSNROk__z = HbOk & O3Ok & HaOk & N2Ok
    maskFluxOk__z = (Ha >= 0) & (Hb >= 0) & (O3 >= 0) & (N2 >= 0)
    
    print 'Zones masked by SNR cut: %d' % (~maskLinesSNROk__z).sum()
    print 'Zones masked by flux positive cut: %d' % (~maskFluxOk__z).sum()
    
    maskBPT__z = None
    
    if BPTLowS06:
        L = Lines()
        N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
        O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
        maskBPT__z = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
    ##########################
    ##########################
    ##########################
            
    ##########################
    ####### STARLIGHT ########
    ##########################        
    x__tZz = K.popx / K.popx.sum(axis = 1).sum(axis = 0)
    
    # pop_young
    age_young_max = 1e8
    flag_young__t = K.ageBase <= age_young_max
    x_young__z = x__tZz[flag_young__t, :].sum(axis = 1).sum(axis = 0)
    
    tau_V__Tz   = np.ma.zeros((N_T, K.N_zone))
    SFR__Tz     = np.ma.zeros((N_T, K.N_zone))
    SFRSD__Tz   = np.ma.zeros((N_T, K.N_zone))
    
    SFRSD__Tyx  = np.ma.zeros((N_T, K.N_y, K.N_x))
    tau_V__Tyx  = np.ma.zeros((N_T, K.N_y, K.N_x))
    
    atau_V__Tr  = np.ma.zeros((N_T, NRbins))
    aSFRSD__Tr  = np.ma.zeros((N_T, NRbins))
    npts__Tr    = np.ma.zeros((N_T, NRbins))
    
    integrated_SFR__T   = np.ma.zeros((N_T))
    integrated_SFRSD__T = np.ma.zeros((N_T))
    
    
    for iT, tSF in enumerate(tSF__T):
        flag__t = K.ageBase <= tSF

        tau_V__z = K.tau_V__z
        
        if xOkMin >= 0.:
            # Compute xOk "raw" image
            xOk__z = x__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
            mask_notOk__z = (xOk__z < xOkMin) | (tau_V__z < tauVOkMin)
            print 'age %.2fMyr: Zones masked by min pop young cut: %d' % (tSF / 1e6, mask_notOk__z.sum())
            
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # if tauVOkMin >= 0:
        #     mask_notOk__z |= (tau_V__Tz[iT] < tauVOkMin)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        
        SFR__z = K.Mini__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0) / tSF
        SFRSD__z = SFR__z / K.zoneArea_pc2

        tau_V__Tz[iT]   = np.ma.masked_array(tau_V__z, mask = mask_notOk__z, dtype = K.tau_V__z.dtype)
        SFR__Tz[iT]     = np.ma.masked_array(SFR__z, mask = mask_notOk__z)
        SFRSD__Tz[iT]   = np.ma.masked_array(SFRSD__z, mask = mask_notOk__z)

        tau_V__Tyx[iT] = K.zoneToYX(tau_V__Tz[iT], extensive = False)
        SFRSD__Tyx[iT] = K.zoneToYX(SFR__Tz[iT], extensive = True)
            
        integrated_SFR__T = SFR__Tz[iT].sum()
        integrated_SFRSD__T = integrated_SFR__T / K.zoneArea_pc2.sum()
                    
        atau_V__Tr[iT], npts__Tr[iT] = K.radialProfile(tau_V__Tyx[iT], Rbin__r, rad_scale = K.HLR_pix, return_npts=True)
        aSFRSD__Tr[iT] = K.radialProfile(SFRSD__Tyx[iT], Rbin__r, rad_scale = K.HLR_pix)
    ##########################
    ##########################
    ##########################    
    
    ##########################
    ########## tau_V #########
    ##########################
    mask_tmp = maskFluxOk__z & maskLinesSNROk__z
    
    print 'union of lines and negative fluxes masks: %d' % mask_tmp.sum()
    
    K.EL._forceMask = ~mask_tmp #Changing global EL mask
    
    if BPTLowS06:
        K.EL._forceMask = ~(mask_tmp & maskBPT__z) #Changing global EL mask
    
    maskOkTauVNeb = (K.EL.tau_V_neb__z >= tauVNebOkMin)
    print 'Zones masked by min tauVNeb cut: %d' % (~maskOkTauVNeb).sum()
    maskOkTauVNeb_err = (K.EL.tau_V_neb_err__z <= tauVNebErrMax)
    print 'Zones masked by min tauVNeb cut: %d' % (~maskOkTauVNeb_err).sum()
    
    K.EL._forceMask = ~(mask_tmp & maskOkTauVNeb & maskOkTauVNeb_err) #Changing global EL mask
    print 'union of lines, negative fluxes and tauVNeb masks: %d' % K.EL._forceMask.sum() 
    
    N_zones_tau_V = len(K.EL.tau_V_neb__z.compressed())
    
    if (N_zones_tau_V == 0):
        exit('<<< %s: tau_V_neb for %d zones' % (galName, N_zones_tau_V))
    
    tau_V_neb__z = K.EL.tau_V_neb__z
    tau_V_neb_err__z = K.EL.tau_V_neb_err__z
    tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
    tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
    atau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
    
    ##########################
    #### intrinsic Ha Lum ####
    ##########################
    L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
    L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun
    
    # L_int_Ha__Lz intrinsic Ha luminosity 
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    eHa = np.ma.exp(qCCM['6563'] * tau_V_neb__z)
    L_obs_HaHb__z = L_obs__Lz[i_Ha, :] / L_obs__Lz[i_Hb, :]
    L_int_Ha__z = L_obs__Lz[i_Ha, :] * eHa
    
    integrated_eHa = np.ma.exp(K.EL._qCCM['6563'] * K.EL.integrated_tau_V_neb)
    integrated_L_obs_Ha__L = K.EL._F_to_L(K.EL.integrated_flux) / L_sun
    integrated_L_int_Ha = integrated_L_obs_Ha__L[i_Ha] * integrated_eHa
    
    # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
    a = L_obs_err__Lz[i_Ha, :]
    b = q * L_obs_HaHb__z * L_obs_err__Lz[i_Hb, :]
    L_int_Ha_err__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
    ##########################
    ##########################
    ##########################
    
    ##########################
    #### SFR and SigmaSFR ####
    ##########################
    # 3.17 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter        
    SFR_Ha__z = 3.13 * L_int_Ha__z / (1.e8)
    SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
    
    integrated_SFR_Ha = 3.13 * integrated_L_int_Ha / (1.e8)
    integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()
    
    SFRSD_Ha__yx = K.zoneToYX(SFR_Ha__z, extensive = True)
    aSFRSD_Ha__r = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
    ##########################
    ##########################
    ##########################
    
    t_calc = time.clock()
    print 'calc: elapsed time: %.2f' % (t_calc - t_init_gal)
    
    if args.plot == True:
        H = H5SFRData(h5fname)
        
        if (len(np.where(H.califaIDs == galName)[0]) == 0):
            exit('<<< plot: %s: no data.' % galName)
            
        aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
        aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
        tau_V__Trg = H.get_data_h5('tau_V__Trg')
        tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')

        tmp = np.copy(K.EL._forceMask)
        K.EL._forceMask = np.zeros_like(tmp, dtype = np.bool)
        EW_Ha__z = K.EL.EW[i_Ha, :]
        EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
        EW_Hb__z = K.EL.EW[i_Hb, :]
        EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
        tauVNeb_err__z = K.EL.tau_V_neb_err__z
        tauVNeb_err__yx = K.zoneToYX(tauVNeb_err__z, extensive = False)
        K.EL._forceMask = np.copy(tmp)

        NRows = 3
        NCols = 4
        
        f, axArr = plt.subplots(NRows, NCols)
        f.set_size_inches((NCols * 5.34, NRows * 5.))
        
        for ax in f.axes:
            ax.set_axis_off()
        
        iT = 11
        age = tSF__T[iT]

        ax = axArr[0, 0]
        ax.set_axis_on()
        galimg = plt.imread(galaxyImgFile)[::-1,:,:]
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.imshow(galimg, origin = 'lower')
        DrawHLRCircleInSDSSImage(ax, K.HLR_pix, pa, ba)
        
        ax = axArr[0, 1]
        ax.set_axis_on()
        xlabel = r'EW(H$\alpha$) [$\AA$]'
        im = ax.imshow(EW_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(xlabel, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)

        ax = axArr[0, 2]
        ax.set_axis_on()
        xlabel = r'EW(H$\beta$) [$\AA$]'
        im = ax.imshow(EW_Hb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(xlabel, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)

        ax = axArr[0, 3]
        ax.set_axis_on()
        xlabel = r'$\epsilon\tau_V^{neb}$'
        im = ax.imshow(tauVNeb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 1)
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(xlabel, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)

        ax = axArr[1, 0]
        ax.set_axis_on()
        x = np.ma.log10(tau_V__Trg[iT].flatten())
        y = np.ma.log10(aSFRSD__Trg[iT].flatten() * 1e6)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xlabel = r'$\log\ \tau_V^{\star}(R)$'
        ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
        xlim = [-2., 0.5]
        ylim = [-3.5, 1]
        sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
        nBox = 20
        dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
        X = x[~mask]
        Y = y[~mask]
        aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
        xbinCenter  = aux[0]
        xMedian     = aux[1]
        xMean       = aux[2]
        xStd        = aux[3]
        yMedian     = aux[4]
        yMean       = aux[5]
        yStd        = aux[6]
        nInBin      = aux[7]
        xPrc        = aux[8]
        yPrc        = aux[9]
        ax.plot(xMedian, yMedian, 'k', lw = 2)
        ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
        ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
        txt = '%.2f Myr' % (age / 1e6)
        plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        a_ols, b_ols = plotOLSbisectorAxis(ax, X, Y, 0.98, 0.02, 14)
        ##########################
        x = np.ma.log10(atau_V__Tr[iT])
        y = np.ma.log10(aSFRSD__Tr[iT] * 1e6)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
        # vmax = xr.mean() + 2. * xr.std()
        # vmin = xr.mean() - 2. * xr.std()
        # ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        sc = ax.scatter(xm, ym, c = RbinCenter__r, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = RbinFin, vmin = RbinIni)
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
        x = np.ma.log10(tau_V__Tyx[iT])
        y = np.ma.log10(SFRSD__Tyx[iT] * 1e6)
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
        print (~mask).sum()
        plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        ax.set_title(xlabel, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)
        
        ax = axArr[1, 2]
        ax.set_axis_on()
        x = np.ma.log10(tau_V__Tyx[iT])
        y = np.ma.log10(SFRSD__Tyx[iT] * 1e6)
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        xlabel = r'$\log\ \tau_V^\star$' 
        ylim = [-3.5, 1]
        vmin = np.log10(tauVOkMin)
        im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
        #im = ax.imshow(x, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(xlabel, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)

        ax = axArr[1, 3]
        ax.set_axis_on()
        xlabel = r'$x_Y [\%]$'
        x_young__yx = K.zoneToYX(x_young__z, extensive = False)
        im = ax.imshow(100. * x_young__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(xlabel, y=-0.15)
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
        xlim = [-2., 0.5]
        ylim = [-3.5, 1]
        sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
        nBox = 20
        dxBox       = (xm.max() - xm.min()) / (nBox - 1.)
        X = x[~mask]
        Y = y[~mask]
        aux         = calcRunningStats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
        xbinCenter  = aux[0]
        xMedian     = aux[1]
        xMean       = aux[2]
        xStd        = aux[3]
        yMedian     = aux[4]
        yMean       = aux[5]
        yStd        = aux[6]
        nInBin      = aux[7]
        xPrc        = aux[8]
        yPrc        = aux[9]
        ax.plot(xMedian, yMedian, 'k', lw = 2)
        ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
        ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
        txt = '%.2f Myr' % (age / 1e6)
        plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        a_ols, b_ols = plotOLSbisectorAxis(ax, X, Y, 0.98, 0.02, 14)
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
        sc = ax.scatter(xm, ym, c = RbinCenter__r, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = RbinFin, vmin = RbinIni)
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
        print (~mask).sum() 
        plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
        ax.set_title(label, y=-0.15)
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
        ax.set_title(label, y=-0.15)
        ax.grid()
        f.colorbar(ax = ax, mappable = im, use_gridspec = True)

        ax = axArr[2, 3]
        ax.set_axis_on()
        x = tau_V_neb__z
        y = tau_V__Tz[iT]
        mask = x.mask | y.mask
        xm = np.ma.masked_array(x, mask = mask)
        ym = np.ma.masked_array(y, mask = mask)
        deltaTau__z = xm - ym
        deltaTau__yx = K.zoneToYX(deltaTau__z, extensive = False)
        label = r'$\delta\ \tau_V$' 
        im = ax.imshow(deltaTau__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
        DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
        ax.set_title(label, y=-0.15)
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

        f.suptitle(r'%s - morph:%s  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (galName, tipos, xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax))
        f.subplots_adjust(left=0.07, bottom=0.1, right=0.99, wspace = 0.1, top=0.9)
        f.savefig('%s_mosaic.png' % galName)
        plt.close(f)
        t_plot = time.clock()
        print 'plot: elapsed time: %.2f' % (t_plot - t_calc)
        print 'total: elapsed time: %.2f' % (t_plot - t_init_gal)