#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
from rgbread import read_rgb_fits
from morph_type import get_morph, morph_number
import numpy as np
from pycasso import fitsQ3DataCube
import pyfits
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
from get_morfologia import get_morfologia
import os

plot = True
#plot = False
debug = False
#debug = True

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size']       = 16
mpl.rcParams['axes.labelsize']  = 22
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['font.family']     = 'sans-serif'

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
    
galaxiesListFile    = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
baseCode            = 'Bgsd6e'
versionSuffix       = 'px1_q043.d14a'
#versionSuffix       = 'v20_q043.d14a'
superFitsDir        = '/Volumes/backupzeira/CALIFA/q043/px1/'
#superFitsDir        = CALIFAWorkDir + 'gal_fits/' + versionSuffix + '/'

#emLinesFitsDir      = CALIFAWorkDir + 'superfits/' + versionSuffix + '/'
emLinesFitsDir      = CALIFAWorkDir + 'rgb-gas/' + versionSuffix + '/'
imgDir              = CALIFAWorkDir + 'images/'

Zsun = 0.019
Lsun = 3.826e33 # erg/s
qCCM = {
    '4861' : 1.16427,
    '5007' : 1.12022,
    '6563' : 0.81775,
    '6583' : 0.81466,
}

read_lines = ['4861', '5007', '6563', '6583']        
read_lines = None        

f               = open(galaxiesListFile, 'r')
listOfPrefixes  = f.readlines()
f.close()

if debug:
    #listOfPrefixes = listOfPrefixes[0:20]        # Q&D tests ...
    listOfPrefixes = ['K0026\n']
    
N_gals = len(listOfPrefixes)

# SFR-time-scale array (index __T)
tSF__T = np.array([ 10.01 , 25.2 , 63.2, 100.1 , 158.6 , 199.6 , 1400.2 ]) * 1.e6
N_T = len(tSF__T)

mask_xOk = True

# Def smallest light fraction (in the flag__t-ageMax age-range) deemed to be Ok for our stats ...
xOkMin = 0.05

# Minimum tauV to be taken seriously ...
tauVOkMin = 0.05

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

def radialProfileWeighted(v__yx, w__yx, bins, rad_scale, func_radialProfile = None):
    v__r = None

    if func_radialProfile:
        w__r = func_radialProfile(w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v_w__r = func_radialProfile(v__yx * w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v__r = v_w__r / w__r

    return v__r

def Mpc_to_cm(dist):
    # google: 1 Mpc = 3.08567758e24 cm
    return dist * 3.08567758e24 

def tauV_to_AV(tauV, err_tauV):
    c = np.log10(np.exp(1)) / 0.4
    AV = c * tauV
    err_AV = c * err_tauV
    
    return AV, err_AV

def flux_to_LSol(flux, distance):
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / Lsun

def calc_Lobs(f_obs__Lz, distance_Mpc):
    '''
    Calculate luminosity using 
    L_\lambda = 4 \pi d^2 F_\lambda
    Distance in Mpc due to flux_to_LSol() uses this unit.
    
    Lacerda@Saco - 25/Jun/2014
    '''
    solidAngle = 4. * np.pi * distance_Mpc
    
    Lobs__Lz        = np.ma.zeros(f_obs__Lz.shape)
    err_Lobs__Lz    = np.ma.zeros(f_obs__Lz.shape)
    
    for line in range(f_obs__Lz.shape[0]):
        Lobs__Lz[line]      = flux_to_LSol(f_obs__Lz[line, :], K.distance_Mpc)
        err_Lobs__Lz[line]  = flux_to_LSol(err_f_obs__Lz[line, :], K.distance_Mpc)
        
    return Lobs__Lz, err_Lobs__Lz

def calc_tauVNeb(K, f_obs__Lz, err_f_obs__Lz, lines):
    '''
    Calculate Balmer optical depth (tau_V).
    
    Lacerda@Saco - 25/Jun/2014
    '''
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    f_int_HaHb = 2.86
    f_obs_HaHb__z = f_obs__Lz[i_Ha, :] / f_obs__Lz[i_Hb, :]
    q = qCCM['4861'] - qCCM['6563']
    #LintHaHb = 2.86
    #LobsHaHb = Lobs['6563'] / Lobs['4861'] # OBS: LHa / LHb = fluxHa / fluxHb 
    
    #tauVNeb__z = np.ma.log(LobsHaHb / LintHaHb) / q
    lnHaHb = np.ma.log(f_obs_HaHb__z / f_int_HaHb)
    tauVNeb__z = lnHaHb / q 
    
    #err_tauVNeb__z = np.sqrt((err_Lobs['6563'] / Lobs['6563']) ** 2.0 + (err_Lobs['4861'] / Lobs['4861']) ** 2.0) / (qCCM['4861'] - qCCM['6563'])
    a = err_f_obs__Lz[i_Ha, :] / f_obs__Lz[i_Ha, :]
    b = err_f_obs__Lz[i_Hb, :] / f_obs__Lz[i_Hb, :]
    err_tauVNeb__z = np.sqrt(a ** 2.0 + b ** 2.0) / q
    
    a                 = err_f_obs__Lz[i_Ha, :]
    b                 = (f_obs__Lz[i_Ha, :] / f_obs__Lz[i_Hb, :]) * err_f_obs__Lz[i_Hb, :]
    err_f_obs_HaHb__z  = np.ma.sqrt(a ** 2. + b ** 2.) / f_obs__Lz[i_Hb, :]
        
    return tauVNeb__z, err_tauVNeb__z, f_obs_HaHb__z, err_f_obs_HaHb__z

def calc_Lint_Ha(Lobs__Lz, err_Lobs__Lz, tauVNeb__z,lines):
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    
    eHa = np.ma.exp(qCCM['6563'] * tauVNeb__z)
    LobsHaHb = Lobs__Lz[i_Ha, :] / Lobs__Lz[i_Hb, :]

    Lint_Ha__z = Lobs__Lz[i_Ha, :] * eHa
    
    a = err_Lobs__Lz[i_Ha, :]
    b = q * LobsHaHb * err_Lobs__Lz[i_Hb, :]
    err_Lint_Ha__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
    
    return Lint_Ha__z, err_Lint_Ha__z

def calc_logZNeb(K, f_obs__Lz, err_f_obs__Lz, tauVNeb__z, lines):
    '''
    Calculates Z_Neb using Asari et al (2007)
    Z_Neb = log((O/H) / (O/H)_Sol) = - 0.14 - 0.25 log([OIII]5007 / [NII]6583)
    
    O3 and N2 are nebular dust corrected using Balmer optical depth (tauVNeb).
    
    [OIII]5007_corrected = [OIII]5007_obs * exp(tauVNeb * q_[OIII]5007)
    [NII]6583_corrected  = [NII]6583_obs * exp(tauVNeb * q_[NII]6583)
    Lacerda@Saco - 25/Jun/2014
    '''
    i_O3 = lines.index('5007')
    i_N2 = lines.index('6583')
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')

    O3cor__z        = f_obs__Lz[i_O3, :] * np.ma.exp(qCCM['5007'] * tauVNeb__z)
    N2cor__z        = f_obs__Lz[i_N2, :] * np.ma.exp(qCCM['4861'] * tauVNeb__z)
    O3N2__z         = O3cor__z / N2cor__z
    
    e               = np.ma.exp(tauVNeb__z * (qCCM['5007'] - qCCM['6583'])) / f_obs__Lz[i_N2, :]
    a               = err_f_obs__Lz[i_O3, :]
    b               = (f_obs__Lz[i_O3, :] / f_obs__Lz[i_N2, :]) * err_f_obs__Lz[i_N2, :]
    q               = (qCCM['5007'] - qCCM['6583']) / (qCCM['4861'] - qCCM['6563'])
    c               = (f_obs__Lz[i_O3, :] / f_obs__Lz[i_Ha, :]) * err_f_obs__Lz[i_Ha, :]
    d               = (f_obs__Lz[i_O3, :] / f_obs__Lz[i_Hb, :]) * err_f_obs__Lz[i_Hb, :]
    err_O3N2__z     = e * np.ma.sqrt(a ** 2. + b ** 2. + q ** 2. * (c ** 2.0 + d ** 2.))

    # TODO: remember to correct the errors O3cor and N2cor 
    logZNeb__z      = - 0.14 - (0.25 * np.log10(O3N2__z))
    err_logZNeb__z  = 0.25 * err_O3N2__z / (np.log(10.) * O3N2__z)

    return logZNeb__z, err_logZNeb__z, O3N2__z, err_O3N2__z 
    
def calcRunningStats(x,y, dxBox=0.3, xbinIni=8.5, xbinFin=12, xbinStep=0.05):
    '''
    Statistics of x & y in equal size x-bins (dx-box).
    Note the mery small default xbinStep, so we have overlaping boxes.. so running stats..

    Cid@Lagoa -
    '''

    # Def x-bins
    xbin = np.arange(xbinIni, xbinFin + xbinStep, xbinStep)
    xbinCenter = (xbin[:-1] + xbin[1:]) / 2.0
    Nbins = len(xbinCenter)

    # Reset in-bin stats arrays
    xMedian , xMean , xStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    yMedian , yMean , yStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    nInBin = np.zeros(Nbins)

    # fill up in x & y stats for each x-bin
    for ixBin in np.arange(Nbins):
        isInBin = (np.abs(x - xbinCenter[ixBin]) <= dxBox/2.)
        xx , yy = x[isInBin] , y[isInBin]
        xMedian[ixBin] , xMean[ixBin] , xStd[ixBin]= np.median(xx) , xx.mean() , xx.std()
        yMedian[ixBin] , yMean[ixBin] , yStd[ixBin]= np.median(yy) , yy.mean() , yy.std()
        nInBin[ixBin] = isInBin.sum()

    return xbinCenter, xMedian, xMean, xStd, yMedian, yMean, yStd, nInBin

def gaussSmooth_YofX(x, y, FWHM):
    '''
    Sloppy function to return the gaussian-smoothed version of an y(x) relation.
    Cid@Lagoa - 07/June/2014
    '''

    sig = FWHM / np.sqrt(8. * np.log(2))
    xS , yS = np.zeros_like(x),  np.zeros_like(x)
    w__ij = np.zeros( (len(x),len(x)) )
    for i in np.arange(len(x)):
        # for j in np.arange(len(x)):
        #     w__ij[i,j] = np.exp( -0.5 * ((x[j] - x[i]) / sig)**2  )

        w__ij[i,:] = np.exp( -0.5 * ((x - x[i]) / sig)**2  )
        w__ij[i,:] = w__ij[i,:] / w__ij[i,:].sum()

        xS[i] = (w__ij[i,:] * x).sum()
        yS[i] = (w__ij[i,:] * y).sum()

    return xS , yS

def calcYofXStats_EqNumberBins(x, y, nPerBin = 25):
    '''
    This gives the statistics of y(x) for x-bins of variable width, but all containing
    the same number of points.
    We 1st sort x, and the y accordingly. Then we compute the median, mean and std
    of x & y in contiguous x-bins in x defined to have nPerBin points each

    Cid@Lagoa - 05/June/2014
    '''

    ind_sx = np.argsort(x)
    xS , yS = x[ind_sx] , y[ind_sx]

    Nbins = len(x) - nPerBin + 1
    xMedian , xMean , xStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    yMedian , yMean , yStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    nInBin = np.zeros(Nbins)

    for ixBin in np.arange(0, Nbins):
        xx , yy = xS[ixBin:ixBin + nPerBin] , yS[ixBin:ixBin + nPerBin]
        xMedian[ixBin] , xMean[ixBin] , xStd[ixBin] = np.median(xx) , xx.mean() , xx.std()
        yMedian[ixBin] , yMean[ixBin] , yStd[ixBin] = np.median(yy) , yy.mean() , yy.std()
        nInBin[ixBin] = len(xx)
    return xMedian, xMean, xStd, yMedian, yMean , yStd, nInBin
    
##############################################################################################
##############################################################################################
##############################################################################################

def nebular_implot(K, tauVNeb__yx, err_tauVNeb__yx, f_obs_HaHb__yx, err_f_obs_HaHb__yx, Lint_Ha__yx, err_Lint_Ha__yx, logZNeb__yx, err_logZNeb__yx, fileName):
    # Setup bins for Radial profiles (in units of HLR)
    RbinIni             = 0.0 
    RbinFin             = 2.0
    RbinStep            = 0.1
    Rbin__r             = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
    RbinCenter__r       = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins              = len(RbinCenter__r)
    
    f, axArr            = plt.subplots(4, 4)
    f.set_size_inches(24, 20)
    
    for ax in f.axes:
        ax.set_axis_off()
        
    ax                  = axArr[0, 0]
    ax.set_axis_on()
    galaxyImgFile       = imgDir + K.califaID + '.jpg'
    galimg              = plt.imread(galaxyImgFile)
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    ax.imshow(galimg)
    
    ax                  = axArr[0, 1]
    ax.set_axis_on()
    tauVNeb__r          = K.radialProfile(tauVNeb__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, tauVNeb__r, 'o-k')
    #ax.tick_params(axis='x', pad=30)
    ax.set_title(r'$\tau_V^{neb}(HLR)$')
    
    ax                  = axArr[0, 2]
    ax.set_axis_on()
    Lint_Ha__r          = K.radialProfile(Lint_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, Lint_Ha__r, 'o-k')
    ax.set_title(r'$L_{H\alpha}^{int}(HLR)$')

    ax                  = axArr[0, 3]
    ax.set_axis_on()
    Lobn__yx            = K.zoneToYX(K.Lobn__z, extensive = True)
    #logZNeb__r          = radialProfileWeighted(logZNeb__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
    logZNeb__r          = K.radialProfile(logZNeb__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, logZNeb__r, 'o-k')
    #ax.set_title(r'$\langle \log\ Z_{neb}\rangle_L (HLR)$')
    ax.set_title(r'$\log\ Z_{neb}(HLR)$')

    ax                  = axArr[1, 1]
    ax.set_axis_on()
    ax.set_title(r'$\tau_V^{neb}$')
    sigma               = tauVNeb__yx.std()
    mean                = tauVNeb__yx.mean()
    vmax                = mean + 2. * sigma
    vmin                = mean - 2. * sigma
    im                  = ax.imshow(tauVNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[2, 1]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (\tau_V^{neb})$')
    sigma               = err_tauVNeb__yx.std()
    mean                = err_tauVNeb__yx.mean()
    vmax                = mean + 2. * sigma
    im                  = ax.imshow(err_tauVNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[3, 1]
    ax.set_axis_on()
    ax.set_title(r'$F^{H\alpha}_{H\beta} / \epsilon(F^{H\alpha}_{H\beta})$')
    signalToNoise       = np.abs(f_obs_HaHb__yx) / np.abs(err_f_obs_HaHb__yx) 
    sigma               = signalToNoise.std()
    mean                = signalToNoise.mean()
    vmax                = mean + 2. * sigma
    im                  = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)

    ax                  = axArr[1, 2]
    ax.set_axis_on()
    ax.set_title(r'$L_{H\alpha}^{int}$')
    sigma               = Lint_Ha__yx.std()
    mean                = Lint_Ha__yx.mean()
    vmax                = mean + sigma
    vmin                = mean - sigma
    im                  = ax.imshow(Lint_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[2, 2]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (L_{H\alpha}^{int})$')
    sigma               = err_Lint_Ha__yx.std()
    mean                = err_Lint_Ha__yx.mean()
    vmax                = mean + sigma
    im                  = ax.imshow(err_Lint_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[3, 2]
    ax.set_axis_on()
    ax.set_title(r'$L_{H\alpha}^{int} / \epsilon(L_{H\alpha}^{int})$')
    signalToNoise       = np.abs(Lint_Ha__yx) / np.abs(err_Lint_Ha__yx) 
    sigma               = signalToNoise.std()
    mean                = signalToNoise.mean()
    vmax                = mean + sigma
    im                  = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)

    ax                  = axArr[1, 3]
    ax.set_axis_on()
    ax.set_title(r'$\log\ Z_{neb}$')
    sigma               = logZNeb__yx.std()
    mean                = logZNeb__yx.mean()
    vmax                = mean + sigma
    vmin                = mean - sigma
    im                  = ax.imshow(logZNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[2, 3]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (log\ Z_{neb})$')
    sigma               = err_logZNeb__yx.std()
    mean                = err_logZNeb__yx.mean()
    vmax                = mean + sigma
    im                  = ax.imshow(err_logZNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[3, 3]
    ax.set_axis_on()
    ax.set_title(r'$\log\ Z_{neb} / \epsilon(log\ Z_{neb})$')
    signalToNoise       = np.abs(logZNeb__yx) / np.abs(err_logZNeb__yx) 
    sigma               = signalToNoise.std()
    mean                = signalToNoise.mean()
    vmax                = mean + sigma
    im                  = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
#    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    f.savefig(fileName)
    plt.close(f)

def plot_3din2d_scatter(x, y, z, xlabel, ylabel, zlabel, fname = None):
    f       = plt.figure(dpi = 96)
    f.set_size_inches(21.88,12.5)
    ax      = f.gca()
    scat    = ax.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)
    
    nBox        = 16
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    ax.plot(xMean, yMean, 'ob-', lw = 2)
    ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'b')
        
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # nPerBin     = len(x) / 350
    # aux         = calcYofXStats_EqNumberBins(x, y, nPerBin)
    # xMedian     = aux[0]
    # xMean       = aux[1]
    # xStd        = aux[2]
    # yMedian     = aux[3]
    # yMean       = aux[4]
    # yStd        = aux[5]
    # nInBin      = aux[6]
    # #ax.plot(xMean, yMean, 'or-', lw = 2)
    # #ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'r')
    # FWHM        = 0.4
    # xS, yS      = gaussSmooth_YofX(xMedian, yMedian, FWHM)
    # ax.plot(xS, yS, 'or--', lw = 0.5)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cb = f.colorbar(ax = ax, mappable = scat)
    cb.set_label(zlabel)
    
    if fname:
        f.savefig(fname)
    else:
        f.show()

def plot_3din2d_scatter_age(x, y, z, xlabel, ylabel, zlabel, age, fname = None):
    f       = plt.figure(dpi = 96)
    f.set_size_inches(21.88,12.5)
    ax      = f.gca()
    scat    = ax.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)
    
    nBox        = 16
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    ax.plot(xMean, yMean, 'ob-', lw = 2)
    ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'b')
        
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # nPerBin     = len(x) / 350
    # aux         = calcYofXStats_EqNumberBins(x, y, nPerBin)
    # xMedian     = aux[0]
    # xMean       = aux[1]
    # xStd        = aux[2]
    # yMedian     = aux[3]
    # yMean       = aux[4]
    # yStd        = aux[5]
    # nInBin      = aux[6]
    # #ax.plot(xMean, yMean, 'or-', lw = 2)
    # #ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'r')
    # FWHM        = 0.4
    # xS, yS      = gaussSmooth_YofX(xMedian, yMedian, FWHM)
    # ax.plot(xS, yS, 'or--', lw = 0.5)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cb = f.colorbar(ax = ax, mappable = scat)
    cb.set_label(zlabel)
    
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    
    if fname:
        f.savefig(fname)
    else:
        f.show()


if __name__ == '__main__':
    #########################################################################        
    ################################## CID ##################################
    ######################################################################### 
    ALL_morfType_GAL__g     = np.ma.zeros((N_gals))
    ALL_at_flux_GAL__g      = np.ma.zeros((N_gals))
    ALL_Mcor_GAL__g         = np.ma.zeros((N_gals))
    ALL_McorSD_GAL__g       = np.ma.zeros((N_gals))
    
    # SK-law *__Trg arrays

    #########################################################################        
    #########################################################################        

    _ALL_tauV__g                = []
    _ALL_tauVNeb__g             = []
    _ALL_tauVNeb_mask__g        = []
    _ALL_logZNeb__g             = []
    _ALL_logZNeb_mask__g        = []
    _ALL_err_tauVNeb__g         = []
    _ALL_Lint_Ha__g             = []
    _ALL_Lint_Ha_mask__g        = []
    _ALL_SFR_Ha__g              = []
    _ALL_SFRSD_Ha__g            = []
    _ALL_tauVNeb2mu__g          = []
    _ALL_McorSD_GAL_zones__g    = []
    _ALL_morfType_GAL_zones__g  = []
    _ALL_WHa__g                 = []
    
    _ALL_tauV__Tg               = []
    _ALL_tauV_mask__Tg          = []
    _ALL_SFR__Tg                = []
    _ALL_SFRSD__Tg              = []

    for iT in range(N_T):
        _ALL_tauV__Tg.append([])
        _ALL_tauV_mask__Tg.append([])
        _ALL_SFR__Tg.append([])
        _ALL_SFRSD__Tg.append([])
        
    for iGal in np.arange(N_gals):
        galName         = listOfPrefixes[iGal][:-1]

        CALIFASuffix    = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        CALIFAFitsFile  = superFitsDir + galName + CALIFASuffix
        emLinesSuffix   = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
        
        if not os.path.isfile(CALIFAFitsFile):
            continue
        
        K = fitsQ3DataCube(CALIFAFitsFile)
        
        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        # read FITSFILE containing galaxy emission lines measured by R.G.B.
        # read_rgb_fits returns False if emLinesFitsFile does not exists.
        read = read_rgb_fits(emLinesFitsFile, read_lines)
        
        if not read:
            ALL_morfType_GAL__g[iGal]           = np.ma.masked
            ALL_at_flux_GAL__g[iGal]            = np.ma.masked
            ALL_Mcor_GAL__g[iGal]               = np.ma.masked
            ALL_McorSD_GAL__g[iGal]             = np.ma.masked
            continue
        
        tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
        ALL_morfType_GAL__g[iGal] = tipo
        aux = np.ones_like(K.Mcor__z) * ALL_morfType_GAL__g[iGal] 
        _ALL_morfType_GAL_zones__g.append(aux)
        
        print '>>> Doing' , iGal , galName , 'hubtyp=', ALL_morfType_GAL__g[iGal], '|  Nzones=' , K.N_zone
        
        #########################################################################        
        ################################## CID ##################################
        ######################################################################### 
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        ALL_McorSD_GAL__g[iGal] = (K.zoneToYX(K.Mcor__z, extensive=True)).mean()
        aux = np.ones_like(K.Mcor__z) * ALL_McorSD_GAL__g[iGal]
        _ALL_McorSD_GAL_zones__g.append(aux)

        # Store galaxy mass & McorSD radial profile. OBS: The McorSD profile is never weighted!
        ALL_Mcor_GAL__g[iGal]   = K.Mcor_tot.sum()

        # Compute & store galaxy-wide at_flux
        numerator__z   = K.Lobn__tZz.sum(axis=1).sum(axis=0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis=1).sum(axis=0)
        ALL_at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()
        
        for iT,tSF in enumerate(tSF__T):
            #--------------------------------------------------------------------------
            # Define mask to pick only populations younger than the input tSF in the computation of SFR & SFRSD.
            flag__t = K.ageBase <= tSF
        
            # SRFSD "raw" image
            # Note that we are NOT dezonifying SFR__z, since it will be compared to the un-dezonifiable tauV!
            aux         = K.Mini__tZz[flag__t,:,:].sum(axis=1).sum(axis=0) / tSF
            SFR__z      = np.ma.masked_array(aux)
            SFRSD__z    = SFR__z / K.zoneArea_pc2
            
            tauV__z     = np.ma.masked_array(K.tau_V__z)
                            
            if mask_xOk:
                # Compute xOk "raw" image
                x__tZz  =  K.popx / K.popx.sum(axis=1).sum(axis=0)
                xOk__z  = x__tZz[flag__t,:,:].sum(axis=1).sum(axis=0)
                
                maskNotOk__z = (xOk__z < xOkMin) | (tauV__z < tauVOkMin)

                tauV__z[maskNotOk__z]   = np.ma.masked
                SFR__z[maskNotOk__z]    = np.ma.masked
                SFRSD__z[maskNotOk__z]  = np.ma.masked
                
            _ALL_tauV__Tg[iT].append(tauV__z.data)
            _ALL_tauV_mask__Tg[iT].append(tauV__z.mask)
            _ALL_SFR__Tg[iT].append(SFR__z)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z)
                
        #########################################################################
        #########################################################################
        #########################################################################
        
        _ALL_tauV__g.append(np.ma.masked_array(K.tau_V__z))

        f_obs__Lz       = read[0]
        err_f_obs__Lz   = read[1]
        fwhm__Lz        = read[2]
        err_fwhm__Lz    = read[3]
        ew__Lz          = read[4]
        lines           = read[5]

        ##################### MASK ######################
        # minimum value of f_lz / err_f_lz
        minSNR = 3.
        
        i_Hb = lines.index('4861')
        i_O3 = lines.index('5007')
        i_Ha = lines.index('6563')
        i_N2 = lines.index('6583')
        
        HbOk = (f_obs__Lz[i_Hb, :] / err_f_obs__Lz[i_Hb, :]) >= minSNR
        O3Ok = (f_obs__Lz[i_O3, :] / err_f_obs__Lz[i_O3, :]) >= minSNR
        HaOk = (f_obs__Lz[i_Ha, :] / err_f_obs__Lz[i_Ha, :]) >= minSNR
        N2Ok = (f_obs__Lz[i_N2, :] / err_f_obs__Lz[i_N2, :]) >= minSNR
        
        maskOk = HbOk & O3Ok & HaOk & N2Ok
        
        f_obs__Lz[:, ~maskOk]       = np.ma.masked
        err_f_obs__Lz[:, ~maskOk]   = np.ma.masked
        ##################################################
        
        _ALL_WHa__g.append(ew__Lz[i_Ha, :])
        
        aux             = calc_Lobs(f_obs__Lz, K.distance_Mpc)
        Lobs__Lz        = aux[0]
        err_Lobs__Lz    = aux[1]
        
        aux                 = calc_tauVNeb(K, f_obs__Lz, err_f_obs__Lz, lines)
        tauVNeb__z          = aux[0] 
        err_tauVNeb__z      = aux[1] 
        f_obs_HaHb__z       = aux[2] 
        err_f_obs_HaHb__z   = aux[3]
        tauVNeb__yx         = K.zoneToYX(tauVNeb__z, extensive = False)
        err_tauVNeb__yx     = K.zoneToYX(err_tauVNeb__z, extensive = False)
        f_obs_HaHb__yx      = K.zoneToYX(f_obs_HaHb__z, extensive = True)
        err_f_obs_HaHb__yx  = K.zoneToYX(err_f_obs_HaHb__z, extensive = True)
        
        _ALL_tauVNeb__g.append(tauVNeb__z.data)
        _ALL_tauVNeb_mask__g.append(tauVNeb__z.mask)
        _ALL_err_tauVNeb__g.append(err_tauVNeb__z)
        
        tauVNeb2mu__yx      = tauVNeb__yx / K.McorSD__yx
                
        aux             = tauV_to_AV(tauVNeb__z, err_tauVNeb__z)
        AVNeb__z        = aux[0]
        err_AVNeb__z    = aux[1]         
        AVNeb__yx       = K.zoneToYX(AVNeb__z, extensive = False)
        err_AVNeb__yx   = K.zoneToYX(err_AVNeb__z, extensive = False)
    
        aux             = calc_Lint_Ha(Lobs__Lz, err_Lobs__Lz, tauVNeb__z, lines)
        Lint_Ha__z      = aux[0]
        err_Lint_Ha__z  = aux[1] 
        Lint_Ha__yx     = K.zoneToYX(Lint_Ha__z, extensive = True)
        err_Lint_Ha__yx = K.zoneToYX(err_Lint_Ha__z, extensive = True)
        
        _ALL_Lint_Ha__g.append(Lint_Ha__z.data)
        _ALL_Lint_Ha_mask__g.append(Lint_Ha__z.mask)
        
        # 3.17 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter        
        SFR_Ha__z       = 3.17 * Lint_Ha__z / (1.e8)
        SFRSD_Ha__z     = SFR_Ha__z / K.zoneArea_pc2
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        
        aux                         = calc_logZNeb(K, f_obs__Lz, err_f_obs__Lz, tauVNeb__z, lines)
        logZNeb__z                  = aux[0]
        err_logZNeb__z              = aux[1]
        O3N2__z                     = aux[2]
        err_O3N2__z                 = aux[3]
        logZNeb__yx                 = K.zoneToYX(logZNeb__z, extensive = False)
        err_logZNeb__yx             = K.zoneToYX(err_logZNeb__z, extensive = False)
        O3N2__yx                    = K.zoneToYX(O3N2__z, extensive = True)
        err_O3N2__yx                = K.zoneToYX(err_O3N2__z, extensive = True)
        
        _ALL_logZNeb__g.append(logZNeb__z.data)
        _ALL_logZNeb_mask__g.append(logZNeb__z.mask)
        
        if plot:
            nebular_implot(K, 
                           tauVNeb__yx, err_tauVNeb__yx, 
                           f_obs_HaHb__yx, err_f_obs_HaHb__yx, 
                           Lint_Ha__yx, err_Lint_Ha__yx,
                           logZNeb__yx, err_logZNeb__yx,  
                           K.califaID + '_' + versionSuffix + '_nebular.png')
    
    aux                         = np.hstack(np.asarray(_ALL_tauVNeb__g))
    auxMask                     = np.hstack(np.asarray(_ALL_tauVNeb_mask__g))
    ALL_tauVNeb__g              = np.ma.masked_array(aux, mask = auxMask)
    ALL_tauV__g                 = np.ma.masked_array(np.hstack(np.asarray(_ALL_tauV__g)), mask = auxMask)
    ALL_err_tauVNeb__g          = np.ma.masked_array(np.hstack(np.asarray(_ALL_err_tauVNeb__g)), mask = auxMask)
    ALL_WHa__g                  = np.ma.masked_array(np.hstack(np.asarray(_ALL_WHa__g)), mask = auxMask)
    
    aux                         = np.hstack(np.asarray(_ALL_logZNeb__g))
    auxMask                     = np.hstack(np.asarray(_ALL_logZNeb_mask__g))
    ALL_logZNeb__g              = np.ma.masked_array(aux, mask = auxMask)

    aux                         = np.hstack(np.asarray(_ALL_Lint_Ha__g))
    auxMask                     = np.hstack(np.asarray(_ALL_Lint_Ha_mask__g))
    ALL_Lint_Ha__g              = np.ma.masked_array(aux, mask = auxMask)
    ALL_SFR_Ha__g               = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR_Ha__g)), mask = auxMask)
    ALL_SFRSD_Ha__g             = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_Ha__g)), mask = auxMask)
    
    ALL_McorSD_GAL_zones__g     = np.ma.masked_array(np.hstack(np.asarray(_ALL_McorSD_GAL_zones__g)))
    ALL_morfType_GAL_zones__g   = np.ma.masked_array(np.hstack(np.asarray(_ALL_morfType_GAL_zones__g)))
    
    ALL_tauV__Tg     = []
    ALL_SFR__Tg      = []
    ALL_SFRSD__Tg    = []
    
    for iT in range(N_T):
        aux     = np.hstack(np.asarray(_ALL_tauV__Tg[iT]))
        auxMask = np.hstack(np.asarray(_ALL_tauV_mask__Tg[iT]))
        ALL_tauV__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_SFR__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR__Tg[iT])), mask = auxMask))
        ALL_SFRSD__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD__Tg[iT])), mask = auxMask))

    SKzero = np.log10(2.5e-4) # -3.60206
    SKslope = 1.4
    Sigma_gas = (np.log10(ALL_SFRSD_Ha__g * 1.e6) - SKzero) / SKslope
    logO_H = ALL_logZNeb__g + np.log10(4.9e-4)
    c = 0
    #c = np.log(0.2)
    logDGR = c + np.log10(ALL_tauVNeb__g) - Sigma_gas

    if plot:
       
        #                    np.log10(ALL_McorSD_GAL_zones__g[m]),
        #                    r'$\log\ \mu^{galaxy}_\star\ [M_\odot pc^{-2}]$'
        #                    12 + logO_H[m],
        #                    r'$12 + \log\ ($O/H$)$',
        #                    np.log10(ALL_McorSD_GAL_zones__g[m]),
        #                    r'$\log\ \mu^{galaxy}_\star\ [M_\odot pc^{-2}]$',

        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0) & (ALL_morfType_GAL_zones__g > 8.5)
        #m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0)
        plot_3din2d_scatter(ALL_logZNeb__g[m], logDGR[m], ALL_morfType_GAL_zones__g[m],
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            'logZNeb_logDGR_morphType.png')

        m = (ALL_err_tauVNeb__g<0.15) & (np.log10(ALL_tauVNeb__g) > -1) & (np.log10(ALL_tauV__g) > -1)
        #& (ALL_morfType_GAL_zones__g > 8.5)
        #m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0) & (ALL_tauV__g > 0)
        plot_3din2d_scatter(np.log10(ALL_tauVNeb__g[m]), np.log10(ALL_SFRSD_Ha__g[m]), 10. ** ALL_logZNeb__g[m],
                            r'$\log\ \tau_V^{neb}$', r'$\log\ \Sigma_{SFR}\ [M_\odot yr^{-1} pc^{-2}]$', r'$Z_{neb}\ [Z_\odot]$',
                            'logtauVNeb_logSFRSD_ZNeb.png')

        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g >= 0)
        plot_3din2d_scatter(ALL_logZNeb__g[m], logDGR[m], np.log10(ALL_WHa__g[m]),
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            'logZNeb_logDGR_logWHa.png')
        
        m &= (ALL_tauV__g >= 0)
        plot_3din2d_scatter(ALL_tauVNeb__g[m], ALL_tauV__g[m], 10. ** ALL_logZNeb__g[m],
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                            'tauV_tauVNeb_ZNeb.png')
        plot_3din2d_scatter(ALL_tauV__g[m], ALL_tauVNeb__g[m], np.log10(ALL_SFRSD_Ha__g[m]), 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}\ [M_\odot yr^{-1} pc^{-2}]$',
                            'tauV_tauVNeb_logSFRSD.png')
        plot_3din2d_scatter(ALL_tauV__g[m], ALL_tauVNeb__g[m], np.log10(ALL_WHa__g[m]), 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            'tauV_tauVNeb_logWHa.png')
        
        m &= (ALL_morfType_GAL_zones__g > 8.5)
        plot_3din2d_scatter(ALL_tauVNeb__g[m], ALL_tauV__g[m], ALL_morfType_GAL_zones__g[m],
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'Morphology Type',
                            'tauV_tauVNeb_morphType.png')
        
        for iT, tSF in enumerate(tSF__T):
            print iT, tSF
            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g >= tauVOkMin)
            x = ALL_SFR_Ha__g
            y = np.log10(ALL_SFR__Tg[iT])
            z = 10. ** ALL_logZNeb__g
            plot_3din2d_scatter_age(np.log10(x[m]), y[m], z[m],
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, 'logSFRHa_logSFR_ZNeb_age_%d.png' % iT)

            m &= (ALL_morfType_GAL_zones__g > 8.5)
            z = ALL_morfType_GAL_zones__g
            plot_3din2d_scatter_age(np.log10(x[m]), y[m], z[m],
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'Morphology Type',
                                    tSF, 'logSFRHa_logSFR_morphType_age_%d.png' % iT)
            
            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g >= 0) & (ALL_tauV__Tg[iT] >= 0)
            plot_3din2d_scatter_age(ALL_tauV__Tg[iT][m], ALL_tauVNeb__g[m], 10. ** ALL_logZNeb__g[m],
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, 'tauV_tauVNeb_ZNeb_age_%d.png' % iT)
            plot_3din2d_scatter_age(ALL_tauV__Tg[iT][m], ALL_tauVNeb__g[m], np.log10(ALL_SFRSD_Ha__g[m]), 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, 'tauV_tauVNeb_logSFRSD_age_%d.png' % iT)
            plot_3din2d_scatter_age(ALL_tauV__Tg[iT][m], ALL_tauVNeb__g[m], np.log10(ALL_WHa__g[m]), 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                                    tSF, 'tauV_tauVNeb_logWHa_age_%d.png' % iT)


##############################################################################################
##############################################################################################
##############################################################################################
