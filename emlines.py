#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import numpy as np
from pycasso import fitsQ3DataCube
import matplotlib as mpl
from matplotlib import pyplot as plt
from get_morfologia import get_morfologia
from scipy import stats as st
import os
from lines import *
import errno

nebular_galaxy_plot = True
#nebular_galaxy_plot = False
plot = False
#plot = True
debug = False
#debug = True
save_npz = False
#save_npz = True

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size']       = 30
mpl.rcParams['axes.labelsize']  = 30
mpl.rcParams['axes.titlesize']  = 32
mpl.rcParams['xtick.labelsize'] = 26
mpl.rcParams['ytick.labelsize'] = 26 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

tauFilteredDir = 'tauVMasked'

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
    
galaxiesListFile    = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
#galaxiesListFile    = CALIFAWorkDir + 'listAll.txt'
baseCode            = 'Bgsd6e'
#versionSuffix       = 'px1_q043.d14a'
versionSuffix       = 'v20_q043.d14a'
#superFitsDir        = '/Volumes/backupzeira/CALIFA/q043/px1/'
superFitsDir        = CALIFAWorkDir + 'gal_fits/' + versionSuffix + '/'

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

f               = open(galaxiesListFile, 'r')
listOfPrefixes  = f.readlines()
f.close()

RbinIni , RbinFin , RbinStep = 0.0 , 2.0 , 0.1
Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins = len(RbinCenter__r)

if debug:
    listOfPrefixes = listOfPrefixes[0:20]        # Q&D tests ...
    #listOfPrefixes = ['K0026\n']
    
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
def create_dir(dir):
    try:
        os.makedirs(dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print "\nBE CAREFUL! Directory %s already exists." % dir

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

def calc_Lint_Ha(L_obs__Lz, L_obs_err__Lz, tau_V_neb__z,lines):
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    
    eHa = np.ma.exp(qCCM['6563'] * tau_V_neb__z)
    LobsHaHb = L_obs__Lz[i_Ha, :] / L_obs__Lz[i_Hb, :]

    Lint_Ha__z = L_obs__Lz[i_Ha, :] * eHa
    
    a = L_obs_err__Lz[i_Ha, :]
    b = q * LobsHaHb * L_obs_err__Lz[i_Hb, :]
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

def nebular_implot(K, tau_V_neb__yx, err_tauVNeb__yx, f_obs_HaHb__yx, err_f_obs_HaHb__yx, Lint_Ha__yx, err_Lint_Ha__yx, logZ_neb__yx, logZ_neb_err__yx, fileName):
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
    tau_V_neb__r          = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, tau_V_neb__r, 'o-k')
    #ax.tick_params(axis='x', pad=30)
    ax.set_title(r'$\tau_V^{neb}(R)$')
    
    ax                  = axArr[0, 2]
    ax.set_axis_on()
    Lint_Ha__r          = K.radialProfile(Lint_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, Lint_Ha__r, 'o-k')
    ax.set_title(r'$L_{H\alpha}^{int}(R)$')

    ax                  = axArr[0, 3]
    ax.set_axis_on()
    #Lobn__yx            = K.zoneToYX(K.Lobn__z, extensive = True)
    #logZNeb__r          = radialProfileWeighted(logZ_neb__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
    logZNeb__r          = K.radialProfile(logZ_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, logZNeb__r, 'o-k')
    #ax.set_title(r'$\langle \log\ Z_{neb}\rangle_L (HLR)$')
    ax.set_title(r'$\log\ Z_{neb}(R)$')

    ax                  = axArr[1, 1]
    ax.set_axis_on()
    ax.set_title(r'$\tau_V^{neb}$')
    sigma               = tau_V_neb__yx.std()
    mean                = tau_V_neb__yx.mean()
    vmax                = mean + 2. * sigma
    vmin                = mean - 2. * sigma
    im                  = ax.imshow(tau_V_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
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
    sigma               = logZ_neb__yx.std()
    mean                = logZ_neb__yx.mean()
    vmax                = mean + sigma
    vmin                = mean - sigma
    im                  = ax.imshow(logZ_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[2, 3]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (log\ Z_{neb})$')
    sigma               = logZ_neb_err__yx.std()
    mean                = logZ_neb_err__yx.mean()
    vmax                = mean + sigma
    im                  = ax.imshow(logZ_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax                  = axArr[3, 3]
    ax.set_axis_on()
    ax.set_title(r'$\log\ Z_{neb} / \epsilon(log\ Z_{neb})$')
    signalToNoise       = np.abs(logZ_neb__yx) / np.abs(logZ_neb_err__yx) 
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
        
    plt.close(f)

def plot_3din2d_scatter_age(x, y, z, xlabel, ylabel, zlabel, age, fname = None):
    f       = plt.figure(dpi = 96)
    f.set_size_inches(21.88,12.5)
    ax      = f.gca()
    scat    = ax.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = 's: %.2f' % rhoSpearman

    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)

    ax.text(0.90, 0.90, txt,
            fontsize = 28,
            transform = ax.transAxes,
            verticalalignment = 'top',
            bbox = textbox)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # 
    # nBox        = 16
    # dxBox       = (x.max() - x.min()) / (nBox - 1.)
    # aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
    # xbinCenter  = aux[0]
    # xMedian     = aux[1]
    # xMean       = aux[2]
    # xStd        = aux[3]
    # yMedian     = aux[4]
    # yMean       = aux[5]
    # yStd        = aux[6]
    # nInBin      = aux[7]
    # ax.plot(xMean, yMean, 'ob-', lw = 2)
    # ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'b')
    #     
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
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

    plt.close(f)

#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
def calc_alogZ_Stuff(K, flag__t, Rbin__r, weiRadProf = False, xOkMin = 0.10):
    '''
    Compute average logZ (alogZ_*) both for each zone (*__z) and the galaxy-wide
    average (*_GAL, computed a la GD14).

    Only st-pops satisfying the input flag__t (ageBase-related) mask are considered!
    This allows us to compute alogZ_* for, say, < 2 Gyr or 1--7 Gyr populations,
    as well as the whole-age range (using flag__t = True for all base ages)
    with a same function and saving the trouble of keeping separate variables for the same thing:-)

    Radial profiles of alogZ_mass and alogZ_flux are also computed. They are/are-not weighted
    (by Mcor__yx & Lobn__yx, respectively) according to weiRadProf.

    ==> return alogZ_mass_GAL, alogZ_flux_GAL, isOkFrac_GAL , alogZ_mass__r, alogZ_flux__r

    Cid@Lagoa - 05/Jun/2014


!!HELP!! ATT: My way of computing alogZ_*__z produces nan', which are ugly but harmless.
    I tried to fix it using masked arrays:

    alogZ_mass__z  = np.ma.masked_array( numerator__z/(denominator__z+0e-30) , mask = (denominator__z == 0))

    but this did not work!

    Cid@Lagoa - 20/Jun/2014
    '''

    #--------------------------------------------------------------------------
    # Initialization

    # Define log of base metallicities **in solar units** for convenience
    logZBase__Z = np.log10(K.metBase / Zsun)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Define alogZ_****__z: flux & mass weighted average logZ for each zone

    # ==> alogZ_mass__z - ATT: There may be nan's here depending on flag__t!
    numerator__z = np.tensordot(K.Mcor__tZz[flag__t, :, :] , logZBase__Z , (1, 0)).sum(axis = 0)
    denominator__z = K.Mcor__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    alogZ_mass__z = numerator__z / denominator__z

    # ==> alogZ_flux__z - ATT: There may be nan's here depending on flag__t!
    numerator__z = np.tensordot(K.Lobn__tZz[flag__t, :, :] , logZBase__Z , (1, 0)).sum(axis = 0)
    denominator__z = K.Lobn__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    alogZ_flux__z = numerator__z / denominator__z
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Def galaxy-wide averages of alogZ in light & mass, but **discards** zones having
    # too little light fractions in the ages given by flag__t

    # Define Ok flag: Zones with light fraction x < xOkMin are not reliable for alogZ (& etc) estimation!
    x__tZz = K.popx / K.popx.sum(axis = 1).sum(axis = 0)
    xOk__z = x__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    isOk__z = xOk__z > xOkMin

    # Fraction of all zones which are Ok in the isOk__z sense. Useful to censor galaxies whose
    # galaxy-wide averages are based on too few zones (hence not representative)
    # OBS: isOkFrac_GAL is not used in this function, but it's returned to be used by the caller
    isOkFrac_GAL = (1.0 * isOk__z.sum()) / (1.0 * K.N_zone)

    # Galaxy wide averages of logZ - ATT: Only isOk__z zones are considered in this averaging!
    numerator__z = K.Mcor__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0) * alogZ_mass__z
    denominator__z = K.Mcor__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    alogZ_mass_GAL = numerator__z[isOk__z].sum() / denominator__z[isOk__z].sum()

    numerator__z = K.Lobn__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0) * alogZ_flux__z
    denominator__z = K.Lobn__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    alogZ_flux_GAL = numerator__z[isOk__z].sum() / denominator__z[isOk__z].sum()
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Radial profiles of alogZ_mass & alogZ_flux. The decision to use proper weighted averaging
    # with the bins (by Mcor__yx for alogZ_mass__r and by Lobn__yx for alogZ_flux__r) is taken
    # cf the input weiRadProf boolean switch.

    # Define aY_****__yx images - needed to compute radial profiles
    alogZ_mass__yx = K.zoneToYX(alogZ_mass__z, extensive = False)
    alogZ_flux__yx = K.zoneToYX(alogZ_flux__z, extensive = False)

    # Compute radial profiles weighted or not.
    # OBS: Our defs of Mcor__yx & Lobn__yx below return SD's (ie, Msun/pc2 and Lsun/A/pc2).
    #      This is NOT what I originally intended, but it works for the current weighting purposes
    if weiRadProf:
        Mcor__yx = K.zoneToYX(K.Mcor__z, extensive = True)
        Lobn__yx = K.zoneToYX(K.Lobn__z, extensive = True)
        alogZ_mass__r = radialProfileWeighted(alogZ_mass__yx, Mcor__yx, Rbin__r, K.HLR_pix, K.radialProfile)
        alogZ_flux__r = radialProfileWeighted(alogZ_flux__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
    else:
        alogZ_mass__r = K.radialProfile(alogZ_mass__yx, Rbin__r, rad_scale = K.HLR_pix)
        alogZ_flux__r = K.radialProfile(alogZ_flux__yx, Rbin__r, rad_scale = K.HLR_pix)
    #--------------------------------------------------------------------------

    return alogZ_mass_GAL, alogZ_flux_GAL, isOkFrac_GAL , alogZ_mass__r, alogZ_flux__r
#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

def plotSFR(x,y,xlabel,ylabel,xlim,ylim,age,fname):
    f = plt.figure(dpi = 96)
    f.set_size_inches(21.88,12.5)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    if fname:
        f.savefig(fname)
    else:
        f.show()
    plt.close(f)

def plotTau(x,y,xlabel,ylabel,xlim,ylim,age,fname):
    f = plt.figure(dpi = 96)
    f.set_size_inches(21.88,12.5)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    #ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    if fname:
        f.savefig(fname)
    else:
        f.show()
    plt.close(f)

if __name__ == '__main__':
    create_dir(tauFilteredDir)
    
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
    _ALL_tau_V_neb__g             = []
    _ALL_tau_V_neb_mask__g        = []
    _ALL_logZNeb__g             = []
    _ALL_logZNeb_mask__g        = []
    _ALL_err_tauVNeb__g         = []
    _ALL_Lint_Ha__g             = []
    _ALL_Lint_Ha_mask__g        = []
    _ALL_SFR_Ha__g              = []
    _ALL_SFRSD_Ha__g            = []
    _ALL_tauVNeb2mu__g          = []
    _ALL_Mcor__g      = []
    _ALL_McorSD_GAL_zones__g    = []
    _ALL_morfType_GAL_zones__g  = []
    _ALL_WHa__g                 = []
    
    _ALL_tauV__Tg               = []
    _ALL_tauV_mask__Tg          = []
    _ALL_SFR__Tg                = []
    _ALL_SFRSD__Tg              = []

    ALL_aSFRSD_Ha__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_alogSFRSD_Ha__rg        = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD__Trg             = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogSFRSD__Trg          = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_tau_V_neb__rg             = np.ma.zeros((NRbins, N_gals))
    ALL_tau_V__Trg               = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_morfType_GAL_zones__rg  = np.ma.zeros((NRbins, N_gals))
 
    ALL_alogZ_mass_GAL__Tg      = np.zeros((N_T, N_gals))
    ALL_alogZ_flux_GAL__Tg      = np.zeros((N_T, N_gals))
    ALL_isOkFrac_GAL__Tg        = np.zeros((N_T, N_gals))
    ALL_alogZ_mass__Trg         = np.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_flux__Trg         = np.zeros((N_T, NRbins, N_gals))
    
    N_zones = 0

    #metBase = np.array([  9.00000000e-05,   3.70000000e-04,   9.30000000e-04,  3.70000000e-03,   7.56000000e-03,   1.90000000e-02, 3.15300000e-02])                  
    #ALL_Rs_histo__Z             = np.ma.masked_all(metBase.shape) 

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
        
        if not (os.path.isfile(CALIFAFitsFile) and os.path.isfile(emLinesFitsFile)):
            ALL_morfType_GAL__g[iGal]           = np.ma.masked
            ALL_at_flux_GAL__g[iGal]            = np.ma.masked
            ALL_Mcor_GAL__g[iGal]               = np.ma.masked
            ALL_McorSD_GAL__g[iGal]             = np.ma.masked
            continue
        
        K = fitsQ3DataCube(CALIFAFitsFile)
        
        N_zones += K.N_zone
        
        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        # read FITSFILE containing galaxy emission lines measured by R.G.B.
        # read_rgb_fits returns False if emLinesFitsFile does not exists.
        #read = read_rgb_fits(emLinesFitsFile, read_lines)
        K.loadEmLinesDataCube(emLinesFitsFile)
        
        tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
        ALL_morfType_GAL__g[iGal] = tipo
        aux = np.ones_like(K.Mcor__z) * ALL_morfType_GAL__g[iGal] 
        _ALL_morfType_GAL_zones__g.append(aux)
        ALL_morfType_GAL_zones__rg[:, iGal] = np.ones_like(RbinCenter__r) * ALL_morfType_GAL__g[iGal]
        
        print '>>> Doing' , iGal , galName , 'hubtyp=', ALL_morfType_GAL__g[iGal], '|  Nzones=' , K.N_zone
        
        #########################################################################        
        ################################## CID ##################################
        ######################################################################### 
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        ALL_McorSD_GAL__g[iGal] = K.McorSD__yx.mean()
        aux = np.ones_like(K.Mcor__z) * ALL_McorSD_GAL__g[iGal]
        _ALL_McorSD_GAL_zones__g.append(aux)
        _ALL_Mcor__g.append(K.Mcor__z)
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

            weiRadProf = True
                            
            aux = calc_alogZ_Stuff(K, flag__t, Rbin__r, weiRadProf, xOkMin = xOkMin)
            ALL_alogZ_mass_GAL__Tg[iT, iGal] = aux[0]
            ALL_alogZ_flux_GAL__Tg[iT, iGal] = aux[1]
            ALL_isOkFrac_GAL__Tg[iT, iGal] = aux[2]
            ALL_alogZ_mass__Trg[iT, :, iGal] = aux[3]
            ALL_alogZ_flux__Trg[iT, :, iGal] = aux[4]

            SFRSD__yx = K.zoneToYX(SFR__z, extensive = True)
            aSFRSD__r = K.radialProfile(SFRSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            alogSFRSD__r = K.radialProfile(np.log10(SFRSD__yx), Rbin__r, rad_scale = K.HLR_pix)
            ALL_aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
            ALL_alogSFRSD__Trg[iT, :, iGal] = alogSFRSD__r
            
            _ALL_tauV__Tg[iT].append(tauV__z.data)
            _ALL_tauV_mask__Tg[iT].append(tauV__z.mask)
            
            tauV__yx = K.zoneToYX(tauV__z, extensive=False)
            ALL_tau_V__Trg[iT, :, iGal] = K.radialProfile(tauV__yx, Rbin__r, rad_scale = K.HLR_pix)
            
            _ALL_SFR__Tg[iT].append(SFR__z)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z)
                
        #########################################################################
        #########################################################################
        #########################################################################
        
        _ALL_tauV__g.append(np.ma.masked_array(K.tau_V__z))

        f_obs__Lz       = K.EL.flux
        err_f_obs__Lz   = K.EL.eflux
        fwhm__Lz        = K.EL.fwhm
        err_fwhm__Lz    = K.EL.efwhm
        ew__Lz          = K.EL.EW
        lines           = K.EL.lines

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
        #N2HaOk = (np.log10(f_obs__Lz[i_N2, :] / f_obs__Lz[i_Ha, :]) <= -0.4)
        
        N2Ha = np.log10(K.EL.N2_obs__z/K.EL.Ha_obs__z)
        O3Hb = np.log10(K.EL.O3_obs__z/K.EL.Hb_obs__z)
        
        maskOk = HbOk & O3Ok & HaOk & N2Ok 
        #& N2HaOk
        L = Lines()
        maskOk = HbOk & O3Ok & HaOk & N2Ok & L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # f = plt.figure()
        # ax = f.gca() 
        # ax.scatter(N2Ha, O3Hb, c = 'r')
        # ax.scatter(N2Ha[m], O3Hb[m], c = 'b')
        # ax.scatter(N2Ha[maskOk], O3Hb[maskOk], c = 'y')
        # ax.plot(L.x['S06'], L.y['S06'], label = 'S06')
        # plt.axis([-1.0, 1.0, -1.0, 1.0])
        # f.savefig('%s-testeBPT.png' % K.califaID)
        # plt.close(f)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        
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
        tau_V_neb__yx         = K.zoneToYX(tauVNeb__z, extensive = False)
        err_tauVNeb__yx     = K.zoneToYX(err_tauVNeb__z, extensive = False)
        f_obs_HaHb__yx      = K.zoneToYX(f_obs_HaHb__z, extensive = True)
        err_f_obs_HaHb__yx  = K.zoneToYX(err_f_obs_HaHb__z, extensive = True)
        
        _ALL_tau_V_neb__g.append(tauVNeb__z.data)
        _ALL_tau_V_neb_mask__g.append(tauVNeb__z.mask)
        _ALL_err_tauVNeb__g.append(err_tauVNeb__z)
        
        maskNotOkTauVNeb__z = tauVNeb__z <= 0
        #maskNotOkTauVNeb__z = tauVNeb__z < 0.05
        tauVNeb__z[maskNotOkTauVNeb__z] = np.ma.masked
        tauVNeb_masked__yx       = K.zoneToYX(tauVNeb__z, extensive = False)
        ALL_tau_V_neb__rg[:, iGal] = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix) 
        
        tauVNeb2mu__yx      = tau_V_neb__yx / K.McorSD__yx
                
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
        
        SFRSD_Ha__yx    = K.zoneToYX(SFR_Ha__z, extensive = True)
        aSFRSD_Ha__r    = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        alogSFRSD_Ha__r = K.radialProfile(np.log10(SFRSD_Ha__yx), Rbin__r, rad_scale = K.HLR_pix)
        ALL_aSFRSD_Ha__rg[:, iGal]    = aSFRSD_Ha__r
        ALL_alogSFRSD_Ha__rg[:, iGal] = alogSFRSD_Ha__r
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        
        aux                         = calc_logZNeb(K, f_obs__Lz, err_f_obs__Lz, tauVNeb__z, lines)
        logZNeb__z                  = aux[0]
        err_logZNeb__z              = aux[1]
        O3N2__z                     = aux[2]
        err_O3N2__z                 = aux[3]
        logZ_neb__yx                 = K.zoneToYX(logZNeb__z, extensive = False)
        logZ_neb_err__yx             = K.zoneToYX(err_logZNeb__z, extensive = False)
        O3N2__yx                    = K.zoneToYX(O3N2__z, extensive = True)
        err_O3N2__yx                = K.zoneToYX(err_O3N2__z, extensive = True)
        
        _ALL_logZNeb__g.append(logZNeb__z.data)
        _ALL_logZNeb_mask__g.append(logZNeb__z.mask)
        
        if nebular_galaxy_plot:
            nebular_implot(K, 
                           tau_V_neb__yx, err_tauVNeb__yx, 
                           f_obs_HaHb__yx, err_f_obs_HaHb__yx, 
                           Lint_Ha__yx, err_Lint_Ha__yx,
                           logZ_neb__yx, logZ_neb_err__yx,  
                           K.califaID + '_' + versionSuffix + '_nebular.png')
    
    aux                         = np.hstack(np.asarray(_ALL_tau_V_neb__g))
    auxMask                     = np.hstack(np.asarray(_ALL_tau_V_neb_mask__g))
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
    
    ALL_Mcor__g       = np.ma.masked_array(np.hstack(np.asarray(_ALL_Mcor__g)))
    ALL_McorSD_GAL_zones__g     = np.ma.masked_array(np.hstack(np.asarray(_ALL_McorSD_GAL_zones__g)))
    ALL_morfType_GAL_zones__g   = np.ma.masked_array(np.hstack(np.asarray(_ALL_morfType_GAL_zones__g)))
    
    ALL_tau_V__Tg     = []
    ALL_SFR__Tg      = []
    ALL_SFRSD__Tg    = []
    
    for iT in range(N_T):
        aux     = np.hstack(np.asarray(_ALL_tauV__Tg[iT]))
        auxMask = np.hstack(np.asarray(_ALL_tauV_mask__Tg[iT]))
        ALL_tau_V__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_SFR__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR__Tg[iT])), mask = auxMask))
        ALL_SFRSD__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD__Tg[iT])), mask = auxMask))\

    SKzero = np.log10(1.6e-4)
    SKslope = 1.4
    Sigma_gas = (np.log10(ALL_SFRSD_Ha__g * 1.e6) - SKzero) / SKslope
    logO_H = ALL_logZNeb__g + np.log10(4.9e-4)
    c = 0
    #c = np.log(0.2)
    logDGR = c + np.log10(ALL_tauVNeb__g) - Sigma_gas
    
    print 'NZones: %d' % N_zones
    print ALL_McorSD_GAL_zones__g.shape
    
    if save_npz:
        D = {        
             'ALL_morfType_GAL__g'         : ALL_morfType_GAL__g, 
             'ALL_at_flux_GAL__g'          : ALL_at_flux_GAL__g,
             'ALL_Mcor_GAL__g'             : ALL_Mcor_GAL__g, 
             'ALL_McorSD_GAL__g'           : ALL_McorSD_GAL__g,
             'ALL_tauVNeb__g'              : ALL_tauVNeb__g,
             'ALL_tauV__g'                 : ALL_tauV__g,
             'ALL_err_tauVNeb__g'          : ALL_err_tauVNeb__g,
             'ALL_WHa__g'                  : ALL_WHa__g,
             'ALL_logZNeb__g'              : ALL_logZNeb__g,
             'ALL_Lint_Ha__g'              : ALL_Lint_Ha__g,
             'ALL_SFR_Ha__g'               : ALL_SFR_Ha__g,
             'ALL_SFRSD_Ha__g'             : ALL_SFRSD_Ha__g,
             'ALL_Mcor__g'       : ALL_Mcor__g,
             'ALL_McorSD_GAL_zones__g'     : ALL_McorSD_GAL_zones__g,
             'ALL_morfType_GAL_zones__g'   : ALL_morfType_GAL_zones__g,
             'ALL_morfType_GAL_zones__rg'  : ALL_morfType_GAL_zones__rg,
             'ALL_tau_V__Tg'                : ALL_tau_V__Tg,
             'ALL_SFR__Tg'                 : ALL_SFR__Tg,
             'ALL_SFRSD__Tg'               : ALL_SFRSD__Tg,
             'ALL_aSFRSD_Ha__rg'           : ALL_aSFRSD_Ha__rg,
             'ALL_alogSFRSD_Ha__rg'        : ALL_alogSFRSD_Ha__rg,
             'ALL_aSFRSD__Trg'             : ALL_aSFRSD__Trg,
             'ALL_alogSFRSD__Trg'          : ALL_alogSFRSD__Trg,
             'ALL_alogZ_mass_GAL__Tg'      : ALL_alogZ_mass_GAL__Tg,
             'ALL_alogZ_flux_GAL__Tg'      : ALL_alogZ_flux_GAL__Tg,
             'ALL_isOkFrac_GAL__Tg'        : ALL_isOkFrac_GAL__Tg,
             'ALL_alogZ_mass__Trg'         : ALL_alogZ_mass__Trg,
             'ALL_alogZ_flux__Trg'         : ALL_alogZ_flux__Trg,
             'ALL_tau_V__Trg'               : ALL_tau_V__Trg,
             'ALL_tau_V_neb__rg'             : ALL_tau_V_neb__rg,
        }
        np.savez_compressed("emlines.npz", D)
        
    if plot:
        ###################
        # tauV not masked #
        ###################
        m = (ALL_morfType_GAL_zones__g > 8.5)
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            'logZNeb_logDGR_morphType.png')
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            'logZNeb_logDGR_logWHa.png')
        x = np.log10(ALL_Mcor__g[m])
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ M_\star\ [M_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            'logMcor_logDGR_morphType.png')
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = 10. ** ALL_logZNeb__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                            'tauV_tauVNeb_ZNeb.png')
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_SFRSD_Ha__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                            'tauV_tauVNeb_logSFRSDNeb.png')
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            'tauV_tauVNeb_logWHa.png')
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'Morphology Type',
                            'tauV_tauVNeb_morphType.png')

        ##################
        # tauVNeb masked #
        ##################
        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0) & (ALL_morfType_GAL_zones__g > 8.5)
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            '%s/logZNeb_logDGR_morphType.png' % tauFilteredDir)
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            '%s/logZNeb_logDGR_logWHa.png' % tauFilteredDir)
        x = np.log10(ALL_Mcor__g[m])
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ M_\star\ [M_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            '%s/logMcor_logDGR_morphType.png' % tauFilteredDir)

        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = 10. ** ALL_logZNeb__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                            '%s/tauV_tauVNeb_ZNeb.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_SFRSD_Ha__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                            '%s/tauV_tauVNeb_logSFRSDNeb.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            '%s/tauV_tauVNeb_logWHa.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'Morphology Type',
                            '%s/tauV_tauVNeb_morphType.png' % tauFilteredDir)

        #######################################################################################
        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0.05) & (ALL_morfType_GAL_zones__g > 8.5)
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            '%s/logZNeb_logDGR_morphType2.png' % tauFilteredDir)
        x = ALL_logZNeb__g[m]
        y = logDGR[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ Z_{neb}\ [Z_\odot]$', r'$\log\ \delta_{DGR}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            '%s/logZNeb_logDGR_logWHa2.png' % tauFilteredDir)
        x = np.log10(ALL_Mcor__g[m])
        y = logDGR[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\log\ M_\star\ [M_\odot]$', r'$\log\ \delta_{DGR}$', r'Morphology type',
                            '%s/logMcor_logDGR_morphType2.png' % tauFilteredDir)

        m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0.05)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = 10. ** ALL_logZNeb__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                            '%s/tauV_tauVNeb_ZNeb2.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_SFRSD_Ha__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                            '%s/tauV_tauVNeb_logSFRSDNeb2.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = np.log10(ALL_WHa__g[m])
        plot_3din2d_scatter(x, y, z, 
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                            '%s/tauV_tauVNeb_logWHa2.png' % tauFilteredDir)
        x = ALL_tauV__g[m]
        y = ALL_tauVNeb__g[m]
        z = ALL_morfType_GAL_zones__g[m]
        plot_3din2d_scatter(x, y, z,
                            r'$\tau_V$', r'$\tau_V^{neb}$', r'Morphology Type',
                            '%s/tauV_tauVNeb_morphType2.png' % tauFilteredDir)

        for iT, tSF in enumerate(tSF__T):
            ######################
            # tauVNeb not masked #
            ######################
            x = np.log10(ALL_SFR_Ha__g)
            y = np.log10(ALL_SFR__Tg[iT])
            z = 10. ** ALL_logZNeb__g
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, 'logSFRNeb_logSFR_ZNeb_age_%sMyr.png' % str(tSF / 1.e6))
            m = (ALL_morfType_GAL_zones__g > 8.5)
            x = np.log10(ALL_SFR_Ha__g[m])
            y = np.log10(ALL_SFR__Tg[iT][m])
            z = ALL_morfType_GAL_zones__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'Morphology Type',
                                    tSF, 'logSFRNeb_logSFR_morphType_age_%sMyr.png' % str(tSF / 1.e6))
            x = ALL_tau_V__Tg[iT]
            y = ALL_tauVNeb__g
            z = 10. ** ALL_logZNeb__g
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, 'tauV_tauVNeb_ZNeb_age_%sMyr.png' % str(tSF / 1.e6))
            x = ALL_tau_V__Tg[iT]
            y = ALL_tauVNeb__g
            z = np.log10(ALL_SFRSD_Ha__g)
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, 'tauV_tauVNeb_logSFRSDNeb_age_%sMyr.png' % str(tSF / 1.e6))
            x = ALL_tau_V__Tg[iT]
            y = ALL_tauVNeb__g
            z = np.log10(ALL_SFRSD__Tg[iT])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^\star\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, 'tauV_tauVNeb_logSFRSD_age_%sMyr.png' % str(tSF / 1.e6))
            x = ALL_tau_V__Tg[iT]
            y = ALL_tauVNeb__g
            z = np.log10(ALL_WHa__g)
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                                    tSF, 'tauV_tauVNeb_logWHa_age_%sMyr.png' % str(tSF / 1.e6))

            ##################
            # tauVNeb masked #
            ##################
            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0)
            x = np.log10(ALL_SFR_Ha__g[m])
            y = np.log10(ALL_SFR__Tg[iT][m])
            z = 10. ** ALL_logZNeb__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, '%s/logSFRHa_logSFR_ZNeb_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))
            m &= (ALL_morfType_GAL_zones__g > 8.5)
            x = np.log10(ALL_SFR_Ha__g[m])
            y = np.log10(ALL_SFR__Tg[iT][m])
            z = ALL_morfType_GAL_zones__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'Morphology Type',
                                    tSF, '%s/logSFRHa_logSFR_morphType_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))

            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0)
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = 10. ** ALL_logZNeb__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, '%s/tauV_tauVNeb_ZNeb_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_SFRSD_Ha__g[m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, '%s/tauV_tauVNeb_logSFRSDNeb_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_SFRSD__Tg[iT][m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^\star\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, '%s/tauV_tauVNeb_logSFRSD_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_WHa__g[m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                                    tSF, '%s/tauV_tauVNeb_logWHa_age_%sMyr.png' % (tauFilteredDir, str(tSF / 1.e6)))

            #######################################################################################
            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0.05)
            x = np.log10(ALL_SFR_Ha__g[m])
            y = np.log10(ALL_SFR__Tg[iT][m])
            z = 10. ** ALL_logZNeb__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, '%s/logSFRHa_logSFR_ZNeb_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))
            m &= (ALL_morfType_GAL_zones__g > 8.5)
            x = np.log10(ALL_SFR_Ha__g[m])
            y = np.log10(ALL_SFR__Tg[iT][m])
            z = ALL_morfType_GAL_zones__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\log\ SFR_{neb}$', r'$\log\ SFR_\star$', r'Morphology Type',
                                    tSF, '%s/logSFRHa_logSFR_morphType_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))

            m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0.05)
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = 10. ** ALL_logZNeb__g[m]
            plot_3din2d_scatter_age(x, y, z,
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$Z_{neb}\ [Z_\odot]$',
                                    tSF, '%s/tauV_tauVNeb_ZNeb_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_SFRSD_Ha__g[m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, '%s/tauV_tauVNeb_logSFRSDNeb_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_SFRSD__Tg[iT][m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ \Sigma_{SFR}^\star\ [M_\odot yr^{-1} pc^{-2}]$',
                                    tSF, '%s/tauV_tauVNeb_logSFRSD_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))
            x = ALL_tau_V__Tg[iT][m]
            y = ALL_tauVNeb__g[m]
            z = np.log10(ALL_WHa__g[m])
            plot_3din2d_scatter_age(x, y, z, 
                                    r'$\tau_V$', r'$\tau_V^{neb}$', r'$\log\ W_{H\alpha}\ [\AA]$',
                                    tSF, '%s/tauV_tauVNeb_logWHa_age_%sMyr2.png' % (tauFilteredDir, str(tSF / 1.e6)))
            
            x = ALL_alogSFRSD__Trg[iT, :, :].flatten()
            y = ALL_alogSFRSD_Ha__rg[:, :].flatten()
            xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr.png' % str(tSF / 1.e6)
            xlim = np.percentile(x, [1, 100 * (x.flatten().shape[0] - x.mask.sum()) / x.flatten().shape[0] - 1])
            ylim = np.percentile(y, [1, 100 * (y.flatten().shape[0] - y.mask.sum()) / y.flatten().shape[0] - 1])
            plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
            
            x = np.log10(ALL_SFR__Tg[iT])
            y = np.log10(ALL_SFR_Ha__g)
            xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
            ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
            fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
            xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
            ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
            plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
            
            x = np.log10(ALL_SFR__Tg[iT])
            y = np.log10(ALL_SFR_Ha__g)
            xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
            ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
            fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
            xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
            ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
            plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
            
            x = ALL_tau_V__Trg[iT, :, :].flatten()
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
            
            x = ALL_tau_V_neb__rg.flatten()
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
            
            x = np.log10(ALL_tau_V__Trg[iT, :, :].flatten())
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\log\ \tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
            
            x = np.log10(ALL_tau_V_neb__rg[:, :].flatten())
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\log\ \tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 