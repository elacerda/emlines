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
from lines import *
import os

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

plot = True
#plot = False
debug = False
#debug = True
BPTLowS06 = True
#BPTLowS06 = False

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

RColor = [ 'r', 'y', 'b', 'k' ]
RRange = [  .5,  1., 1.5, 2.  ]

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
tauVNebOkMin = 0.05

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

def radialProfileWeighted(v__yx, w__yx, bins, rad_scale, func_radialProfile = None):
    v__r = None

    if func_radialProfile:
        w__r = func_radialProfile(w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v_w__r = func_radialProfile(v__yx * w__yx, bin_r = bins, mode = 'sum', rad_scale = rad_scale)
        v__r = v_w__r / w__r

    return v__r

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

def Mpc_to_cm(dist):
    # google: 1 Mpc = 3.08567758e24 cm
    return dist * 3.08567758e24 

def flux_to_LSol(flux, distance):
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / Lsun

def calc_Lobs(f_obs__Lz, err_f_obs__Lz, distance_Mpc):
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

def calc_Lint_Ha(Lobs__Lz, err_Lobs__Lz, tauVNeb__z, lines):
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

def plotCid(x1, x1label, y1, y1label, x2, x2label, y2, y2label, fname):
    f, axArr = plt.subplots(1, 2)
    f.set_dpi(96)
    f.set_size_inches(20,10)    
    ax = axArr[0]
    scat = ax.scatter(x1, y1, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.set_xlabel(x1label)
    ax.set_ylabel(y1label)
    ax = axArr[1]
    scat = ax.scatter(x2, y2, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.set_xlabel(x2label)
    ax.set_ylabel(y2label)
    plt.tight_layout()
    f.savefig(fname)
    
def plotRunningStatsAxis(ax, x, y, color):
    nBox        = 25
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
    ax.plot(xMean, yMean, 'o-', c = color, lw = 2)
    ax.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = color)

if __name__ == '__main__':
    #########################################################################        
    ################################## CID ##################################
    ######################################################################### 
    ALL_morfType_GAL__g     = np.ma.zeros((N_gals))
    ALL_at_flux_GAL__g      = np.ma.zeros((N_gals))
    ALL_Mcor_GAL__g         = np.ma.zeros((N_gals))
    ALL_McorSD_GAL__g       = np.ma.zeros((N_gals))
    #########################################################################
    #########################################################################
    #########################################################################

    #########################################################################
    ################## VARIABLES WITH shape N_Zone * N_gals #################
    #########################################################################
    _ALL_morfType_GAL_zones__g  = []
    _ALL_Mcor_GAL_zones__g      = []
    _ALL_McorSD_GAL_zones__g    = []
    _ALL_tauVNeb__g             = []
    _ALL_tauVNeb_err__g         = []
    _ALL_tauVNeb_mask__g        = []
    _ALL_SFRNeb__g              = []
    _ALL_SFRSDNeb__g            = []
    _ALL_SFR_Ha__g              = []
    _ALL_SFRSD_Ha__g            = []
    _ALL_L_int_Ha__g            = []
    _ALL_L_int_Ha_err__g        = []
    _ALL_L_int_Ha_mask__g       = []
    #########################################################################
    #########################################################################
    #########################################################################
    
    #########################################################################
    ############### VARIABLES WITH shape N_T * N_Zone * N_gals ##############
    #########################################################################
    _ALL_tauV__Tg               = []
    _ALL_tauV_mask__Tg          = []
    _ALL_SFR__Tg                = []
    _ALL_SFRSD__Tg              = []
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_morfType_GAL_zones__rg  = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD_Ha__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_alogSFRSD_Ha__rg        = np.ma.zeros((NRbins, N_gals))
    ALL_tauVneb__rg             = np.ma.zeros((NRbins, N_gals))
    ALL_alogZ_mass_GAL__Tg      = np.ma.zeros((N_T, N_gals))
    ALL_alogZ_flux_GAL__Tg      = np.ma.zeros((N_T, N_gals))
    ALL_isOkFrac_GAL__Tg        = np.ma.zeros((N_T, N_gals))
    ALL_aSFRSD__Trg             = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogSFRSD__Trg          = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_tauV__Trg               = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_mass__Trg         = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_flux__Trg         = np.ma.zeros((N_T, NRbins, N_gals))
    
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
        
        f_obs__Lz = K.EL.flux
        f_obs_err__Lz = K.EL.flux
        L_obs__Lz, L_obs_err__Lz = calc_Lobs(K.EL.flux, K.EL.eflux, K.distance_Mpc)
        L_int_Ha__z, L_int_Ha_err__z = calc_Lint_Ha(L_obs__Lz, L_obs_err__Lz, K.EL.tau_V_neb__z, K.EL.lines)
         
        ##################### MASK ######################
        # minimum value of f_lz / err_f_lz

        minSNR = 3.
        
        i_Hb = K.EL.lines.index('4861')
        i_O3 = K.EL.lines.index('5007')
        i_Ha = K.EL.lines.index('6563')
        i_N2 = K.EL.lines.index('6583')
        
        HbOk = (K.EL.flux[i_Hb, :] / K.EL.eflux[i_Hb, :]) >= minSNR
        O3Ok = (K.EL.flux[i_O3, :] / K.EL.eflux[i_O3, :]) >= minSNR
        HaOk = (K.EL.flux[i_Ha, :] / K.EL.eflux[i_Ha, :]) >= minSNR
        N2Ok = (K.EL.flux[i_N2, :] / K.EL.eflux[i_N2, :]) >= minSNR
        
        N2Ha = np.log10(K.EL.N2_obs__z/K.EL.Ha_obs__z)
        O3Hb = np.log10(K.EL.O3_obs__z/K.EL.Hb_obs__z)
        
        maskOk = HbOk & O3Ok & HaOk & N2Ok
        maskBPT = None
        
        if BPTLowS06:
            L = Lines()
            maskBPT = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
            maskOk &= maskBPT

        f_obs__Lz[:, ~maskOk] = np.ma.masked
        L_obs__Lz[:, ~maskOk] = np.ma.masked
        L_int_Ha__z[~maskOk] = np.ma.masked
        #################################################

        #########################################################################        
        ################################## CID ##################################
        ######################################################################### 
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        ALL_McorSD_GAL__g[iGal] = K.McorSD__yx.mean()
        aux = np.ones_like(K.Mcor__z) * ALL_McorSD_GAL__g[iGal]
        _ALL_McorSD_GAL_zones__g.append(aux)
        _ALL_Mcor_GAL_zones__g.append(K.Mcor__z)
        ALL_Mcor_GAL__g[iGal]   = K.Mcor_tot.sum()

        # Compute & store galaxy-wide at_flux
        numerator__z   = K.Lobn__tZz.sum(axis=1).sum(axis=0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis=1).sum(axis=0)
        ALL_at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()
        
        for iT, tSF in enumerate(tSF__T):
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
                 
                if BPTLowS06:
                    maskNotOk__z |= ~L.maskBelowlinebpt('S06', N2Ha, O3Hb)

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
            alogSFRSD__r = K.radialProfile(np.log10(SFRSD__yx * 1e6), Rbin__r, rad_scale = K.HLR_pix)
            ALL_aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
            ALL_alogSFRSD__Trg[iT, :, iGal] = alogSFRSD__r

            _ALL_tauV__Tg[iT].append(tauV__z.data)
            _ALL_tauV_mask__Tg[iT].append(tauV__z.mask)
            
            tauV__yx = K.zoneToYX(tauV__z, extensive=False)
            ALL_tauV__Trg[iT, :, iGal] = K.radialProfile(tauV__yx, Rbin__r, rad_scale = K.HLR_pix)
            
            _ALL_SFR__Tg[iT].append(SFR__z)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z)
        #########################################################################
        #########################################################################
        #########################################################################
        
        maskOkTauVNeb = (K.EL.tau_V_neb__z >= tauVNebOkMin) & (K.EL.tau_V_neb_err__z <= 0.15) 
        
        tauVNeb__z = K.EL.tau_V_neb__z
        tauVNeb__z[~maskOkTauVNeb] = np.ma.masked
        tauVNeb_err__z = K.EL.tau_V_neb_err__z
        tauVNeb_err__z[~maskOkTauVNeb] = np.ma.masked

        _ALL_tauVNeb__g.append(tauVNeb__z.data)
        _ALL_tauVNeb_err__g.append(tauVNeb_err__z.data)
        _ALL_tauVNeb_mask__g.append(tauVNeb__z.mask)
        
        _ALL_L_int_Ha__g.append(L_int_Ha__z.data)
        _ALL_L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        _ALL_L_int_Ha_mask__g.append(L_int_Ha__z.mask)

        tauVNeb__yx = K.zoneToYX(tauVNeb__z, extensive = False)
        ALL_tauVneb__rg[:, iGal] = K.radialProfile(tauVNeb__yx, Rbin__r, rad_scale = K.HLR_pix)
        
        # 3.17 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter        
        SFR_Ha__z       = 3.17 * L_int_Ha__z / (1.e8)
        SFRSD_Ha__z     = SFR_Ha__z / K.zoneArea_pc2
        
        SFRSD_Ha__yx    = K.zoneToYX(SFR_Ha__z, extensive = True)
        aSFRSD_Ha__r    = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        alogSFRSD_Ha__r = K.radialProfile(np.log10(SFRSD_Ha__yx * 1e6), Rbin__r, rad_scale = K.HLR_pix)
        ALL_aSFRSD_Ha__rg[:, iGal]    = aSFRSD_Ha__r
        ALL_alogSFRSD_Ha__rg[:, iGal] = alogSFRSD_Ha__r
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        
        K.close()
        
    aux = np.hstack(np.asarray(_ALL_tauVNeb__g))
    auxMask = np.hstack(np.asarray(_ALL_tauVNeb_mask__g))
    ALL_tauVNeb__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_tauVNeb_err__g = np.ma.masked_array(np.hstack(np.asarray(_ALL_tauVNeb_err__g)), mask = auxMask)

    aux                         = np.hstack(np.asarray(_ALL_L_int_Ha__g))
    auxMask                     = np.hstack(np.asarray(_ALL_L_int_Ha_mask__g))
    ALL_L_int_Ha__g             = np.ma.masked_array(aux, mask = auxMask)
    ALL_SFR_Ha__g               = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR_Ha__g)), mask = auxMask)
    ALL_SFRSD_Ha__g             = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_Ha__g)), mask = auxMask)
    
    ALL_Mcor_GAL_zones__g       = np.ma.masked_array(np.hstack(np.asarray(_ALL_Mcor_GAL_zones__g)))
    ALL_McorSD_GAL_zones__g     = np.ma.masked_array(np.hstack(np.asarray(_ALL_McorSD_GAL_zones__g)))
    ALL_morfType_GAL_zones__g   = np.ma.masked_array(np.hstack(np.asarray(_ALL_morfType_GAL_zones__g)))
    
    ALL_tauV__Tg     = []
    ALL_SFR__Tg      = []
    ALL_SFRSD__Tg    = []

    for iT,tSF in enumerate(tSF__T):
        aux     = np.hstack(np.asarray(_ALL_tauV__Tg[iT]))
        auxMask = np.hstack(np.asarray(_ALL_tauV_mask__Tg[iT]))
        ALL_tauV__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_SFR__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR__Tg[iT])), mask = auxMask))
        ALL_SFRSD__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD__Tg[iT])), mask = auxMask))
        
    if plot:
        for iT,tSF in enumerate(tSF__T):
            x1 = np.log10(ALL_SFR__Tg[iT])
            x1label = r'$\log\ \mathrm{SFR}_\star\ [M_\odot yr^{-1}]$' 
            y1 = np.log10(ALL_SFR_Ha__g)
            y1label = r'$\log\ \mathrm{SFR}_{neb}\ [M_\odot yr^{-1}]$' 
            x2 = ALL_alogSFRSD__Trg[iT, :, :].flatten()
            x2label = r'$\log\ \Sigma_{\mathrm{SFR}}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            y2 = ALL_alogSFRSD_Ha__rg.flatten()
            y2label = r'$\log\ \Sigma_{\mathrm{SFR}}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            fname = 'SFReSFRSD_%sMyr.png' % str(tSF / 1.e6)
            plotCid(x1, x1label, y1, y1label, x2, x2label, y2, y2label, fname)

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
             
            x = ALL_tauV__Trg[iT, :, :].flatten()
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
             
            x = ALL_tauVneb__rg.flatten()
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
             
            x = np.log10(ALL_tauV__Trg[iT, :, :].flatten())
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\log\ \tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
             
            x = np.log10(ALL_tauVneb__rg[:, :].flatten())
            y = ALL_alogSFRSD_Ha__rg.flatten() - ALL_alogSFRSD__Trg[iT, :, :].flatten()
            xlabel = r'$\log\ \tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
         
            ###################### RADIUS COLOR ######################
            xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
                    
                x = ALL_alogSFRSD__Trg[iT, RMask, :].flatten()
                y = ALL_alogSFRSD_Ha__rg[RMask, :].flatten()
                scat = ax.scatter(x, y, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, x, y, RColor[iR])
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
            rhoSpearman, pvalSpearman = st.spearmanr(x, y)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
     
            xlabel = r'$\tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
                x = ALL_tauV__Trg[iT, RMask, :].flatten()
                y = ALL_alogSFRSD_Ha__rg[RMask, :].flatten() - ALL_alogSFRSD__Trg[iT, RMask, :].flatten()
                scat = ax.scatter(x, y, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, x, y, RColor[iR])
            rhoSpearman, pvalSpearman = st.spearmanr(x, y)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
             
            xlabel = r'$\tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
                x = ALL_tauVneb__rg[RMask, :].flatten()
                y = ALL_alogSFRSD_Ha__rg[RMask, :].flatten() - ALL_alogSFRSD__Trg[iT, RMask, :].flatten()
                scat = ax.scatter(x, y, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, x, y, RColor[iR])
            rhoSpearman, pvalSpearman = st.spearmanr(x, y)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
             
            xlabel = r'$\log\ \tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
                x = np.log10(ALL_tauV__Trg[iT, RMask, :].flatten())
                y = ALL_alogSFRSD_Ha__rg[RMask, :].flatten() - ALL_alogSFRSD__Trg[iT, RMask, :].flatten()
                scat = ax.scatter(x, y, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, x, y, RColor[iR])
            rhoSpearman, pvalSpearman = st.spearmanr(x, y)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
             
            xlabel = r'$\log\ \tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
                x = np.log10(ALL_tauVneb__rg[RMask, :].flatten())
                y = ALL_alogSFRSD_Ha__rg[RMask, :].flatten() - ALL_alogSFRSD__Trg[iT, RMask, :].flatten()
                scat = ax.scatter(x, y, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, x, y, RColor[iR])
            rhoSpearman, pvalSpearman = st.spearmanr(x, y)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
