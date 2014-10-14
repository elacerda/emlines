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
import scipy.optimize as so
from lines import *
import os
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
from matplotlib.ticker import MultipleLocator, MaxNLocator
from astropy import units as u

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#nebular_plot = True
nebular_plot = False
#plot = True
plot = False
debug = False
#debug = True
#BPTLowS06 = True
BPTLowS06 = False
#hdf5 = False
hdf5 = True

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
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
    listOfPrefixes = listOfPrefixes[1:10]        # Q&D tests ...
    #listOfPrefixes = ['K0277\n']
    
N_gals = len(listOfPrefixes)

# SFR-time-scale array (index __T)
base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
tSF__T = base.ageBase
N_T = base.nAges
tSF_to_plot = [0, 10, 14, 17, 20, 23, 26, 29, 32, 35, 39 ]

mask_xOk = True

# Def smallest light fraction (in the flag__t-ageMax age-range) deemed to be Ok for our stats ...
xOkMin = 0.05

# Minimum tauV to be taken seriously ...
tauVOkMin = 0.05
tauVNebOkMin = 0.05
tauVNebErrMax = 0.15

def get_attrib_h5(h5, attrib):
    if any([ attrib in s for s in h5['masked/mask'].keys() ]):
        node = '/masked/data/' + attrib
        ds = h5[node]
        if type(ds) == h5py.Dataset:
            data = h5.get('/masked/data/' + attrib).value
            mask = h5.get('/masked/mask/' + attrib).value
            arr = np.ma.masked_array(data, mask = mask)
        else:
            arr = []
            tSF__T = h5.get('/data/tSF__T').value
            for iT,tSF in enumerate(tSF__T):
                group = '%s/%d' % (attrib, iT)
                data = h5.get('/masked/data/' + group).value
                mask = h5.get('/masked/mask/' + group).value
                arr.append(np.ma.masked_array(data, mask = mask))
        return arr 
    else:
        return h5.get('/data/' + attrib).value

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level
 
def density_contour(xdata, ydata, binsx, binsy, ax=None, **contour_kwargs):
    """ Create a density contour plot.
 
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """    
    nbins_x = len(binsx) - 1
    nbins_y = len(binsy) - 1

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=[binsx,binsy],  normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
 
    pdf = (H*(x_bin_sizes*y_bin_sizes))
 
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [one_sigma, two_sigma, three_sigma]
 
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
 
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
 
    return contour

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
    f.set_size_inches(10,8)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    yxlabel = r'$%s /\ %s $ ' % (ylabel.split('[')[0].strip('$ '), xlabel.split('[')[0].strip('$ '))
    txt = '%s mean:%.3f  median:%.3f  $\sigma(y/x)$:%.3f  Rs: %.2f' % (yxlabel, (y/x).mean(), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    print age, txt
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.03, 0.97, txt, fontsize = 16, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
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
    f.set_size_inches(10,8)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    #ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    yxlabel = r'$%s /\ %s $ ' % (ylabel.split('[')[0].strip('$ '), xlabel.split('[')[0].strip('$ '))
    txt = '%s mean:%.3f  median:%.3f  $\sigma(y/x)$:%.3f  Rs: %.2f' % (yxlabel, (y/x).mean(), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.3, 0.97, txt, fontsize = 15, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
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
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / L_sun

def calc_Lint_Ha(L_obs__Lz, L_obs_err__Lz, tau_V_neb__z, lines):
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    
    eHa = np.ma.exp(qCCM['6563'] * tau_V_neb__z)
    L_obs_HaHb__z = L_obs__Lz[i_Ha, :] / L_obs__Lz[i_Hb, :]

    L_int_Ha__z = L_obs__Lz[i_Ha, :] * eHa
    
    a = L_obs_err__Lz[i_Ha, :]
    b = q * L_obs_HaHb__z * L_obs_err__Lz[i_Hb, :]
    L_int_Ha_err__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
    
    return L_int_Ha__z, L_int_Ha_err__z

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
    plt.close(f)
 
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
    _ALL_tau_V_neb__g           = []
    _ALL_tau_V_neb_err__g       = []
    _ALL_tau_V_neb_mask__g      = []
    _ALL_SFR_Ha__g              = []
    _ALL_SFRSD_Ha__g            = []
    _ALL_SFRSD_Ha_kpc__g        = []
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
    _ALL_SFRSD_kpc__Tg          = []
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_morfType_GAL_zones__rg  = np.ma.zeros((NRbins, N_gals))
    ALL_tau_V_neb__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_alogZ_mass_GAL__Tg      = np.ma.zeros((N_T, N_gals))
    ALL_alogZ_flux_GAL__Tg      = np.ma.zeros((N_T, N_gals))
    ALL_isOkFrac_GAL__Tg        = np.ma.zeros((N_T, N_gals))

    ALL_aSFRSD_Ha__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD_Ha_kpc__rg       = np.ma.zeros((NRbins, N_gals))
    
    ALL_aSFRSD__Trg             = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_aSFRSD_kpc__Trg         = np.ma.zeros((N_T, NRbins, N_gals))
    
    ALL_tauV__Trg               = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_mass__Trg         = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_flux__Trg         = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_McorSD_GAL__rg          = np.ma.zeros((NRbins, N_gals))
    
    ALL_zones_tau_V = 0
    ALL_zones_Lum   = 0
    ALL_gals        = 0
    ALL_zones       = 0
    
    for iT in range(N_T):
        _ALL_tauV__Tg.append([])
        _ALL_tauV_mask__Tg.append([])
        _ALL_SFR__Tg.append([])
        _ALL_SFRSD__Tg.append([])
        _ALL_SFRSD_kpc__Tg.append([])
        
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
            print '<<< %s galaxy: miss files' & galName 
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
        ALL_McorSD_GAL__rg[:, iGal] = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)
        
        ##########################
        ###### MASK EmLines ######
        ##########################        
        # minimum value of f_lz / err_f_lz
        minSNR = 3.
        
        i_Hb = K.EL.lines.index('4861')
        i_O3 = K.EL.lines.index('5007')
        i_Ha = K.EL.lines.index('6563')
        i_N2 = K.EL.lines.index('6583')
        
        HbOk = np.array((K.EL.flux[i_Hb, :] / K.EL.eflux[i_Hb, :]) >= minSNR, dtype = np.bool)
        O3Ok = np.array((K.EL.flux[i_O3, :] / K.EL.eflux[i_O3, :]) >= minSNR, dtype = np.bool)
        HaOk = np.array((K.EL.flux[i_Ha, :] / K.EL.eflux[i_Ha, :]) >= minSNR, dtype = np.bool)
        N2Ok = np.array((K.EL.flux[i_N2, :] / K.EL.eflux[i_N2, :]) >= minSNR, dtype = np.bool)
        
        N2Ha = np.ma.log10(K.EL.N2_obs__z/K.EL.Ha_obs__z)
        O3Hb = np.ma.log10(K.EL.O3_obs__z/K.EL.Hb_obs__z)
        
        maskOk = HbOk & O3Ok & HaOk & N2Ok
        maskBPT = None
        
        if BPTLowS06:
            L = Lines()
            maskBPT = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        ##########################
        ##########################
        ##########################
        
        K.EL._forceMask = ~maskOk
        
        if len(K.EL.flux[0,:].compressed()) == 0:
            ALL_morfType_GAL__g[iGal]           = np.ma.masked
            ALL_at_flux_GAL__g[iGal]            = np.ma.masked
            ALL_Mcor_GAL__g[iGal]               = np.ma.masked
            ALL_McorSD_GAL__g[iGal]             = np.ma.masked
            print '<<< %s galaxy: no minSNR (%.1f) in 4 BPT lines' % (galName, minSNR)
            continue
        
        print '>>> Doing' , iGal , galName , 'hubtyp=', ALL_morfType_GAL__g[iGal], '|  Nzones=' , K.N_zone
        ALL_zones += K.N_zone
        ALL_gals += 1

        K.EL._forceMask = None
        
        ##########################
        ####### STARLIGHT ########
        ##########################        
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
            aux             = K.Mini__tZz[flag__t,:,:].sum(axis=1).sum(axis=0) / tSF
            SFR__z          = np.ma.masked_array(aux)
            SFRSD__z        = SFR__z / K.zoneArea_pc2
            SFRSD_kpc__z    = SFRSD__z * 1e6
            
            tauV__z     = np.ma.masked_array(K.tau_V__z)
                            
            if mask_xOk:
                # Compute xOk "raw" image
                x__tZz  =  K.popx / K.popx.sum(axis=1).sum(axis=0)
                xOk__z  = x__tZz[flag__t,:,:].sum(axis=1).sum(axis=0)
                
                maskNotOk__z = (xOk__z < xOkMin) | (tauV__z < tauVOkMin) 
                 
                if BPTLowS06:
                    maskNotOk__z |= ~maskBPT

                tauV__z[maskNotOk__z]   = np.ma.masked
                SFR__z[maskNotOk__z]    = np.ma.masked
                SFRSD__z[maskNotOk__z]  = np.ma.masked
                SFRSD_kpc__z[maskNotOk__z]  = np.ma.masked

            weiRadProf = True
                            
            aux = calc_alogZ_Stuff(K, flag__t, Rbin__r, weiRadProf, xOkMin = xOkMin)
            ALL_alogZ_mass_GAL__Tg[iT, iGal] = aux[0]
            ALL_alogZ_flux_GAL__Tg[iT, iGal] = aux[1]
            ALL_isOkFrac_GAL__Tg[iT, iGal] = aux[2]
            ALL_alogZ_mass__Trg[iT, :, iGal] = aux[3]
            ALL_alogZ_flux__Trg[iT, :, iGal] = aux[4]

            aSFRSD__yx      = K.zoneToYX(SFR__z, extensive = True)
            aSFRSD__r       = K.radialProfile(aSFRSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            aSFRSD_kpc__r       = K.radialProfile(aSFRSD__yx * 1e6, Rbin__r, rad_scale = K.HLR_pix)

            
            ALL_aSFRSD__Trg[iT, :, iGal]        = aSFRSD__r
            ALL_aSFRSD_kpc__Trg[iT, :, iGal]    = aSFRSD_kpc__r

            _ALL_tauV__Tg[iT].append(tauV__z.data)
            _ALL_tauV_mask__Tg[iT].append(tauV__z.mask)
            
            tauV__yx = K.zoneToYX(tauV__z, extensive=False)
            ALL_tauV__Trg[iT, :, iGal] = K.radialProfile(tauV__yx, Rbin__r, rad_scale = K.HLR_pix)
            
            _ALL_SFR__Tg[iT].append(SFR__z)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z)
            _ALL_SFRSD_kpc__Tg[iT].append(SFRSD__z)
        ##########################
        ##########################
        ##########################    

        ##########################
        ########## tau_V #########
        ##########################
        if BPTLowS06:
            K.EL._forceMask = ~(maskOk & maskBPT) #Changing global EL mask
        else:
            K.EL._forceMask = ~maskOk #Changing global EL mask
            
        maskOkTauVNeb = (K.EL.tau_V_neb__z >= tauVNebOkMin) & (K.EL.tau_V_neb_err__z <= tauVNebErrMax)
        mask_temp = maskOk & maskOkTauVNeb
        K.EL._forceMask = ~mask_temp #Changing global EL mask 
        
        N_zones_tau_V = len(K.EL.tau_V_neb__z.compressed())
        print 'tauV calculated for %d zones (maskOK and maskOkTauVNeb)' % N_zones_tau_V
        ALL_zones_tau_V += N_zones_tau_V
        
        tau_V_neb__z                = K.EL.tau_V_neb__z
        tau_V_neb_err__z            = K.EL.tau_V_neb_err__z
        tau_V_neb__yx               = K.zoneToYX(tau_V_neb__z, extensive = False)
        tau_V_neb_err__yx           = K.zoneToYX(tau_V_neb_err__z, extensive = False)
        tau_V_neb__r                = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)

        ALL_tau_V_neb__rg[:, iGal]  = tau_V_neb__r
        _ALL_tau_V_neb__g.append(tau_V_neb__z.data)
        _ALL_tau_V_neb_err__g.append(tau_V_neb_err__z.data)
        _ALL_tau_V_neb_mask__g.append(tau_V_neb__z.mask)
        ##########################
        ##########################
        ##########################
        
        ##########################
        #### SFR and SigmaSFR ####
        ##########################
        N_zones_Lum = len(K.EL.flux[0,:].compressed())
        print 'Luminosity and products: minSNR mask: using %d zones' % N_zones_Lum
            
        ALL_zones_Lum += N_zones_Lum
            
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux)  / L_sun
        L_int_Ha__z, L_int_Ha_err__z = calc_Lint_Ha(L_obs__Lz, L_obs_err__Lz, K.EL.tau_V_neb__z, K.EL.lines)
        
        _ALL_L_int_Ha__g.append(L_int_Ha__z.data)
        _ALL_L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        _ALL_L_int_Ha_mask__g.append(L_int_Ha__z.mask)
                
        # 3.17 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter        
        SFR_Ha__z       = 3.17 * L_int_Ha__z / (1.e8)
        SFRSD_Ha__z     = SFR_Ha__z / K.zoneArea_pc2
        SFRSD_Ha_kpc__z = SFRSD_Ha__z * 1e6
        
        SFRSD_Ha__yx        = K.zoneToYX(SFR_Ha__z, extensive = True)
        aSFRSD_Ha__r        = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        aSFRSD_Ha_kpc__r    = K.radialProfile(SFRSD_Ha__yx * 1e6, Rbin__r, rad_scale = K.HLR_pix)

        ALL_aSFRSD_Ha__rg[:, iGal]          = aSFRSD_Ha__r
        ALL_aSFRSD_Ha_kpc__rg[:, iGal]      = aSFRSD_Ha_kpc__r
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        _ALL_SFRSD_Ha_kpc__g.append(SFRSD_Ha__z)
        ##########################
        ##########################
        ##########################
        
        if nebular_plot:
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
            ax.plot(RbinCenter__r, tau_V_neb__r, 'o-k')
            #ax.tick_params(axis='x', pad=30)
            ax.set_title(r'$\tau_V^{neb}(R)$')
            
            ax                  = axArr[0, 2]
            ax.set_axis_on()
            L_int_Ha__yx        = K.zoneToYX(L_int_Ha__z, extensive = True)
            L_int_Ha__r         = K.radialProfile(L_int_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
            ax.plot(RbinCenter__r, L_int_Ha__r, 'o-k')
            ax.set_title(r'$L_{H\alpha}^{int}(R)$')
        
            ax                  = axArr[0, 3]
            ax.set_axis_on()
            #Lobn__yx            = K.zoneToYX(K.Lobn__z, extensive = True)
            #logZNeb__r          = radialProfileWeighted(logZ_neb__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
            logZ_neb__yx         = K.zoneToYX(K.EL.logZ_neb_S06__z, extensive = False)
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
            sigma               = tau_V_neb_err__yx.std()
            mean                = tau_V_neb_err__yx.mean()
            vmax                = mean + 2. * sigma
            im                  = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax                  = axArr[3, 1]
            ax.set_axis_on()
            ax.set_title(r'$F^{H\alpha}_{H\beta} / \epsilon(F^{H\alpha}_{H\beta})$')
            HaHb__yx            = K.zoneToYX(K.EL.HaHb__z, extensive = True)
            err_HaHb__yx        = K.zoneToYX(K.EL.HaHb_err__z, extensive = True)
            signalToNoise       = np.abs(HaHb__yx) / np.abs(err_HaHb__yx) 
            sigma               = signalToNoise.std()
            mean                = signalToNoise.mean()
            vmax                = mean + 2. * sigma
            im                  = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
        
            ax                  = axArr[1, 2]
            ax.set_axis_on()
            ax.set_title(r'$L_{H\alpha}^{int}$')
            sigma               = L_int_Ha__yx.std()
            mean                = L_int_Ha__yx.mean()
            vmax                = mean + sigma
            vmin                = mean - sigma
            im                  = ax.imshow(L_int_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax                  = axArr[2, 2]
            ax.set_axis_on()
            ax.set_title(r'$\epsilon (L_{H\alpha}^{int})$')
            L_int_Ha_err__yx    = K.zoneToYX(L_int_Ha_err__z, extensive = True)
            sigma               = L_int_Ha_err__yx.std()
            mean                = L_int_Ha_err__yx.mean()
            vmax                = mean + sigma
            im                  = ax.imshow(L_int_Ha_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax                  = axArr[3, 2]
            ax.set_axis_on()
            ax.set_title(r'$L_{H\alpha}^{int} / \epsilon(L_{H\alpha}^{int})$')
            signalToNoise       = np.abs(L_int_Ha__yx) / np.abs(L_int_Ha_err__yx) 
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
            logZ_neb_err__yx    = K.zoneToYX(K.EL.logZ_neb_S06_err__z, extensive = False)
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
            f.savefig(K.califaID + '_' + versionSuffix + '_nebular.png')
            plt.close(f)

        K.close()

    print 'Total of %d galaxies (%d zones): %s zones for tau_V and %d zones for SFR' % (ALL_gals, ALL_zones, ALL_zones_tau_V, ALL_zones_Lum)
    
    aux = np.hstack(np.asarray(_ALL_tau_V_neb__g))
    auxMask = np.hstack(np.asarray(_ALL_tau_V_neb_mask__g))
    ALL_tauVNeb__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_tauVNeb_err__g = np.ma.masked_array(np.hstack(np.asarray(_ALL_tau_V_neb_err__g)), mask = auxMask)

    aux                         = np.hstack(np.asarray(_ALL_L_int_Ha__g))
    auxMask                     = np.hstack(np.asarray(_ALL_L_int_Ha_mask__g))
    ALL_L_int_Ha__g             = np.ma.masked_array(aux, mask = auxMask)
    ALL_SFR_Ha__g               = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR_Ha__g)), mask = auxMask)
    ALL_SFRSD_Ha__g             = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_Ha__g)), mask = auxMask)
    ALL_SFRSD_Ha_kpc__g         = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_Ha_kpc__g)), mask = auxMask)
    
    ALL_Mcor_GAL_zones__g       = np.ma.masked_array(np.hstack(np.asarray(_ALL_Mcor_GAL_zones__g)))
    ALL_McorSD_GAL_zones__g     = np.ma.masked_array(np.hstack(np.asarray(_ALL_McorSD_GAL_zones__g)))
    ALL_morfType_GAL_zones__g   = np.ma.masked_array(np.hstack(np.asarray(_ALL_morfType_GAL_zones__g)))

    ALL_tauV__Tg        = []
    ALL_SFR__Tg         = []
    ALL_SFRSD__Tg       = []
    ALL_SFRSD_kpc__Tg   = []
    correl_SFR = np.ones_like(tSF__T)
    correl_SFRSD__rT = np.ones((1 + len(RRange), tSF__T.shape[0]))
    correl_SFRSD_kpc__rT = np.ones((1 + len(RRange), tSF__T.shape[0])) 

    for iT,tSF in enumerate(tSF__T):
        aux     = np.hstack(np.asarray(_ALL_tauV__Tg[iT]))
        auxMask = np.hstack(np.asarray(_ALL_tauV_mask__Tg[iT]))
        ALL_tauV__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_SFR__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFR__Tg[iT])), mask = auxMask))
        ALL_SFRSD__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD__Tg[iT])), mask = auxMask))
        ALL_SFRSD_kpc__Tg.append(np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_kpc__Tg[iT])), mask = auxMask))

        x = np.ma.log10(ALL_SFR__Tg[iT])
        y = np.ma.log10(ALL_SFR_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_SFR[iT] = st.spearmanr(xm, ym)[0]

        x = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
        y = np.ma.log10(ALL_aSFRSD_Ha__rg[:, :].flatten())
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_SFRSD__rT[0, iT] = st.spearmanr(xm, ym)[0]
        
        for iR, RUp in enumerate(RRange):
            if iR == 0:
                RMask = RbinCenter__r <= RUp
            else:
                RDown =  RRange[iR - 1]
                RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                
            iiR = iR + 1

            x = np.ma.log10(ALL_aSFRSD__Trg[iT, RMask, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten())
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            correl_SFRSD__rT[iiR, iT] = st.spearmanr(xm, ym)[0]

        x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
        y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[:, :].flatten())
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_SFRSD_kpc__rT[0, iT] = st.spearmanr(xm, ym)[0]
        
        for iR, RUp in enumerate(RRange):
            if iR == 0:
                RMask = RbinCenter__r <= RUp
            else:
                RDown =  RRange[iR - 1]
                RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                
            iiR = iR + 1

            x = ALL_aSFRSD_kpc__Trg[iT, RMask, :].flatten()
            y = ALL_aSFRSD_Ha_kpc__rg[RMask, :].flatten()
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            correl_SFRSD_kpc__rT[iiR, iT] = st.spearmanr(xm, ym)[0]








    if hdf5:
        import h5py
        
            
                
        if BPTLowS06:       
            filename = 'SFR_BelowS06.h5'
        else:
            filename = 'SFR.h5'
            
        file = h5py.File(filename, 'w')
        
        D = {
            '/masked/data/ALL_morfType_GAL__g' : ALL_morfType_GAL__g.data,
            '/masked/data/ALL_at_flux_GAL__g' : ALL_at_flux_GAL__g.data,
            '/masked/data/ALL_Mcor_GAL__g' : ALL_Mcor_GAL__g.data,  
            '/masked/data/ALL_McorSD_GAL__g' : ALL_McorSD_GAL__g.data,  
            '/masked/data/ALL_morfType_GAL_zones__rg' : ALL_morfType_GAL_zones__rg.data,  
            '/masked/data/ALL_tau_V_neb__rg' : ALL_tau_V_neb__rg.data, 
            '/masked/data/ALL_alogZ_mass_GAL__Tg' : ALL_alogZ_mass_GAL__Tg.data, 
            '/masked/data/ALL_alogZ_flux_GAL__Tg' : ALL_alogZ_flux_GAL__Tg.data, 
            '/masked/data/ALL_isOkFrac_GAL__Tg' : ALL_isOkFrac_GAL__Tg.data, 
            '/masked/data/ALL_aSFRSD_Ha__rg' : ALL_aSFRSD_Ha__rg.data, 
            '/masked/data/ALL_aSFRSD_Ha_kpc__rg' : ALL_aSFRSD_Ha_kpc__rg.data, 
            '/masked/data/ALL_aSFRSD__Trg' : ALL_aSFRSD__Trg.data, 
            '/masked/data/ALL_aSFRSD_kpc__Trg' : ALL_aSFRSD_kpc__Trg.data, 
            '/masked/data/ALL_tauV__Trg' : ALL_tauV__Trg.data, 
            '/masked/data/ALL_alogZ_mass__Trg' : ALL_alogZ_mass__Trg.data, 
            '/masked/data/ALL_alogZ_flux__Trg' : ALL_alogZ_flux__Trg.data, 
            '/masked/data/ALL_McorSD_GAL__rg' : ALL_McorSD_GAL__rg.data, 
            '/masked/data/ALL_tauVNeb__g' :  ALL_tauVNeb__g.data, 
            '/masked/data/ALL_tauVNeb_err__g' : ALL_tauVNeb_err__g.data, 
            '/masked/data/ALL_L_int_Ha__g' : ALL_L_int_Ha__g.data, 
            '/masked/data/ALL_SFR_Ha__g' : ALL_SFR_Ha__g.data, 
            '/masked/data/ALL_SFRSD_Ha__g' : ALL_SFRSD_Ha__g.data, 
            '/masked/data/ALL_SFRSD_Ha_kpc__g' : ALL_SFRSD_Ha_kpc__g.data, 
            '/masked/data/ALL_Mcor_GAL_zones__g' : ALL_Mcor_GAL_zones__g.data, 
            '/masked/data/ALL_McorSD_GAL_zones__g' : ALL_McorSD_GAL_zones__g.data, 
            '/masked/data/ALL_morfType_GAL_zones__g' : ALL_morfType_GAL_zones__g.data, 
            '/masked/mask/ALL_morfType_GAL__g' : ALL_morfType_GAL__g.mask, 
            '/masked/mask/ALL_at_flux_GAL__g' : ALL_at_flux_GAL__g.mask, 
            '/masked/mask/ALL_Mcor_GAL__g' : ALL_Mcor_GAL__g.mask,  
            '/masked/mask/ALL_McorSD_GAL__g' : ALL_McorSD_GAL__g.mask,  
            '/masked/mask/ALL_morfType_GAL_zones__rg' : ALL_morfType_GAL_zones__rg.mask,  
            '/masked/mask/ALL_tau_V_neb__rg' : ALL_tau_V_neb__rg.mask, 
            '/masked/mask/ALL_alogZ_mass_GAL__Tg' : ALL_alogZ_mass_GAL__Tg.mask, 
            '/masked/mask/ALL_alogZ_flux_GAL__Tg' : ALL_alogZ_flux_GAL__Tg.mask, 
            '/masked/mask/ALL_isOkFrac_GAL__Tg' : ALL_isOkFrac_GAL__Tg.mask, 
            '/masked/mask/ALL_aSFRSD_Ha__rg' : ALL_aSFRSD_Ha__rg.mask, 
            '/masked/mask/ALL_aSFRSD_Ha_kpc__rg' : ALL_aSFRSD_Ha_kpc__rg.mask, 
            '/masked/mask/ALL_aSFRSD__Trg' : ALL_aSFRSD__Trg.mask, 
            '/masked/mask/ALL_aSFRSD_kpc__Trg' : ALL_aSFRSD_kpc__Trg.mask, 
            '/masked/mask/ALL_tauV__Trg' : ALL_tauV__Trg.mask, 
            '/masked/mask/ALL_alogZ_mass__Trg' : ALL_alogZ_mass__Trg.mask, 
            '/masked/mask/ALL_alogZ_flux__Trg' : ALL_alogZ_flux__Trg.mask, 
            '/masked/mask/ALL_McorSD_GAL__rg' : ALL_McorSD_GAL__rg.mask, 
            '/masked/mask/ALL_tauVNeb__g' :  ALL_tauVNeb__g.mask, 
            '/masked/mask/ALL_tauVNeb_err__g' : ALL_tauVNeb_err__g.mask, 
            '/masked/mask/ALL_L_int_Ha__g' : ALL_L_int_Ha__g.mask, 
            '/masked/mask/ALL_SFR_Ha__g' : ALL_SFR_Ha__g.mask, 
            '/masked/mask/ALL_SFRSD_Ha__g' : ALL_SFRSD_Ha__g.mask, 
            '/masked/mask/ALL_SFRSD_Ha_kpc__g' : ALL_SFRSD_Ha_kpc__g.mask, 
            '/masked/mask/ALL_Mcor_GAL_zones__g' : ALL_Mcor_GAL_zones__g.mask, 
            '/masked/mask/ALL_McorSD_GAL_zones__g' : ALL_McorSD_GAL_zones__g.mask, 
            '/masked/mask/ALL_morfType_GAL_zones__g' : ALL_morfType_GAL_zones__g.mask, 
            '/data/correl_SFR' : correl_SFR, 
            '/data/correl_SFRSD__rT' : correl_SFRSD__rT, 
            '/data/correl_SFRSD_kpc__rT' : correl_SFRSD_kpc__rT, 
            '/data/ALL_zones_tau_V' : ALL_zones_tau_V, 
            '/data/ALL_zones_Lum' : ALL_zones_Lum, 
            '/data/ALL_gals' : ALL_gals, 
            '/data/ALL_zones' : ALL_zones, 
            '/data/RbinIni' : RbinIni, 
            '/data/RbinFin' : RbinFin, 
            '/data/RbinStep' : RbinStep, 
            '/data/Rbin__r' : Rbin__r, 
            '/data/RbinCenter__r' : RbinCenter__r, 
            '/data/NRbins' : NRbins, 
            '/data/RColor' : RColor, 
            '/data/RRange' : RRange, 
            '/data/N_gals' : N_gals, 
            '/data/tSF__T' : tSF__T, 
            '/data/xOkMin' : xOkMin, 
            '/data/tauVOkMin' : tauVOkMin, 
            '/data/tauVNebOkMin' : tauVNebOkMin, 
            '/data/tauVNebErrMax' : tauVNebErrMax,
        }

        for iT,tSF in enumerate(tSF__T):
            D['/masked/data/ALL_tauV__Tg/%d' % iT] = ALL_tauV__Tg[iT].data
            D['/masked/data/ALL_SFR__Tg/%d' % iT] = ALL_SFR__Tg[iT].data
            D['/masked/data/ALL_SFRSD__Tg/%d' % iT] = ALL_SFRSD__Tg[iT].data
            D['/masked/data/ALL_SFRSD_kpc__Tg/%d' % iT] = ALL_SFRSD_kpc__Tg[iT].data
            D['/masked/mask/ALL_tauV__Tg/%d' % iT] = ALL_tauV__Tg[iT].mask
            D['/masked/mask/ALL_SFR__Tg/%d' % iT] = ALL_SFR__Tg[iT].mask
            D['/masked/mask/ALL_SFRSD__Tg/%d' % iT] = ALL_SFRSD__Tg[iT].mask
            D['/masked/mask/ALL_SFRSD_kpc__Tg/%d' % iT] = ALL_SFRSD_kpc__Tg[iT].mask

        for k in D.keys():
            try:
                file.create_dataset(k, data = D[k], compression='gzip', compression_opts=4)
            except TypeError:
                file.create_dataset(k, data = D[k])
    
    if plot:
        NRows = 4
        NCols = 3
          
        f, axArr = plt.subplots(NRows, NCols)
        f.set_dpi(96)
        f.set_size_inches(10,12)
           
        for ax in f.axes:
            ax.set_axis_off()
      
        xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
       
        NAxes = len(f.axes) 
               
        k = 0
          
        for i in range(0, NRows):
            for j in range(0, NCols):
                ax = axArr[i, j]
     
                if k < len(tSF_to_plot):
                    iT = tSF_to_plot[k]
                  
                if i < (NRows - 1) or j < (NCols - 1):
                    ax.set_axis_on()
                    x = np.ma.log10(ALL_SFR__Tg[iT])
                    y = np.ma.log10(ALL_SFR_Ha__g)
                    mask = x.mask | y.mask
                    xm = x[~mask]
                    ym = y[~mask]
                    age = tSF__T[iT]
                    print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age/1e6, mask.sum(), len(x), len(x) - mask.sum())
                    xran = [-6,0]
                    yran = [-6,0]
                    scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.3)
                    binsx = np.linspace(-6.,0., 101)
                    binsy = np.linspace(min(ym),max(ym), 101)
                    density_contour(xm, ym, binsx, binsy, ax=ax)
                    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
                    txt = 'age: %.2d Myr - $R_S$: %.2f' % ((age / 1e6), correl_SFR[iT])
                    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
                    ax.text(0.05, 0.92, txt, fontsize = 12, 
                            transform = ax.transAxes, 
                            verticalalignment = 'top', horizontalalignment= 'left', 
                            bbox = textbox)
                    ax.grid()
                    ax.set_xlim(xran)
                    ax.set_ylim(yran)
                    if j == 0 and i == 1:
                        ax.set_ylabel(ylabel)
                    if j == 1 and i == (NRows - 1):
                        ax.set_xlabel(xlabel)
                    k += 1
                elif i == (NRows - 1) and j == (NCols - 1):
                    ax.set_axis_on()
                    ax.plot(np.log10(tSF__T), correl_SFR, 'k-', label = r'$R_S$')
                    ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
                    ax.legend(fontsize=12, frameon=False)
                    ax.xaxis.set_major_locator(MultipleLocator(1))
                    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.set_ylim([0., 1.])
                    ax.grid(which = 'minor')
       
        f.savefig('SFR_allages.png')
        plt.close(f)
      
        NCols = 3 
        NRows = 4
 
        pos_y_ini = 0.38
        pos_step = 0.09
        Rfontsize = 10
          
        f, axArr = plt.subplots(NRows, NCols)
        f.set_dpi(96)
        f.set_size_inches(10,12)
           
        for ax in f.axes:
            ax.set_axis_off()
       
        xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
           
        NAxes = len(f.axes) 
        k = 0
               
        for i in range(0, NRows):
            for j in range(0, NCols):
                ax = axArr[i, j]
 
                if k < len(tSF_to_plot):
                    iT = tSF_to_plot[k]
                 
                if i < (NRows - 1) or j < (NCols - 1):
                    ax.set_axis_on()
                    age = tSF__T[iT]
                    n_mask = n_tot = 0
                 
                    for iR, RUp in enumerate(RRange):
                        if iR == 0:
                            RMask = RbinCenter__r <= RUp
                            txt = 'R <= %.1f HLR' % RUp
                        else:
                            RDown =  RRange[iR - 1]
                            RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                            txt = '%.1f < R <= %.1f HLR' % (RDown, RUp)
                             
                        c = RColor[iR] 
                        x = np.ma.log10(ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                        y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten())
                        mask = x.mask | y.mask
                        xm = x[~mask]
                        ym = y[~mask]
                        n_mask += mask.sum()
                        n_tot += len(x)
                         
                        if i == 0 and j == 0:
                            pos_y = pos_y_ini - (iR * pos_step)
                            textbox = dict(alpha = 0.)
                            ax.text(0.05, pos_y, txt,
                                    fontsize = Rfontsize, color = c, 
                                    transform = ax.transAxes, 
                                    va = 'top', ha = 'left', 
                                    bbox = textbox)
                        scat = ax.scatter(xm, ym, c = c, marker = 'o', s = 1., edgecolor = 'none', alpha = 1.)
                    print 'SigmaSFR x SigmaSFR_Ha Age: %.2f Myr: masked %d points of %d (Total: %d)' % (age/1e6, n_mask, n_tot, n_tot - n_mask)
                             
                    ax.legend(loc = 'lower left', fontsize=12, frameon=False)
                    age = tSF__T[iT]
                    xran = [-3.5, 1.]
                    yran = [-3.5, 1.]
                    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                    # binsx = np.linspace(-4.5, 1., 51)
                    # binsy = np.linspace(min(ym),max(ym), 51)
                    # density_contour(xm, ym, binsx, binsy, ax=ax)
                    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
                    txt = 'age: %.2d Myr - $R_S$: %.2f' % ((age / 1e6), correl_SFRSD__rT[0, iT])
                    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
                    ax.text(0.05, 0.92, txt, fontsize = 12, 
                            transform = ax.transAxes, 
                            verticalalignment = 'top', horizontalalignment= 'left', 
                            bbox = textbox)
                    ax.grid()
                    ax.set_xlim(xran)
                    ax.set_ylim(yran)
                    ax.xaxis.set_major_locator(MultipleLocator(1))
                    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
                    ax.yaxis.set_major_locator(MultipleLocator(1))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
                 
                    if j == 0 and i == 1:
                        ax.set_ylabel(ylabel)
                    if j == 1 and i == (NRows - 1):
                        ax.set_xlabel(xlabel)
                    k += 1
                elif i == (NRows - 1) and j == (NCols - 1):
                    ax.set_axis_on()
                    ax.plot(np.log10(tSF__T), correl_SFRSD__rT[0, :], 'k-', label = r'$R_S$')
 
                    for iR, RUp in enumerate(RRange):
                        if iR == 0:
                            RMask = RbinCenter__r <= RUp
                        else:
                            RDown =  RRange[iR - 1]
                            RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                        iiR = iR + 1
                         
                        c = RColor[iR]
                        ax.plot(np.log10(tSF__T), correl_SFRSD__rT[iiR, :], color = c, ls = '--', label = None)
 
                    ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
                    ax.legend(fontsize=12, frameon=False)
                    ax.xaxis.set_major_locator(MultipleLocator(1))
                    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
                    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
                    ax.set_ylim([0., 1.])
                    ax.grid(which = 'minor')
       
        #f.tight_layout()   
        f.savefig('SigmaSFR_allages.pdf')
        plt.close(f)
          
        f = plt.figure()
        ax = f.gca()
        ax.plot(np.log10(tSF__T), correl_SFRSD__rT[0, :], 'k-', label = r'$R_S(\Sigma_{\mathrm{SFR}})$')
        ax.plot(np.log10(tSF__T), correl_SFRSD_kpc__rT[0, :], 'k-', label = r'$R_S(\Sigma_{\mathrm{SFR}}\ \star)$')
        ax.plot(np.log10(tSF__T), correl_SFR, 'k--', label = r'$R_S(\mathrm{SFR})$')
        ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
        ax.set_ylabel(r'$R_S$')
        ax.set_ylim([0, 1.])
        ax.grid()
        ax.legend(loc = 'best', fontsize = 12)
        f.tight_layout()
        f.savefig('Rs_SFR.pdf')
        plt.close(f)
  
        for iT,tSF in enumerate(tSF__T):
            x1 = np.ma.log10(ALL_SFR__Tg[iT])
            x1label = r'$\log\ \mathrm{SFR}_\star\ [M_\odot yr^{-1}]$' 
            y1 = np.ma.log10(ALL_SFR_Ha__g)
            y1label = r'$\log\ \mathrm{SFR}_{neb}\ [M_\odot yr^{-1}]$' 
            x2 = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
            x2label = r'$\log\ \Sigma_{\mathrm{SFR}}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            y2 = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten())
            y2label = r'$\log\ \Sigma_{\mathrm{SFR}}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            fname = 'SFReSFRSD_%sMyr.png' % str(tSF / 1.e6)
            mask1 = ~(x1.mask | y1.mask)
            x1m = x1[mask1]
            y1m = y1[mask1]
            mask2 = ~(x2.mask | y2.mask)
            x2m = x2[mask2]
            y2m = y2[mask2]
            plotCid(x1m, x1label, y1m, y1label, x2m, x2label, y2m, y2label, fname)

            x = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg[:, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
            xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
            ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
            xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} pc^{-2}]$' 
            fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr1.png' % str(tSF / 1.e6)
            plotSFR(xm,ym,xlabel,ylabel,xlim,ylim,tSF,fname)

            x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[:, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
            xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
            ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
            xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr2.png' % str(tSF / 1.e6)
            plotSFR(xm,ym,xlabel,ylabel,xlim,ylim,tSF,fname)

            x = np.ma.log10(ALL_SFR__Tg[iT])
            y = np.ma.log10(ALL_SFR_Ha__g)
            xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
            ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
            fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
            xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
            ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
            plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)

            x = ALL_tauV__Trg[iT, :, :].flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten() / ALL_aSFRSD__Trg[iT, :, :].flatten())
            xlabel = r'$\tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr1.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 

            x = ALL_tauV__Trg[iT, :, :].flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
            xlabel = r'$\tau_V^\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr2.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 

            x = ALL_tau_V_neb__rg.flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten() / ALL_aSFRSD__Trg[iT, :, :].flatten())
            xlabel = r'$\tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr1.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 

            x = ALL_tau_V_neb__rg.flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
            xlabel = r'$\tau_V^{neb}(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr2.png' % str(tSF / 1.e6)
            plotTau(x,y,xlabel,ylabel,None,None,tSF,fname)
                         
            ###################### RADIUS COLOR ######################
            xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
            ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
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
                        
                x = np.ma.log10(ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten())
                mask = ~(x.mask | y.mask)
                xm = x[mask]
                ym = y[mask]
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
             
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
 
            x = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
             
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % ((ym/xm).mean(), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
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
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                mask = ~(x.mask | y.mask)
                xm = x[mask]
                ym = y[mask]
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
 
            x = ALL_tauV__Trg[iT, :, :].flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
 
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
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
                 
                x = ALL_tau_V_neb__rg[RMask, :].flatten()
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                mask = ~(x.mask | y.mask)
                xm = x[mask]
                ym = y[mask]
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
                 
            x = ALL_tau_V_neb__rg.flatten()
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
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
 
                x = np.ma.log10(ALL_tauV__Trg[iT, RMask, :].flatten())
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                mask = ~(x.mask | y.mask)
                xm = x[mask]
                ym = y[mask]
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
 
            x = np.ma.log10(ALL_tauV__Trg[iT, :, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
             
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
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
 
                x = np.ma.log10(ALL_tau_V_neb__rg[RMask, :].flatten())
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                mask = ~(x.mask | y.mask)
                xm = x[mask]
                ym = y[mask]
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
 
            x = np.ma.log10(ALL_tau_V_neb__rg.flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
            mask = ~(x.mask | y.mask)
            xm = x[mask]
            ym = y[mask]
             
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
             
            xlabel = r'$\log\ \mathcal{M}_\star(R)$'
            ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
            fname = 'logMcorSD_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
            f = plt.figure(dpi = 96)
            f.set_size_inches(21.88,12.5)
            ax = f.gca()
            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown =  RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
 
                x = np.ma.log10(ALL_McorSD_GAL__rg[RMask, :].flatten())
                y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
                x.mask = y.mask
                xm = x
                ym = y
                scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
                plotRunningStatsAxis(ax, xm, ym, RColor[iR])
 
            x = np.ma.log10(ALL_McorSD_GAL__rg.flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
            x.mask = y.mask
            xm = x
            ym = y
 
            rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
            txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.10, 0.93, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
            ax.grid()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
            f.savefig(fname)
            plt.close(f)
