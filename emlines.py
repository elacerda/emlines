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

#plot = True
plot = False
#debug = False
debug = True

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size'] = 16
mpl.rcParams['font.family'] = 'sans-serif'
    
galaxiesListFile    = '/Users/lacerda/CALIFA/listOf300GalPrefixes.txt'
baseCode            = 'Bgsd6e'
#versionSuffix      = 'px1_q043.d14a'
versionSuffix       = 'v20_q043.d14a'
superFitsDir        = '/Users/lacerda/CALIFA/gal_fits/' + versionSuffix + '/'
#emLinesFitsDir      = '/Volumes/SAMSUNG/CALIFA/superfits/' + versionSuffix + '/'
emLinesFitsDir      = '/Users/lacerda/CALIFA/rgb-gas/' + versionSuffix + '/'

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
    #listOfPrefixes = listOfPrefixes[0:40]        # Q&D tests ...
    listOfPrefixes = ['K0073\n']
    
N_gals = len(listOfPrefixes)

# Setup bins for Radial profiles (in units of HLR)
RbinIni , RbinFin , RbinStep = 0.0 , 2.0 , 0.1
Rbin__r         = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r   = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins          = len(RbinCenter__r)

# Decide whether or not to weight radial profiles averages
weiRadProf = False

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     # ageMax array to compute alogZ stats for st-pops younger than ageMax__U.
#     # To be used in conjunction with flag_t when calling the compute-alogZ-stuff routine
#     ageMax__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
#     ageMax__U = np.array([1.0 , 2.0 , 14.2]) * 1.e9     # Q&D test
# 
#     N_U = len(ageMax__U)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

# SFR-time-scale array (index __T)
tSF__T = np.array([10.01 , 25.2 , 63.2, 100.1 , 158.6 ,199.6 ,1400.2 ]) * 1.e6
N_T = len(tSF__T)

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
    
def tauVNeb_implot(tauVNeb__yx, err_tauVNeb__yx, f_obs_HaHb__yx, err_f_obs_HaHb__yx, fileName):
    f, axArr = plt.subplots(1, 3)
    f.set_size_inches(18, 4)
    
    ax      = axArr[0]
    ax.set_title(r'$\tau_V^{neb}$')
    sigma   = tauVNeb__yx.std()
    mean    = tauVNeb__yx.mean()
    vmax    = mean + 2. * sigma
    vmin    = mean - 2. * sigma
    im      = ax.imshow(tauVNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax      = axArr[1]
    ax.set_title(r'$\epsilon (\tau_V^{neb})$')
    sigma   = err_tauVNeb__yx.std()
    mean    = err_tauVNeb__yx.mean()
    vmax    = mean + 2. * sigma
    im      = ax.imshow(err_tauVNeb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax              = axArr[2]
    ax.set_title(r'$F_{H\beta}^{H\alpha} / \epsilon(F_{H\beta}^{H\alpha})$')
    signalToNoise   = np.abs(f_obs_HaHb__yx) / np.abs(err_f_obs_HaHb__yx) 
    sigma           = signalToNoise.std()
    mean            = signalToNoise.mean()
    vmax            = mean + 2. * sigma
    im              = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    f.savefig(fileName)
    plt.close(f)

def Lint_Ha_implot(Lint_Ha__yx, err_Lint_Ha__yx, fileName):
    f, axArr = plt.subplots(1, 3)
    f.set_size_inches(18, 4)
    
    ax      = axArr[0]
    ax.set_title(r'$L_{H_\alpha}^{int}$')
    sigma   = Lint_Ha__yx.std()
    mean    = Lint_Ha__yx.mean()
    vmax    = mean + sigma
    vmin    = mean - sigma
    im      = ax.imshow(Lint_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax      = axArr[1]
    ax.set_title(r'$\epsilon (L_{H_\alpha}^{int})$')
    sigma   = err_Lint_Ha__yx.std()
    mean    = err_Lint_Ha__yx.mean()
    vmax    = mean + sigma
    im      = ax.imshow(err_Lint_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax              = axArr[2]
    ax.set_title(r'$L_{H_\alpha}^{int} / \epsilon(L_{H_\alpha}^{int})$')
    signalToNoise   = np.abs(Lint_Ha__yx) / np.abs(err_Lint_Ha__yx) 
    sigma           = signalToNoise.std()
    mean            = signalToNoise.mean()
    vmax            = mean + sigma
    im              = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    f.savefig(fileName)
    plt.close(f)

#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

def calc_SKlaw_Stuff(K, tSF, Rbin__r, maskOk=True, xOkMin = 0.05, tauVOkMin = 0.05):
    '''
    Computes radial profiles of a series of quantities relevant to study the SK-law.
    If maskOk == True then we add a spatial mask to the images so that spaxels where
    (a) there is too little light in the populations of age < tSF (ie, those used to defined a SFR), and
    **AND** (b) too little dust. Otherwise all spaxels are considered in the calculations.
    The motivation is to ensure we have enough young stars to make the SFR credible, and similarly for
    the dust optical depth (tauV).

    ==> return aSFRSD__r, atauV__r, atauV2mu__r, aMcorSD__r ,\
           alogSFRSD__r, alogtauV__r, alogtauV2mu__r , alogMcorSD__r , \
           alogZ_mass__r , alogZ_flux__r , at_flux__r

    Cid@Lagoa - 23/Jun/2014
    '''


    #--------------------------------------------------------------------------
    # Define mask to pick only populations younger than the input tSF in the computation of SFR & SFRSD.
    flag__t = K.ageBase <= tSF

    # SRFSD "raw" image
    # Note that we are NOT dezonifying SFR__z, since it will be compared to the un-dezonifiable tauV!
    SFR__z    = K.Mini__tZz[flag__t,:,:].sum(axis=1).sum(axis=0) / tSF
    SFRSD__yx = K.zoneToYX(SFR__z / K.zoneArea_pc2, extensive=False, surface_density=False)

    # Dust optical depth (at V band) and tauV/mu Dust/Stars ratio "raw" images
    tauV__yx    = K.A_V__yx / (2.5 * np.log10(np.exp(1.)))
    tauV2mu__yx = tauV__yx / K.McorSD__yx

    # tauV2SFR Dust/SFR ratio "raw" images - to produce a variable comparable to tauV2mu
    tauV2SFR__yx = tauV__yx / SFRSD__yx

    # Compute alogZ_****_z & define alogZ_****__yx "raw" images. OBS: These are all-ages-Z's
    alogZ_mass__z , alogZ_flux__z = calc_alogZ__z(K)
    alogZ_mass__yx = K.zoneToYX(alogZ_mass__z,extensive=False)
    alogZ_flux__yx = K.zoneToYX(alogZ_flux__z,extensive=False)

    # flux weighted mean age & McorSD = both already pre-computed by PyCASSO, but we my want to mask'em below
    at_flux__yx = K.at_flux__yx.copy()
    McorSD__yx = K.McorSD__yx.copy()
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Mask (or not!) points with too little (< xOkMin) light in the flag__t age-range AND too little extinction.
    if maskOk:

        # Compute xOk "raw" image
        x__tZz  =  K.popx / K.popx.sum(axis=1).sum(axis=0)
        xOk__z  = x__tZz[flag__t,:,:].sum(axis=1).sum(axis=0)
        xOk__yx = K.zoneToYX(xOk__z, extensive=False, surface_density=False)

        # Build mask to keep only pts with xOk__yx >= xOkMin  AND  tauV__yx >= tauVOkMin
        maskNotOk__yx = (xOk__yx < xOkMin) | (tauV__yx < tauVOkMin)

        # Mask raw images
        tauV__yx[maskNotOk__yx]       = np.ma.masked
        SFRSD__yx[maskNotOk__yx]      = np.ma.masked
        tauV2mu__yx[maskNotOk__yx]    = np.ma.masked
        alogZ_mass__yx[maskNotOk__yx] = np.ma.masked
        alogZ_flux__yx[maskNotOk__yx] = np.ma.masked
        at_flux__yx[maskNotOk__yx]    = np.ma.masked
        McorSD__yx[maskNotOk__yx]     = np.ma.masked
        tauV2SFR__yx [maskNotOk__yx]     = np.ma.masked
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Radial profiles of SFRSD, McorSD, tauV & tauV/mu
    SFRSD__r    = K.radialProfile(SFRSD__yx, Rbin__r, rad_scale=K.HLR_pix)
    tauV__r     = K.radialProfile(tauV__yx, Rbin__r, rad_scale=K.HLR_pix)
    tauV2mu__r  = K.radialProfile(tauV2mu__yx, Rbin__r, rad_scale=K.HLR_pix)
    McorSD__r   = K.radialProfile(McorSD__yx, Rbin__r, rad_scale=K.HLR_pix)
    tauV2SFR__r = K.radialProfile(tauV2SFR__yx, Rbin__r, rad_scale=K.HLR_pix)

    # alog*__r versions! (ie, radial-mean of the log(stuff)). ATT: log(0) warnings are harmless
    logSFRSD__r    = K.radialProfile(np.log10(SFRSD__yx), Rbin__r, rad_scale=K.HLR_pix)
    logtauV__r     = K.radialProfile(np.log10(tauV__yx), Rbin__r, rad_scale=K.HLR_pix)
    logtauV2mu__r  = K.radialProfile(np.log10(tauV2mu__yx), Rbin__r, rad_scale=K.HLR_pix)
    logMcorSD__r   = K.radialProfile(np.log10(McorSD__yx), Rbin__r, rad_scale=K.HLR_pix)
    logtauV2SFR__r = K.radialProfile(np.log10(tauV2SFR__yx), Rbin__r, rad_scale=K.HLR_pix)

    # Radial profiles of alogZ's & at's: weighted or not cf the global weiRadProf boolean switch.
    if weiRadProf:
        Mcor__yx = K.zoneToYX(K.Mcor__z, extensive = True)
        Lobn__yx = K.zoneToYX(K.Lobn__z, extensive = True)
        if maskOk:
            Mcor__yx[maskNotOk__yx] = np.ma.masked
            Lobn__yx[maskNotOk__yx] = np.ma.masked
        alogZ_mass__r = radialProfileWeighted(alogZ_mass__yx, Mcor__yx, Rbin__r, K.HLR_pix, K.radialProfile)
        alogZ_flux__r = radialProfileWeighted(alogZ_flux__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
        at_flux__r    = radialProfileWeighted(at_flux__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
    else:
        alogZ_mass__r = K.radialProfile(alogZ_mass__yx, Rbin__r, rad_scale=K.HLR_pix)
        alogZ_flux__r = K.radialProfile(alogZ_flux__yx, Rbin__r, rad_scale=K.HLR_pix)
        at_flux__r    = K.radialProfile(at_flux__yx, Rbin__r, rad_scale=K.HLR_pix)
    #--------------------------------------------------------------------------

    # Pack thing in a dictionary ???
    return SFRSD__r, tauV__r, tauV2mu__r, McorSD__r ,\
           logSFRSD__r, logtauV__r, logtauV2mu__r , logMcorSD__r , \
           alogZ_mass__r , alogZ_flux__r , at_flux__r , \
           tauV2SFR__r , logtauV2SFR__r

def calc_alogZ__z(K, flag__t = None):
    '''
    Compute mass & flux weighted average log ("alog") Z for each zone z.
    Only st-pops satisfying the (optional) input flag__t (ageBase-related) mask are considered.
    Ex: If flag__t == True only for t < 3 Gyr, then we'l compute tha alogz_*__z for stars
    younger than 3 Gry. If flag__t is not provided then all ages are used.

    OBS: Using masked arrays to avoid division by zero when computing alogZ_****__z.

    ==> return alogZ_mass__z, alogZ_flux__z

    Cid@Lagoa - 23/Jun/2014
    '''

    # Define an mask-nothing all-ages mask if flag__t not provided on input
    if flag__t == None:
        flag__t = K.ageBase > 0

    # Define log of base metallicities **in solar units** for convenience
    logZBase__Z = np.log10(K.metBase / Zsun)

    # ==> alogZ_mass__z
    numerator__z   = np.tensordot( K.Mcor__tZz[flag__t,:,:] , logZBase__Z , (1,0) ).sum(axis=0)
    denominator__z = K.Mcor__tZz[flag__t,:,:].sum(axis=1).sum(axis=0)
    numerator__z   = np.ma.masked_array( numerator__z   , mask=denominator__z == 0)
    denominator__z = np.ma.masked_array( denominator__z , mask=denominator__z == 0)
    alogZ_mass__z  = numerator__z / denominator__z

    # ==> alogZ_flux__z
    numerator__z   = np.tensordot( K.Lobn__tZz[flag__t,:,:] , logZBase__Z , (1,0) ).sum(axis=0)
    denominator__z = K.Lobn__tZz[flag__t,:,:].sum(axis=1).sum(axis=0)
    numerator__z   = np.ma.masked_array( numerator__z   , mask=denominator__z == 0)
    denominator__z = np.ma.masked_array( denominator__z , mask=denominator__z == 0)
    alogZ_flux__z  = numerator__z / denominator__z

    return alogZ_mass__z, alogZ_flux__z

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

if __name__ == '__main__':
    #########################################################################        
    ################################## CID ##################################
    ######################################################################### 
    ALL_morfType_GAL__g     = np.ma.zeros((N_gals))
    ALL_Mcor_GAL__g         = np.ma.zeros((N_gals))
    ALL_at_flux_GAL__g      = np.ma.zeros((N_gals))
    ALL_McorSD_GAL__g       = np.ma.zeros((N_gals))
        # Radially dependent *__rg
    ALL_at_flux__rg         = np.ma.zeros((NRbins, N_gals))
    ALL_McorSD__rg          = np.ma.zeros((NRbins, N_gals))
        # SK-law *__Trg arrays
    ALL_SFRSD__Trg          = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_tauV__Trg           = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_tauV2mu__Trg        = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_McorSD__Trg         = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_logSFRSD__Trg       = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_logtauV__Trg        = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_logtauV2mu__Trg     = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_logMcorSD__Trg      = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_alogZ_mass__Trg     = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_alogZ_flux__Trg     = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_at_flux__Trg        = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_tauV2SFR__Trg       = np.ma.zeros((N_T,NRbins, N_gals))
    ALL_logtauV2SFR__Trg    = np.ma.zeros((N_T,NRbins, N_gals))
    #########################################################################        
    #########################################################################        

    ALL_tauVNeb__rg         = np.ma.zeros((NRbins, N_gals))
    ALL_tauVNeb2mu__rg      = np.ma.zeros((NRbins, N_gals))
    ALL_log_tauVNeb__rg     = np.ma.zeros((NRbins, N_gals))
    ALL_log_tauVNeb2mu__rg  = np.ma.zeros((NRbins, N_gals))
    ALL_logZNeb__rg         = np.ma.zeros((NRbins, N_gals))
    
    _ALL_tauV__g                = []
    _ALL_tauVNeb__g             = []
    _ALL_tauVNeb_mask__g        = []
    _ALL_logZNeb__g             = []
    _ALL_logZNeb_mask__g        = []
    _ALL_err_tauVNeb__g         = []
    _ALL_SFRSD_Ha__g            = []
    _ALL_tauVNeb2mu__g          = []
    _ALL_McorSD_GAL_zones__g    = []
    _ALL_morfType_GAL_zones__g  = []
    
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
            ALL_Mcor_GAL__g[iGal]               = np.ma.masked
            ALL_at_flux_GAL__g[iGal]            = np.ma.masked
            ALL_McorSD_GAL__g[iGal]             = np.ma.masked
            ALL_at_flux__rg[:, iGal]            = np.ma.masked
            ALL_McorSD__rg[:, iGal]             = np.ma.masked
            ALL_SFRSD__Trg[:, :, iGal]          = np.ma.masked
            ALL_tauV__Trg[:, :, iGal]           = np.ma.masked
            ALL_tauV2mu__Trg[:, :, iGal]        = np.ma.masked
            ALL_McorSD__Trg[:, :, iGal]         = np.ma.masked
            ALL_logSFRSD__Trg[:, :, iGal]       = np.ma.masked
            ALL_logtauV__Trg[:, :, iGal]        = np.ma.masked
            ALL_logtauV2mu__Trg[:, :, iGal]     = np.ma.masked
            ALL_logMcorSD__Trg[:, :, iGal]      = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal]     = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal]     = np.ma.masked
            ALL_at_flux__Trg[:, :, iGal]        = np.ma.masked
            ALL_tauV2SFR__Trg[:, :, iGal]       = np.ma.masked
            ALL_logtauV2SFR__Trg[:, :, iGal]    = np.ma.masked
            ALL_tauVNeb__rg[:, iGal]            = np.ma.masked
            ALL_tauVNeb2mu__rg[:, iGal]         = np.ma.masked
            ALL_log_tauVNeb__rg[:, iGal]        = np.ma.masked
            ALL_log_tauVNeb2mu__rg[:, iGal]     = np.ma.masked
            ALL_logZNeb__rg[:, iGal]            = np.ma.masked
            continue
        
        tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
        ALL_morfType_GAL__g[iGal] = tipo
        aux = np.ones_like(K.Mcor__z) * ALL_morfType_GAL__g[iGal] 
        _ALL_morfType_GAL_zones__g.append(aux)
        
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # morph = get_morph(califaName = galName)
        # 
        # if morph:
        #     ALL_morfType_GAL__g[iGal] = morph_number(morph['hubtyp'][0], morph['hubsubtyp'][0])
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

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
        ALL_McorSD__rg[:,iGal]  = K.radialProfile(K.McorSD__yx , Rbin__r, rad_scale=K.HLR_pix)

        # Compute & store galaxy-wide at_flux
        numerator__z   = K.Lobn__tZz.sum(axis=1).sum(axis=0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis=1).sum(axis=0)
        ALL_at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()

        # Store at_flux radial profile (weighted or not).
        if weiRadProf:
            Lobn__yx = K.zoneToYX(K.Lobn__z, extensive = True)
            ALL_at_flux__rg[:,iGal] = radialProfileWeighted(K.at_flux__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
        else:
            ALL_at_flux__rg[:,iGal] = K.radialProfile(K.at_flux__yx, Rbin__r, rad_scale=K.HLR_pix)

        # Compute SKlaw things for each tSF__T.
        for iT,tSF in enumerate(tSF__T):
            aux                             = calc_SKlaw_Stuff(K, tSF, Rbin__r, xOkMin=xOkMin, tauVOkMin=tauVOkMin)
            ALL_SFRSD__Trg[iT,:,iGal]       = aux[0]
            ALL_tauV__Trg[iT,:,iGal]        = aux[1]
            ALL_tauV2mu__Trg[iT,:,iGal]     = aux[2]
            ALL_McorSD__Trg[iT,:,iGal]      = aux[3]
            ALL_logSFRSD__Trg[iT,:,iGal]    = aux[4]
            ALL_logtauV__Trg[iT,:,iGal]     = aux[5]
            ALL_logtauV2mu__Trg[iT,:,iGal]  = aux[6]
            ALL_logMcorSD__Trg[iT,:,iGal]   = aux[7]
            ALL_alogZ_mass__Trg[iT,:,iGal]  = aux[8]
            ALL_alogZ_flux__Trg[iT,:,iGal]  = aux[9]
            ALL_at_flux__Trg[iT,:,iGal]     = aux[10]
            ALL_tauV2SFR__Trg[iT,:,iGal]    = aux[11]
            ALL_logtauV2SFR__Trg[iT,:,iGal] = aux[12]
        #########################################################################
        #########################################################################
        #########################################################################

        tauV__z    = K.A_V / (2.5 * np.log10(np.exp(1.)))
        _ALL_tauV__g.append(tauV__z)

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
        
        _ALL_tauVNeb__g.append(tauVNeb__z)
        _ALL_tauVNeb_mask__g.append(tauVNeb__z.mask)
        _ALL_err_tauVNeb__g.append(err_tauVNeb__z)
        
        tauVNeb2mu__yx      = tauVNeb__yx / K.McorSD__yx
        
        tauVNeb__r                          = K.radialProfile(tauVNeb__yx, Rbin__r, rad_scale = K.HLR_pix)
        tauVNeb2mu__r                       = K.radialProfile(tauVNeb2mu__yx, Rbin__r, rad_scale = K.HLR_pix)
        log_tauVNeb__r                      = K.radialProfile(np.log10(tauVNeb__yx), Rbin__r, rad_scale = K.HLR_pix)
        log_tauVNeb2mu__r                   = K.radialProfile(np.log10(tauVNeb2mu__yx), Rbin__r, rad_scale = K.HLR_pix)
        ALL_tauVNeb__rg[:, iGal]            = tauVNeb__r 
        ALL_tauVNeb2mu__rg[:, iGal]         = tauVNeb2mu__r 
        ALL_log_tauVNeb__rg[:, iGal]        = log_tauVNeb__r 
        ALL_log_tauVNeb2mu__rg[:, iGal]     = log_tauVNeb2mu__r
        
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
        
        SFRSD_Ha__z     = 2.05 * Lint_Ha__z / (1.e8 * K.zoneArea_pc2)
        
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

        logZNeb__r                  = K.radialProfile(logZNeb__yx, Rbin__r, rad_scale = K.HLR_pix)
        ALL_logZNeb__rg[:, iGal]    = logZNeb__r 
        
        _ALL_logZNeb__g.append(logZNeb__z)
        _ALL_logZNeb_mask__g.append(logZNeb__z.mask)
        
        if plot:
            tauVNeb_implot(tauVNeb__yx, err_tauVNeb__yx, f_obs_HaHb__yx, err_f_obs_HaHb__yx, K.califaID + '_' + versionSuffix + '-tauVNeb.png')
            Lint_Ha_implot(Lint_Ha__yx, err_Lint_Ha__yx, K.califaID + '_' + versionSuffix + '-Lint.png')
    
    aux                         = np.hstack(np.asarray(_ALL_tauVNeb__g))
    auxMask                     = np.hstack(np.asarray(_ALL_tauVNeb_mask__g))
    ALL_tauVNeb__g              = np.ma.masked_array(aux, mask = auxMask)
    ALL_err_tauVNeb__g          = np.ma.masked_array(np.hstack(np.asarray(_ALL_err_tauVNeb__g)), mask = auxMask)
    ALL_tauV__g                 = np.ma.masked_array(np.hstack(np.asarray(_ALL_tauV__g)), mask = auxMask)
    
    aux                         = np.hstack(np.asarray(_ALL_logZNeb__g))
    auxMask                     = np.hstack(np.asarray(_ALL_logZNeb_mask__g))
    ALL_logZNeb__g              = np.ma.masked_array(aux, mask = auxMask)
    ALL_SFRSD_Ha__g             = np.ma.masked_array(np.hstack(np.asarray(_ALL_SFRSD_Ha__g)), mask = auxMask)
    
    ALL_McorSD_GAL_zones__g     = np.ma.masked_array(np.hstack(np.asarray(_ALL_McorSD_GAL_zones__g)))
    ALL_morfType_GAL_zones__g   = np.ma.masked_array(np.hstack(np.asarray(_ALL_morfType_GAL_zones__g)))
    
    sys.exit('off fiiiii')

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

m = (ALL_err_tauVNeb__g<0.15) & (np.log10(ALL_tauVNeb__g) > -1) & (np.log10(ALL_tauV__g) > -1)
x = np.log10(ALL_tauVNeb__g[m])
y = np.log10(ALL_tauV__g[m])
z = ALL_logZNeb__g[m]
plt.clf()
plt.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)

nBox = 16
dxBox = (x.max() - x.min()) / (nBox - 1.)
aux = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
xbinCenter  = aux[0]
xMedian     = aux[1]
xMean       = aux[2]
xStd        = aux[3]
yMedian     = aux[4]
yMean       = aux[5]
yStd        = aux[6]
nInBin      = aux[7]

plt.plot(xMean, yMean, 'ob-', lw = 2)
plt.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'k')

plt.grid()
plt.ylabel(r'$\log\ \tau_V$')
plt.xlabel(r'$\log\ \tau_V^{neb}$')
cb = plt.colorbar()
cb.set_label(r'$\log\ Z_{neb}\ [Z_\odot]$')

##############################################################################################

m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0) & (ALL_morfType_GAL_zones__g > 8.5)
#m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g > 0)
c = 0
#c = np.log10(0.2)
logDGR = c + np.log10(ALL_tauVNeb__g[m]) - (1./1.4) * np.log10(ALL_SFRSD_Ha__g[m]) 
x = ALL_logZNeb__g[m]
y = logDGR
z = np.log10(ALL_McorSD_GAL_zones__g[m])
#z = ALL_morfType_GAL_zones__g[m]
plt.clf()
plt.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)

nBox = 16
dxBox = (x.max() - x.min()) / (nBox - 1.)
aux = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
xbinCenter  = aux[0]
xMedian     = aux[1]
xMean       = aux[2]
xStd        = aux[3]
yMedian     = aux[4]
yMean       = aux[5]
yStd        = aux[6]
nInBin      = aux[7]

plt.plot(xMean, yMean, 'ob-', lw = 2)
plt.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'k')

plt.grid()
plt.ylabel(r'$\log\ \delta_{DGR}$')
plt.xlabel(r'$\log\ Z_{neb}\ [Z_\odot]$')
cb = plt.colorbar()
cb.set_label(r'$\log\ \mu^{galaxy}_\star\ [M_\odot pc^{-2}]$')
#cb.set_label(r'Morphology type')

##############################################################################################

m = (ALL_err_tauVNeb__g<0.15) & (ALL_tauVNeb__g >= 0) & (ALL_tauV__g >= 0)
x = ALL_tauV__g[m]
y = ALL_tauVNeb__g[m]
z = 10 ** ALL_logZNeb__g[m]
plt.clf()
plt.scatter(x, y, c = z, cmap = 'jet', edgecolor = 'none', alpha = 0.5)

nBox = 16
dxBox = (x.max() - x.min()) / (nBox - 1.)
aux = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
xbinCenter  = aux[0]
xMedian     = aux[1]
xMean       = aux[2]
xStd        = aux[3]
yMedian     = aux[4]
yMean       = aux[5]
yStd        = aux[6]
nInBin      = aux[7]

plt.plot(xMean, yMean, 'ob-', lw = 2)
plt.plot(ALL_tauV__g[m], ALL_tauV__g[m], 'b--' )
#plt.plot(ALL_tauV__g[m], 1.2 * ALL_tauV__g[m]), 'k--')
plt.plot(ALL_tauV__g[m], 2. * ALL_tauV__g[m], 'g--')
plt.errorbar(xMean, yMean, yerr = yStd, xerr = xStd, c = 'k')

plt.grid()
plt.xlabel(r'$\tau_V$')
plt.ylabel(r'$\tau_V^{neb}$')
cb = plt.colorbar()
cb.set_label(r'$Z_{neb}\ [Z_\odot]$')

##############################################################################################