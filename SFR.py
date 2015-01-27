#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import numpy as np
from pycasso import fitsQ3DataCube
import matplotlib as mpl
from get_morfologia import get_morfologia
from scipy import stats as st
from lines import *
import os
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
import time
import sys
import argparse as ap

default = {
    'debug' : False,
    'underS06' : False,
    'hdf5' : None,
    'spiral' : False,
    'weiradprof' : True,
    'minpopx' : 0.05,
    'mintauv' : 0.05,
    'mintauvneb' : 0.05,
    'maxtauvneberr' : 999.,
    'rbinini' : 0.0,
    'rbinfin' : 2.0,
    'rbinstep' : 0.1,
}

def parser_args():
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--spiral', '-S',
                        action = 'store_true',
                        default = default['spiral'])
    parser.add_argument('--underS06', 
                        action = 'store_true',
                        default = default['underS06'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default['weiradprof'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
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
hdf5 = args.hdf5
onlySpiral = args.spiral
weiRadProf = args.weiradprof
mask_xOk = True
# Def smallest light fraction (in the flag__t-ageMax age-range) deemed to be Ok for our stats ...
xOkMin = args.minpopx
# Minimum tauV to be taken seriously ...
tauVOkMin = args.mintauv
tauVNebOkMin = args.mintauvneb
#tauVNebErrMax = 0.15
tauVNebErrMax = args.maxtauvneberr
RbinIni = args.rbinini
RbinFin = args.rbinfin
RbinStep = args.rbinstep

print 'debug: ', debug
print 'BPTLowS06: ', BPTLowS06
print 'hdf5: ', hdf5
print 'onlySpiral: ', onlySpiral
print 'weiRadProf: ', weiRadProf
print 'xOkMin: ', xOkMin
print 'tauVOkMin: ', tauVOkMin
print 'tauVNebOkMin: ', tauVNebOkMin
print 'tauVNebErrMax: ', tauVNebErrMax
print 'RbinIni: ', RbinIni
print 'RbinFin: ', RbinFin
print 'RbinStep: ', RbinStep

CALIFAWorkDir = '/Users/lacerda/CALIFA/'
galaxiesListFile = CALIFAWorkDir + 'listOf300GalPrefixes.txt'
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

f = open(galaxiesListFile, 'r')
listOfPrefixes = f.readlines()
f.close()

Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins = len(RbinCenter__r)

RColor = [ 'r', 'y', 'b', 'k' ]
RRange = [  .5, 1., 1.5, 2.  ]

if debug:
    listOfPrefixes = listOfPrefixes[1:20]        # Q&D tests ...
    listOfPrefixes = ['K0073\n']
    
N_gals = len(listOfPrefixes)

# SFR-time-scale array (index __T)
base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
tSF__T = base.ageBase
N_T = base.nAges

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
    alogZ_mass__z = np.ma.masked_array(numerator__z / denominator__z)

    # ==> alogZ_flux__z - ATT: There may be nan's here depending on flag__t!
    numerator__z = np.tensordot(K.Lobn__tZz[flag__t, :, :] , logZBase__Z , (1, 0)).sum(axis = 0)
    denominator__z = K.Lobn__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
    alogZ_flux__z = np.ma.masked_array(numerator__z / denominator__z)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Def galaxy-wide averages of alogZ in light & mass, but **discards** zones having
    # too little light fractions in the ages given by flag__t
    isOk__z = np.ones_like(K.Mcor__z, dtype = np.bool)
    
    # Define Ok flag: Zones with light fraction x < xOkMin are not reliable for alogZ (& etc) estimation!
    if xOkMin >= 0.:
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

    return alogZ_mass__z, alogZ_flux__z, alogZ_mass_GAL, alogZ_flux_GAL, isOkFrac_GAL , alogZ_mass__r, alogZ_flux__r
#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


if __name__ == '__main__':
    t_init_prog = time.clock()
    ALL_N_zones__g = np.ma.zeros((N_gals))

    #########################################################################        
    ################################## CID ##################################
    ######################################################################### 
    ALL_morfType_GAL__g = np.ma.zeros((N_gals))
    ALL_at_flux_GAL__g  = np.ma.zeros((N_gals))
    ALL_Mcor_GAL__g     = np.ma.zeros((N_gals))
    ALL_McorSD_GAL__g   = np.ma.zeros((N_gals))
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_califaID__rg            = np.ma.zeros((NRbins, N_gals), dtype = '|S5')
    ALL_morfType_GAL_zones__rg  = np.ma.zeros((NRbins, N_gals))
    ALL_Mr_GAL_zones__rg        = np.ma.zeros((NRbins, N_gals))
    ALL_ur_GAL_zones__rg        = np.ma.zeros((NRbins, N_gals))
    ALL_tau_V_neb__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD_Ha__rg           = np.ma.zeros((NRbins, N_gals))
    ALL_McorSD_GAL__rg          = np.ma.zeros((NRbins, N_gals))
    ALL_logZ_neb_S06__rg        = np.ma.zeros((NRbins, N_gals))
    
    ALL_califaID__Trg   = np.ma.zeros((N_T, NRbins, N_gals), dtype = '|S5')
    ALL_aSFRSD__Trg     = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_tau_V__Trg      = np.ma.zeros((N_T, NRbins, N_gals))
    #ALL_at_flux__Trg    = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_mass__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_flux__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    
    ALL_alogZ_mass_GAL__Tg  = np.ma.zeros((N_T, N_gals))
    ALL_alogZ_flux_GAL__Tg  = np.ma.zeros((N_T, N_gals))
    ALL_isOkFrac_GAL__Tg    = np.ma.zeros((N_T, N_gals))
    
    #ALL_integrated_at_flux__Tg      = np.ma.zeros((N_T, N_gals))
    ALL_integrated_SFR__Tg          = np.ma.zeros((N_T, N_gals))
    ALL_integrated_SFRSD__Tg        = np.ma.zeros((N_T, N_gals))
    
    ALL_integrated_SFR_Ha__g        = np.ma.zeros((N_gals))
    ALL_integrated_SFRSD_Ha__g      = np.ma.zeros((N_gals))
    ALL_integrated_L_int_Ha__g      = np.ma.zeros((N_gals))
    
    #########################################################################
    ############# temporary lists with shape N_Zone * N_gals ################
    #########################################################################
    _ALL_califaID_GAL_zones__g = []
    _ALL_morfType_GAL_zones__g = []
    _ALL_Mr_GAL_zones__g = []
    _ALL_ur_GAL_zones__g = []
    _ALL_Mcor__g = []
    _ALL_McorSD__g = []
    _ALL_Mcor_GAL_zones__g = []
    _ALL_McorSD_GAL_zones__g = []
    _ALL_at_flux_GAL_zones__g = []
    _ALL_tau_V_neb__g = []
    _ALL_tau_V_neb_err__g = []
    _ALL_tau_V_neb_mask__g = []
    _ALL_logZ_neb_S06__g = []
    _ALL_logZ_neb_S06_err__g = []
    _ALL_logZ_neb_S06_mask__g= []
    _ALL_SFR_Ha__g = []
    _ALL_SFRSD_Ha__g = []
    _ALL_F_obs_Ha__g = []
    _ALL_L_int_Ha__g = []
    _ALL_L_int_Ha_err__g = []
    _ALL_L_int_Ha_mask__g = []
    _ALL_dist_zone__g = []
    _ALL_x_young__g = []

    _ALL_alogZ_mass__g = []
    
    _ALL_EW_Ha__g = []
    _ALL_EW_Hb__g = []
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # _ALL_at_flux__g = []
    # _ALL_at_mass__g = []
    # _ALL_aZ_mass__g = []
    # _ALL_alogZ_mass__g = []
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    #########################################################################
    #########################################################################
    #########################################################################
    
    #########################################################################
    ####### lists with length N_T, each with shape (N_Zone, N_gals) #########
    #########################################################################
    _ALL_tau_V__Tg = []
    _ALL_tau_V_mask__Tg = []
    
    _ALL_SFR__Tg = []
    _ALL_SFR_mask__Tg = []
    _ALL_SFRSD__Tg = []
    _ALL_SFRSD_mask__Tg = []
    
    _ALL_alogZ_mass__Tg = []
    _ALL_alogZ_mass_mask__Tg = []
    _ALL_alogZ_flux__Tg = []
    _ALL_alogZ_flux_mask__Tg = []
    
    _ALL_at_flux__Tg = []
    _ALL_at_flux_mask__Tg = []

    for iT in range(N_T):
        _ALL_tau_V__Tg.append([])
        _ALL_tau_V_mask__Tg.append([])
        _ALL_SFR__Tg.append([])
        _ALL_SFR_mask__Tg.append([])
        _ALL_SFRSD__Tg.append([])
        _ALL_SFRSD_mask__Tg.append([])
        _ALL_alogZ_mass__Tg.append([])
        _ALL_alogZ_mass_mask__Tg.append([])
        _ALL_alogZ_flux__Tg.append([])
        _ALL_alogZ_flux_mask__Tg.append([])
        _ALL_at_flux__Tg.append([])
        _ALL_at_flux_mask__Tg.append([])
        
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_N_zones_tau_V = 0
    ALL_N_gals = 0
    ALL_N_zones = 0
        
    for iGal in np.arange(N_gals):
        t_init_gal = time.clock()
        galName = listOfPrefixes[iGal][:-1]

        CALIFASuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        CALIFAFitsFile = superFitsDir + galName + CALIFASuffix
        emLinesSuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
        
        # both files
        if not (os.path.isfile(CALIFAFitsFile) and os.path.isfile(emLinesFitsFile)):
            ALL_N_zones__g[iGal] = np.ma.masked
            ALL_morfType_GAL__g[iGal] = np.ma.masked
            ALL_at_flux_GAL__g[iGal] = np.ma.masked
            ALL_Mcor_GAL__g[iGal] = np.ma.masked
            ALL_McorSD_GAL__g[iGal] = np.ma.masked
            ALL_Mr_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_ur_GAL_zones__rg[:, iGal] = np.ma.masked 
            ALL_morfType_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_tau_V_neb__rg[:, iGal] = np.ma.masked
            ALL_logZ_neb_S06__rg[:, iGal] = np.ma.masked
            ALL_alogZ_mass_GAL__Tg[:, iGal] = np.ma.masked
            ALL_alogZ_flux_GAL__Tg[:, iGal] = np.ma.masked
            ALL_isOkFrac_GAL__Tg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD__Trg[:, :, iGal] = np.ma.masked
            ALL_tau_V__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal] = np.ma.masked
            ALL_McorSD_GAL__rg[:, iGal] = np.ma.masked
            ALL_integrated_SFR__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFR_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha__g[iGal] = np.ma.masked
            ALL_integrated_L_int_Ha__g[iGal] = np.ma.masked
            ALL_califaID__rg[:, iGal] = np.ma.masked
            ALL_califaID__Trg[:, :, iGal] = np.ma.masked

            print '<<< %s galaxy: miss files' % galName 
            continue
        
        K = fitsQ3DataCube(CALIFAFitsFile)
        tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
        ALL_morfType_GAL__g[iGal] = tipo
        
        # Only spiral
        if onlySpiral and tipo <= 8: 
            ALL_N_zones__g[iGal] = np.ma.masked
            ALL_morfType_GAL__g[iGal] = np.ma.masked
            ALL_at_flux_GAL__g[iGal] = np.ma.masked
            ALL_Mcor_GAL__g[iGal] = np.ma.masked
            ALL_McorSD_GAL__g[iGal] = np.ma.masked
            ALL_Mr_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_ur_GAL_zones__rg[:, iGal] = np.ma.masked 
            ALL_morfType_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_tau_V_neb__rg[:, iGal] = np.ma.masked
            ALL_logZ_neb_S06__rg[:, iGal] = np.ma.masked
            ALL_alogZ_mass_GAL__Tg[:, iGal] = np.ma.masked
            ALL_alogZ_flux_GAL__Tg[:, iGal] = np.ma.masked
            ALL_isOkFrac_GAL__Tg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD__Trg[:, :, iGal] = np.ma.masked
            ALL_tau_V__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal] = np.ma.masked
            ALL_McorSD_GAL__rg[:, iGal] = np.ma.masked
            ALL_integrated_SFR__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFR_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha__g[iGal] = np.ma.masked
            ALL_integrated_L_int_Ha__g[iGal] = np.ma.masked
            ALL_califaID__rg[:, iGal] = np.ma.masked
            ALL_califaID__Trg[:, :, iGal] = np.ma.masked

            print '<<< %s galaxy: is not a spiral (type: %d)' % (galName, tipo) 
            continue
        
        # read FITSFILE containing galaxy emission lines measured by R.G.B.
        # read_rgb_fits returns False if emLinesFitsFile does not exists.
        #read = read_rgb_fits(emLinesFitsFile, read_lines)
        K.loadEmLinesDataCube(emLinesFitsFile)
        
        # Problem in FITS file
        if K.EL.flux[0, :].sum() == 0.:
            ALL_N_zones__g[iGal] = np.ma.masked
            ALL_morfType_GAL__g[iGal] = np.ma.masked
            ALL_at_flux_GAL__g[iGal] = np.ma.masked
            ALL_Mcor_GAL__g[iGal] = np.ma.masked
            ALL_McorSD_GAL__g[iGal] = np.ma.masked
            ALL_morfType_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_tau_V_neb__rg[:, iGal] = np.ma.masked
            ALL_Mr_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_ur_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_logZ_neb_S06__rg[:, iGal] = np.ma.masked
            ALL_alogZ_mass_GAL__Tg[:, iGal] = np.ma.masked
            ALL_alogZ_flux_GAL__Tg[:, iGal] = np.ma.masked
            ALL_isOkFrac_GAL__Tg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD__Trg[:, :, iGal] = np.ma.masked
            ALL_tau_V__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal] = np.ma.masked
            ALL_McorSD_GAL__rg[:, iGal] = np.ma.masked
            ALL_integrated_SFR__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFR_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha__g[iGal] = np.ma.masked
            ALL_integrated_L_int_Ha__g[iGal] = np.ma.masked
            ALL_califaID__rg[:, iGal] = np.ma.masked
            ALL_califaID__Trg[:, :, iGal] = np.ma.masked
            
            print '<<< %s EmLines FITS problem' % galName
            continue

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
        
        maskLinesSNR__z = HbOk & O3Ok & HaOk & N2Ok
        maskFluxOk__z = (Ha >= 0) & (Hb >= 0) & (O3 >= 0) & (N2 >= 0)
        
        maskBPT__z = None
        
        if BPTLowS06:
            L = Lines()
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
            maskBPT__z = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        ##########################
        ##########################
        ##########################
                
        print '>>> Doing' , iGal , galName , 'hubtyp=', ALL_morfType_GAL__g[iGal], '|  Nzones=' , K.N_zone
        ALL_N_zones__g[iGal] = K.N_zone
        ALL_N_zones += K.N_zone
        ALL_N_gals += 1

        # a fake morfType per zone for all galaxies, creating a stamp for each zone
        aux = np.ones_like(K.Mcor__z) * ALL_morfType_GAL__g[iGal]
        _ALL_morfType_GAL_zones__g.append(aux)
        # The same above but for radial
        ALL_morfType_GAL_zones__rg[:, iGal] = np.ones_like(RbinCenter__r) * ALL_morfType_GAL__g[iGal]

        _ALL_Mr_GAL_zones__g.append(np.ones_like(K.Mcor__z) * np.float(K.masterListData['Mr']))
        ALL_Mr_GAL_zones__rg[:, iGal] = np.ones_like(RbinCenter__r) * np.float(K.masterListData['Mr'])
        
        _ALL_ur_GAL_zones__g.append(np.ones_like(K.Mcor__z) * np.float(K.masterListData['u-r']))
        ALL_ur_GAL_zones__rg[:, iGal] = np.ones_like(RbinCenter__r) * np.float(K.masterListData['u-r'])
        
        _ALL_califaID_GAL_zones__g.append(np.asarray([K.califaID for i in range(K.N_zone)]))
        
        zoneDistHLR = np.sqrt((K.zonePos['x'] - K.x0)**2. + (K.zonePos['y'] - K.y0)**2.) / K.HLR_pix
        _ALL_dist_zone__g.append(zoneDistHLR)
        
        ALL_califaID__rg[:,iGal] = np.asarray([ K.califaID for i in range(NRbins) ])
        ALL_califaID__Trg[...,iGal] = np.asarray([ K.califaID for i in range(N_T * NRbins) ]).reshape(N_T, NRbins)

        ##########################
        ####### STARLIGHT ########
        ##########################        
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        ALL_McorSD_GAL__g[iGal] = K.McorSD__yx.mean()

        # a fake McorSD per zone for all galaxies, creating a stamp for each zone
        aux = np.ones_like(K.Mcor__z) * ALL_McorSD_GAL__g[iGal]
        _ALL_McorSD_GAL_zones__g.append(aux)
        _ALL_Mcor__g.append(K.Mcor__z)
        _ALL_McorSD__g.append(K.Mcor__z / K.zoneArea_pc2)
        ALL_Mcor_GAL__g[iGal] = K.Mcor_tot.sum()
        # McorSD by radius
        ALL_McorSD_GAL__rg[:, iGal] = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)

        # a fake all galaxy mass in stars per zone for all galaxies, creating a stamp for each zone
        aux = np.ones_like(K.Mcor__z) * ALL_Mcor_GAL__g[iGal]
        _ALL_Mcor_GAL_zones__g.append(aux)

        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        ALL_at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()
        
        aux = np.ones_like(K.Mcor__z) * ALL_at_flux_GAL__g[iGal] 
        _ALL_at_flux_GAL_zones__g.append(aux)
        
        # pop_young
        flag_young__t = K.ageBase <= 1e8
        aux = K.popx[flag_young__t, :, :].sum(axis = 1).sum(axis = 0)
        aux /= K.popx.sum(axis = 1).sum(axis = 0)
        x_young__z = np.ma.masked_array(aux)
        _ALL_x_young__g.append(x_young__z)
        
        # alogZ_mass per zone
        _ALL_alogZ_mass__g.append(K.alogZ_mass__z)
        
        for iT, tSF in enumerate(tSF__T):
            #--------------------------------------------------------------------------
            # Define mask to pick only populations younger than the input tSF in the computation of SFR & SFRSD.
            flag__t = K.ageBase <= tSF
            # SRFSD "raw" image
            # Note that we are NOT dezonifying SFR__z, since it will be compared to the un-dezonifiable tauV!
            aux = K.Mini__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0) / tSF
            SFR__z = np.ma.masked_array(aux)
            SFRSD__z = SFR__z / K.zoneArea_pc2
                                        
            tau_V__z = np.ma.masked_array(K.tau_V__z)
                            
            if xOkMin >= 0.:
                # Compute xOk "raw" image
                x__tZz = K.popx / K.popx.sum(axis = 1).sum(axis = 0)
                xOk__z = x__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
                
                maskNotOk__z = (xOk__z < xOkMin) | (tau_V__z < tauVOkMin) 
                 
                tau_V__z[maskNotOk__z] = np.ma.masked
                SFR__z[maskNotOk__z] = np.ma.masked
                SFRSD__z[maskNotOk__z] = np.ma.masked
                
            integrated_SFR = SFR__z.sum()
            integrated_SFRSD = integrated_SFR / K.zoneArea_pc2.sum()
                        
            _ALL_tau_V__Tg[iT].append(tau_V__z.data)
            _ALL_tau_V_mask__Tg[iT].append(tau_V__z.mask)
                            
            ALL_integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL_integrated_SFRSD__Tg[iT, iGal] = integrated_SFRSD
            aSFRSD__yx = K.zoneToYX(SFR__z, extensive = True)
            aSFRSD__r = K.radialProfile(aSFRSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL_aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
            tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
            ALL_tau_V__Trg[iT, :, iGal] = K.radialProfile(tau_V__yx, Rbin__r, rad_scale = K.HLR_pix)
            _ALL_SFR__Tg[iT].append(SFR__z.data)
            _ALL_SFR_mask__Tg[iT].append(SFR__z.mask)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z.data)
            _ALL_SFRSD_mask__Tg[iT].append(SFRSD__z.mask)
            
        for iT, tSF in enumerate(tSF__T):
            #--------------------------------------------------------------------------
            # Define mask to pick only populations younger than the input tSF in the computation of SFR & SFRSD.
            flag__t = K.ageBase <= tSF

            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # #at_flux(t)
            # #at_flux__z = np.ma.masked_array(K.at_flux__z)
            # popx_sumZ = K.popx[flag__t, :, :].sum(axis=1)
            # popx_sum = popx_sumZ.sum(axis=0)
            # popx_sumZ /= popx_sum
            # aux = np.tensordot(popx_sumZ, np.log10(K.ageBase[flag__t]), (0, 0))
            # at_flux__z = np.ma.masked_array(aux)
            # at_flux__z[np.isnan(aux)] = np.ma.masked
            # _ALL_at_flux__Tg[iT].append(at_flux__z.data)
            # _ALL_at_flux_mask__Tg[iT].append(at_flux__z.mask)
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

            aux = calc_alogZ_Stuff(K, flag__t, Rbin__r, weiRadProf, xOkMin = xOkMin)
            _ALL_alogZ_mass__Tg[iT].append(aux[0])
            _ALL_alogZ_mass_mask__Tg[iT].append(aux[0].mask)
            _ALL_alogZ_flux__Tg[iT].append(aux[1])
            _ALL_alogZ_flux_mask__Tg[iT].append(aux[1].mask)
            ALL_alogZ_mass_GAL__Tg[iT, iGal] = aux[2]
            ALL_alogZ_flux_GAL__Tg[iT, iGal] = aux[3]
            ALL_isOkFrac_GAL__Tg[iT, iGal] = aux[4]
            ALL_alogZ_mass__Trg[iT, :, iGal] = aux[5]
            ALL_alogZ_flux__Trg[iT, :, iGal] = aux[6]
                
        ##########################
        ##########################
        ##########################    

        ##########################
        ########## tau_V #########
        ##########################
        mask_tmp = maskFluxOk__z & maskLinesSNR__z
        
        K.EL._forceMask = ~mask_tmp #Changing global EL mask

        if BPTLowS06:
            K.EL._forceMask = ~(mask_tmp & maskBPT__z) #Changing global EL mask
        
        maskOkTauVNeb = (K.EL.tau_V_neb__z >= tauVNebOkMin) & (K.EL.tau_V_neb_err__z <= tauVNebErrMax)
        
        K.EL._forceMask = ~(mask_tmp & maskOkTauVNeb) #Changing global EL mask 
        
        N_zones_tau_V = len(K.EL.tau_V_neb__z.compressed())
        print 'tauV calculated for %d zones (maskOK and maskOkTauVNeb)' % N_zones_tau_V
        ALL_N_zones_tau_V += N_zones_tau_V
        
        tau_V_neb__z = K.EL.tau_V_neb__z
        tau_V_neb_err__z = K.EL.tau_V_neb_err__z
        tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
        tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
        tau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)

        ALL_tau_V_neb__rg[:, iGal] = tau_V_neb__r
        _ALL_tau_V_neb__g.append(tau_V_neb__z.data)
        _ALL_tau_V_neb_err__g.append(tau_V_neb_err__z.data)
        _ALL_tau_V_neb_mask__g.append(tau_V_neb__z.mask)
        ##########################
        ##########################
        ##########################

        ##########################
        ######### Z_neb ##########
        ##########################
        logZ_neb_S06__z = K.EL.logZ_neb_S06__z
        logZ_neb_S06_err__z = K.EL.logZ_neb_S06_err__z
        logZ_neb_S06__yx = K.zoneToYX(K.EL.logZ_neb_S06__z, extensive = False)
        logZ_neb_S06_err__yx = K.zoneToYX(K.EL.logZ_neb_S06_err__z, extensive = False)
        logZ_neb_S06__r = K.radialProfile(logZ_neb_S06__yx, Rbin__r, rad_scale = K.HLR_pix)

        ALL_logZ_neb_S06__rg[:, iGal] = logZ_neb_S06__r
        _ALL_logZ_neb_S06__g.append(logZ_neb_S06__z.data)
        _ALL_logZ_neb_S06_err__g.append(logZ_neb_S06_err__z.data)
        _ALL_logZ_neb_S06_mask__g.append(K.EL.logZ_neb_S06__z.mask)
        ##########################
        ##########################
        ##########################

        ##########################
        ########### EW ###########
        ##########################
        _ALL_EW_Ha__g.append(K.EL.EW[i_Ha, :])
        _ALL_EW_Hb__g.append(K.EL.EW[i_Hb, :])
        ##########################
        ##########################
        ##########################

        ##########################
        #### intrinsic Ha Lum ####
        ##########################
        _ALL_F_obs_Ha__g.append(K.EL.flux[i_Ha, :])
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
        
        _ALL_L_int_Ha__g.append(L_int_Ha__z.data)
        _ALL_L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        _ALL_L_int_Ha_mask__g.append(L_int_Ha__z.mask)
        ALL_integrated_L_int_Ha__g[iGal] = integrated_L_int_Ha
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

        ALL_aSFRSD_Ha__rg[:, iGal] = aSFRSD_Ha__r
        
        ALL_integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha
        ALL_integrated_SFRSD_Ha__g[iGal] = integrated_SFRSD_Ha
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        ##########################
        ##########################
        ##########################
        
        print 'time per galaxy: %s %.2f' % (galName, time.clock() - t_init_gal)
    
    print 'Total of %d galaxies (%d zones): %s zones for tau_V' % (ALL_N_gals, ALL_N_zones, ALL_N_zones_tau_V)
    
    ALL_dist_zone__g = np.ma.masked_array(np.hstack(np.asarray(_ALL_dist_zone__g)))

    aux = np.hstack(_ALL_tau_V_neb__g)
    auxMask = np.hstack(_ALL_tau_V_neb_mask__g)
    ALL_tau_V_neb__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_tau_V_neb_err__g = np.ma.masked_array(np.hstack(_ALL_tau_V_neb_err__g), mask = auxMask)

    aux = np.hstack(_ALL_logZ_neb_S06__g)
    auxMask = np.hstack(_ALL_logZ_neb_S06_mask__g)
    ALL_logZ_neb_S06__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_logZ_neb_S06_err__g = np.ma.masked_array(np.hstack(_ALL_logZ_neb_S06_err__g), mask = auxMask)

    aux = np.hstack(_ALL_L_int_Ha__g)
    auxMask = np.hstack(_ALL_L_int_Ha_mask__g)
    ALL_L_int_Ha__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_F_obs_Ha__g = np.ma.masked_array(np.hstack(_ALL_F_obs_Ha__g), mask = auxMask)
    ALL_SFR_Ha__g = np.ma.masked_array(np.hstack(_ALL_SFR_Ha__g), mask = auxMask)
    ALL_SFRSD_Ha__g = np.ma.masked_array(np.hstack(_ALL_SFRSD_Ha__g), mask = auxMask)
    
    ALL_Mcor__g = np.ma.masked_array(np.hstack(_ALL_Mcor__g))
    ALL_McorSD__g = np.ma.masked_array(np.hstack(_ALL_McorSD__g))
    ALL_Mcor_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_Mcor_GAL_zones__g))
    ALL_McorSD_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_McorSD_GAL_zones__g))
    ALL_morfType_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_morfType_GAL_zones__g))
    ALL_at_flux_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_at_flux_GAL_zones__g))

    ALL_x_young__g = np.ma.masked_array(np.hstack(_ALL_x_young__g))
    ALL_EW_Ha__g = np.ma.masked_array(np.hstack(_ALL_EW_Ha__g))
    ALL_EW_Hb__g = np.ma.masked_array(np.hstack(_ALL_EW_Hb__g))
    
    ALL_alogZ_mass__g = np.ma.masked_array(np.hstack(_ALL_alogZ_mass__g))

    ALL_Mr_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_Mr_GAL_zones__g))
    ALL_ur_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_ur_GAL_zones__g))
    ALL_califaID_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_califaID_GAL_zones__g))

    ALL_tau_V__Tg = []
    ALL_SFR__Tg = []
    ALL_SFRSD__Tg = []
    ALL_alogZ_mass__Tg = []
    ALL_alogZ_flux__Tg = []
    #ALL_at_flux__Tg = []
    
    correl_SFR__T = np.ones_like(tSF__T)
    correl_SFRSD__T = np.ones_like(tSF__T)
    correl_aSFRSD__rT = np.ones((1 + len(RRange), tSF__T.shape[0]))

    for iT, tSF in enumerate(tSF__T):
        aux = np.hstack(_ALL_tau_V__Tg[iT])
        auxMask = np.hstack(_ALL_tau_V_mask__Tg[iT])
        ALL_tau_V__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        
        aux = np.hstack(_ALL_SFR__Tg[iT])
        auxMask = np.hstack(_ALL_SFR_mask__Tg[iT])        
        ALL_SFR__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        
        aux = np.hstack(_ALL_SFRSD__Tg[iT])
        auxMask = np.hstack(_ALL_SFRSD_mask__Tg[iT])
        ALL_SFRSD__Tg.append(np.ma.masked_array(aux, mask = auxMask))

        aux = np.hstack(_ALL_alogZ_mass__Tg[iT])
        auxMask = np.hstack(_ALL_alogZ_mass_mask__Tg[iT])
        #ALL_aZ_mass__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_alogZ_mass__Tg.append(np.ma.masked_array(aux))

        aux = np.hstack(_ALL_alogZ_flux__Tg[iT])
        auxMask = np.hstack(_ALL_alogZ_flux_mask__Tg[iT])
        #ALL_aZ_flux__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_alogZ_flux__Tg.append(np.ma.masked_array(aux))

        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # aux = np.hstack(_ALL_at_flux__Tg[iT])
        # auxMask = np.hstack(_ALL_at_flux_mask__Tg[iT])
        # ALL_at_flux__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

        x = np.ma.log10(ALL_SFR__Tg[iT])
        y = np.ma.log10(ALL_SFR_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_SFR__T[iT] = st.spearmanr(xm, ym)[0]

        x = np.ma.log10(ALL_SFRSD__Tg[iT])
        y = np.ma.log10(ALL_SFRSD_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_SFRSD__T[iT] = st.spearmanr(xm, ym)[0]

        x = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
        y = np.ma.log10(ALL_aSFRSD_Ha__rg[:, :].flatten())
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        correl_aSFRSD__rT[0, iT] = st.spearmanr(xm, ym)[0]
        
        for iR, RUp in enumerate(RRange):
            if iR == 0:
                RMask = RbinCenter__r <= RUp
            else:
                RDown = RRange[iR - 1]
                RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                
            iiR = iR + 1

            x = np.ma.log10(ALL_aSFRSD__Trg[iT, RMask, :].flatten())
            y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten())
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            correl_aSFRSD__rT[iiR, iT] = st.spearmanr(xm, ym)[0]

    if hdf5 != None:
        import h5py
        
        filename = args.hdf5
            
        h5 = h5py.File(filename, 'w')
        
        D = {
            '/masked/data/morfType_GAL__g' : ALL_morfType_GAL__g.data,
            '/masked/data/at_flux_GAL__g' : ALL_at_flux_GAL__g.data,
            '/masked/data/Mcor_GAL__g' : ALL_Mcor_GAL__g.data,
            '/masked/data/McorSD_GAL__g' : ALL_McorSD_GAL__g.data,
            '/masked/data/morfType_GAL_zones__rg' : ALL_morfType_GAL_zones__rg.data,
            '/masked/data/Mr_GAL_zones__rg' : ALL_Mr_GAL_zones__rg.data,
            '/masked/data/ur_GAL_zones__rg' : ALL_ur_GAL_zones__rg.data,
            '/masked/data/tau_V_neb__rg' : ALL_tau_V_neb__rg.data,
            '/masked/data/alogZ_mass_GAL__Tg' : ALL_alogZ_mass_GAL__Tg.data,
            '/masked/data/alogZ_flux_GAL__Tg' : ALL_alogZ_flux_GAL__Tg.data,
            '/masked/data/isOkFrac_GAL__Tg' : ALL_isOkFrac_GAL__Tg.data,
            '/masked/data/aSFRSD_Ha__rg' : ALL_aSFRSD_Ha__rg.data,
            '/masked/data/aSFRSD__Trg' : ALL_aSFRSD__Trg.data,
            '/masked/data/logZ_neb_S06__rg': ALL_logZ_neb_S06__rg.data,
            '/masked/data/tau_V__Trg' : ALL_tau_V__Trg.data,
            '/masked/data/alogZ_mass__Trg' : ALL_alogZ_mass__Trg.data,
            '/masked/data/alogZ_flux__Trg' : ALL_alogZ_flux__Trg.data,
            '/masked/data/McorSD_GAL__rg' : ALL_McorSD_GAL__rg.data,
            '/masked/data/tau_V_neb__g' :  ALL_tau_V_neb__g.data,
            '/masked/data/tau_V_neb_err__g' : ALL_tau_V_neb_err__g.data,
            '/masked/data/L_int_Ha__g' : ALL_L_int_Ha__g.data,
            '/masked/data/F_obs_Ha__g' : ALL_F_obs_Ha__g.data,
            '/masked/data/SFR_Ha__g' : ALL_SFR_Ha__g.data,
            '/masked/data/SFRSD_Ha__g' : ALL_SFRSD_Ha__g.data,
            '/masked/data/Mcor__g' : ALL_Mcor__g.data,
            '/masked/data/McorSD__g' : ALL_McorSD__g.data,
            '/masked/data/Mcor_GAL_zones__g' : ALL_Mcor_GAL_zones__g.data,
            '/masked/data/Mr_GAL_zones__g' : ALL_Mr_GAL_zones__g.data,
            '/masked/data/ur_GAL_zones__g' : ALL_ur_GAL_zones__g.data,
            '/masked/data/califaID_GAL_zones__g' : ALL_califaID_GAL_zones__g.data,            
            '/masked/data/McorSD_GAL_zones__g' : ALL_McorSD_GAL_zones__g.data,
            '/masked/data/morfType_GAL_zones__g' : ALL_morfType_GAL_zones__g.data,
            '/masked/data/at_flux_GAL_zones__g' : ALL_at_flux_GAL_zones__g.data,
            '/masked/data/integrated_SFR__Tg' : ALL_integrated_SFR__Tg.data,
            '/masked/data/integrated_SFRSD__Tg' : ALL_integrated_SFRSD__Tg.data,
            '/masked/data/integrated_SFR_Ha__g' : ALL_integrated_SFR_Ha__g.data,
            '/masked/data/integrated_SFRSD_Ha__g' : ALL_integrated_SFRSD_Ha__g.data,
            '/masked/data/dist_zone__g' : ALL_dist_zone__g.data,
            '/masked/data/N_zones__g' : ALL_N_zones__g.data,
            '/masked/data/califaID__rg' : ALL_califaID__rg.data,
            '/masked/data/califaID__Trg' : ALL_califaID__Trg.data,
            '/masked/data/logZ_neb_S06__g' : ALL_logZ_neb_S06__g.data,
            '/masked/data/logZ_neb_S06_err__g' : ALL_logZ_neb_S06_err__g.data,
            '/masked/data/x_young__g' : ALL_x_young__g.data,
            '/masked/data/alogZ_mass__g': ALL_alogZ_mass__g.data,
            '/masked/data/EW_Ha__g': ALL_EW_Ha__g.data,
            '/masked/data/EW_Hb__g': ALL_EW_Hb__g.data, 
            '/masked/mask/morfType_GAL__g' : ALL_morfType_GAL__g.mask,
            '/masked/mask/at_flux_GAL__g' : ALL_at_flux_GAL__g.mask,
            '/masked/mask/Mcor_GAL__g' : ALL_Mcor_GAL__g.mask,
            '/masked/mask/McorSD_GAL__g' : ALL_McorSD_GAL__g.mask,
            '/masked/mask/morfType_GAL_zones__rg' : ALL_morfType_GAL_zones__rg.mask,
            '/masked/mask/Mr_GAL_zones__rg' : ALL_Mr_GAL_zones__rg.mask,
            '/masked/mask/ur_GAL_zones__rg' : ALL_ur_GAL_zones__rg.mask,
            '/masked/mask/tau_V_neb__rg' : ALL_tau_V_neb__rg.mask,
            '/masked/mask/alogZ_mass_GAL__Tg' : ALL_alogZ_mass_GAL__Tg.mask,
            '/masked/mask/alogZ_flux_GAL__Tg' : ALL_alogZ_flux_GAL__Tg.mask,
            '/masked/mask/isOkFrac_GAL__Tg' : ALL_isOkFrac_GAL__Tg.mask,
            '/masked/mask/aSFRSD_Ha__rg' : ALL_aSFRSD_Ha__rg.mask,
            '/masked/mask/aSFRSD__Trg' : ALL_aSFRSD__Trg.mask,
            '/masked/mask/logZ_neb_S06__rg': ALL_logZ_neb_S06__rg.mask,
            '/masked/mask/tau_V__Trg' : ALL_tau_V__Trg.mask,
            '/masked/mask/alogZ_mass__Trg' : ALL_alogZ_mass__Trg.mask,
            '/masked/mask/alogZ_flux__Trg' : ALL_alogZ_flux__Trg.mask,
            '/masked/mask/McorSD_GAL__rg' : ALL_McorSD_GAL__rg.mask,
            '/masked/mask/tau_V_neb__g' :  ALL_tau_V_neb__g.mask,
            '/masked/mask/tau_V_neb_err__g' : ALL_tau_V_neb_err__g.mask,
            '/masked/mask/L_int_Ha__g' : ALL_L_int_Ha__g.mask,
            '/masked/mask/F_obs_Ha__g' : ALL_F_obs_Ha__g.mask,
            '/masked/mask/SFR_Ha__g' : ALL_SFR_Ha__g.mask,
            '/masked/mask/SFRSD_Ha__g' : ALL_SFRSD_Ha__g.mask,
            '/masked/mask/Mcor__g' : ALL_Mcor__g.mask,
            '/masked/mask/McorSD__g' : ALL_McorSD__g.mask,
            '/masked/mask/Mr_GAL_zones__g' : ALL_Mr_GAL_zones__g.mask,
            '/masked/mask/ur_GAL_zones__g' : ALL_ur_GAL_zones__g.mask,
            '/masked/mask/califaID_GAL_zones__g' : ALL_califaID_GAL_zones__g.mask,
            '/masked/mask/Mcor_GAL_zones__g' : ALL_Mcor_GAL_zones__g.mask,
            '/masked/mask/McorSD_GAL_zones__g' : ALL_McorSD_GAL_zones__g.mask,
            '/masked/mask/morfType_GAL_zones__g' : ALL_morfType_GAL_zones__g.mask,
            '/masked/mask/at_flux_GAL_zones__g' : ALL_at_flux_GAL_zones__g.mask,
            '/masked/mask/integrated_SFR__Tg' : ALL_integrated_SFR__Tg.mask,
            '/masked/mask/integrated_SFRSD__Tg' : ALL_integrated_SFRSD__Tg.mask,
            '/masked/mask/integrated_SFR_Ha__g' : ALL_integrated_SFR_Ha__g.mask,
            '/masked/mask/integrated_SFRSD_Ha__g' : ALL_integrated_SFRSD_Ha__g.mask,
            '/masked/mask/dist_zone__g' : ALL_dist_zone__g.mask,
            '/masked/mask/N_zones__g' : ALL_N_zones__g.mask,
            '/masked/mask/califaID__rg' : ALL_califaID__rg.mask,
            '/masked/mask/califaID__Trg' : ALL_califaID__Trg.mask, 
            '/masked/mask/logZ_neb_S06__g' : ALL_logZ_neb_S06__g.mask,
            '/masked/mask/logZ_neb_S06_err__g' : ALL_logZ_neb_S06_err__g.mask,
            '/masked/mask/x_young__g' : ALL_x_young__g.mask, 
            '/masked/mask/alogZ_mass__g': ALL_alogZ_mass__g.mask,
            '/masked/mask/EW_Ha__g': ALL_EW_Ha__g.mask,
            '/masked/mask/EW_Hb__g': ALL_EW_Hb__g.mask, 
            '/data/correl_SFR__T' : correl_SFR__T,
            '/data/correl_SFRSD__T' : correl_SFRSD__T,
            '/data/correl_aSFRSD__rT' : correl_aSFRSD__rT,
            '/data/zones_tau_V' : ALL_N_zones_tau_V,
            '/data/gals' : ALL_N_gals,
            '/data/Nzones' : ALL_N_zones,
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

        for iT, tSF in enumerate(tSF__T):
            D['/masked/data/tau_V__Tg/%d' % iT] = ALL_tau_V__Tg[iT].data
            D['/masked/data/SFR__Tg/%d' % iT] = ALL_SFR__Tg[iT].data
            D['/masked/data/SFRSD__Tg/%d' % iT] = ALL_SFRSD__Tg[iT].data
            D['/masked/data/alogZ_mass__Tg/%d' % iT] = ALL_alogZ_mass__Tg[iT].data
            D['/masked/data/alogZ_flux__Tg/%d' % iT] = ALL_alogZ_flux__Tg[iT].data
            #D['/masked/data/at_flux__Tg/%d' % iT] = ALL_at_flux__Tg[iT].data
            D['/masked/mask/tau_V__Tg/%d' % iT] = ALL_tau_V__Tg[iT].mask
            D['/masked/mask/SFR__Tg/%d' % iT] = ALL_SFR__Tg[iT].mask
            D['/masked/mask/SFRSD__Tg/%d' % iT] = ALL_SFRSD__Tg[iT].mask
            D['/masked/mask/alogZ_mass__Tg/%d' % iT] = ALL_alogZ_mass__Tg[iT].mask
            D['/masked/mask/alogZ_flux__Tg/%d' % iT] = ALL_alogZ_flux__Tg[iT].mask
            #D['/masked/mask/at_flux__Tg/%d' % iT] = ALL_at_flux__Tg[iT].mask

        for k in D.keys():
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
                
        h5.close()
        
    print 'total time: %.2f' % (time.clock() - t_init_prog)