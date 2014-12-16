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
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#nebular_plot = True
nebular_plot = False
debug = False
#debug = True
#BPTLowS06 = True
BPTLowS06 = False
#hdf5 = False
hdf5 = True

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

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

RbinIni , RbinFin , RbinStep = 0.0 , 2.0 , 0.1
Rbin__r = np.arange(RbinIni, RbinFin + RbinStep, RbinStep)
RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
NRbins = len(RbinCenter__r)

RColor = [ 'r', 'y', 'b', 'k' ]
RRange = [  .5, 1., 1.5, 2.  ]

if debug:
    listOfPrefixes = listOfPrefixes[1:20]        # Q&D tests ...
    #listOfPrefixes = ['K0277\n']
    
N_gals = len(listOfPrefixes)

# SFR-time-scale array (index __T)
base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
tSF__T = base.ageBase
N_T = base.nAges

mask_xOk = True

# Def smallest light fraction (in the flag__t-ageMax age-range) deemed to be Ok for our stats ...
xOkMin = 0.05

# Minimum tauV to be taken seriously ...
tauVOkMin = 0.05
tauVNebOkMin = 0.05
tauVNebErrMax = 0.15


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


if __name__ == '__main__':
    ALL_N_zones__g = np.ma.zeros((N_gals))

    #########################################################################        
    ################################## CID ##################################
    ######################################################################### 
    ALL_morfType_GAL__g = np.ma.zeros((N_gals))
    ALL_at_flux_GAL__g = np.ma.zeros((N_gals))
    ALL_Mcor_GAL__g = np.ma.zeros((N_gals))
    ALL_McorSD_GAL__g = np.ma.zeros((N_gals))
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_morfType_GAL_zones__rg = np.ma.zeros((NRbins, N_gals))
    ALL_Mr_GAL_zones__rg = np.ma.zeros((NRbins, N_gals))
    ALL_ur_GAL_zones__rg = np.ma.zeros((NRbins, N_gals))
    ALL_tau_V_neb__rg = np.ma.zeros((NRbins, N_gals))
    ALL_alogZ_mass_GAL__Tg = np.ma.zeros((N_T, N_gals))
    ALL_alogZ_flux_GAL__Tg = np.ma.zeros((N_T, N_gals))
    ALL_isOkFrac_GAL__Tg = np.ma.zeros((N_T, N_gals))
    ALL_aSFRSD_Ha__rg = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD_Ha_kpc__rg = np.ma.zeros((NRbins, N_gals))
    ALL_aSFRSD__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_aSFRSD_kpc__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_tau_V__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_mass__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_alogZ_flux__Trg = np.ma.zeros((N_T, NRbins, N_gals))
    ALL_McorSD_GAL__rg = np.ma.zeros((NRbins, N_gals))
    
    ALL_integrated_SFR__Tg = np.ma.zeros((N_T, N_gals))
    ALL_integrated_SFRSD__Tg = np.ma.zeros((N_T, N_gals))
    ALL_integrated_SFRSD_kpc__Tg = np.ma.zeros((N_T, N_gals))
    ALL_integrated_SFR_Ha__g = np.ma.zeros((N_gals))
    ALL_integrated_SFRSD_Ha__g = np.ma.zeros((N_gals))
    ALL_integrated_SFRSD_Ha_kpc__g = np.ma.zeros((N_gals))
    ALL_integrated_L_int_Ha__g = np.ma.zeros((N_gals))
    
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
    _ALL_SFR_Ha__g = []
    _ALL_SFRSD_Ha__g = []
    _ALL_SFRSD_Ha_kpc__g = []
    _ALL_F_obs_Ha__g = []
    _ALL_L_int_Ha__g = []
    _ALL_L_int_Ha_err__g = []
    _ALL_L_int_Ha_mask__g = []
    _ALL_dist_zone__g = []
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # _ALL_at_flux__g = []
    # _ALL_aZ_flux__g = []
    # _ALL_alogZ_flux__g = []
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
    _ALL_SFRSD__Tg = []
    _ALL_SFRSD_kpc__Tg = []

    for iT in range(N_T):
        _ALL_tau_V__Tg.append([])
        _ALL_tau_V_mask__Tg.append([])
        _ALL_SFR__Tg.append([])
        _ALL_SFRSD__Tg.append([])
        _ALL_SFRSD_kpc__Tg.append([])
    #########################################################################
    #########################################################################
    #########################################################################

    ALL_N_zones_tau_V = 0
    ALL_N_gals = 0
    ALL_N_zones = 0
        
    for iGal in np.arange(N_gals):
        galName = listOfPrefixes[iGal][:-1]

        CALIFASuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        CALIFAFitsFile = superFitsDir + galName + CALIFASuffix
        emLinesSuffix = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
        
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
            ALL_alogZ_mass_GAL__Tg[:, iGal] = np.ma.masked
            ALL_alogZ_flux_GAL__Tg[:, iGal] = np.ma.masked
            ALL_isOkFrac_GAL__Tg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha_kpc__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD__Trg[:, :, iGal] = np.ma.masked
            ALL_aSFRSD_kpc__Trg[:, :, iGal] = np.ma.masked
            ALL_tau_V__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal] = np.ma.masked
            ALL_McorSD_GAL__rg[:, iGal] = np.ma.masked
            ALL_integrated_SFR__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD_kpc__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFR_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha_kpc__g[iGal] = np.ma.masked
            ALL_integrated_L_int_Ha__g[iGal] = np.ma.masked
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
        
        N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
        O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
        
        maskOk = HbOk & O3Ok & HaOk & N2Ok
        maskBPT = None
        
        if BPTLowS06:
            L = Lines()
            maskBPT = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        ##########################
        ##########################
        ##########################
        
        K.EL._forceMask = ~maskOk
        
        if len(K.EL.flux[0, :].compressed()) == 0:
            ALL_N_zones__g[iGal] = np.ma.masked
            ALL_morfType_GAL__g[iGal] = np.ma.masked
            ALL_at_flux_GAL__g[iGal] = np.ma.masked
            ALL_Mcor_GAL__g[iGal] = np.ma.masked
            ALL_McorSD_GAL__g[iGal] = np.ma.masked
            ALL_morfType_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_tau_V_neb__rg[:, iGal] = np.ma.masked
            ALL_Mr_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_ur_GAL_zones__rg[:, iGal] = np.ma.masked
            ALL_alogZ_mass_GAL__Tg[:, iGal] = np.ma.masked
            ALL_alogZ_flux_GAL__Tg[:, iGal] = np.ma.masked
            ALL_isOkFrac_GAL__Tg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD_Ha_kpc__rg[:, iGal] = np.ma.masked
            ALL_aSFRSD__Trg[:, :, iGal] = np.ma.masked
            ALL_aSFRSD_kpc__Trg[:, :, iGal] = np.ma.masked
            ALL_tau_V__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_mass__Trg[:, :, iGal] = np.ma.masked
            ALL_alogZ_flux__Trg[:, :, iGal] = np.ma.masked
            ALL_McorSD_GAL__rg[:, iGal] = np.ma.masked
            ALL_integrated_SFR__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFRSD_kpc__Tg[:, iGal] = np.ma.masked
            ALL_integrated_SFR_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha__g[iGal] = np.ma.masked
            ALL_integrated_SFRSD_Ha_kpc__g[iGal] = np.ma.masked
            ALL_integrated_L_int_Ha__g[iGal] = np.ma.masked
            
            print '<<< %s galaxy: no minSNR (%.1f) in 4 BPT lines' % (galName, minSNR)
            continue

        tipos, tipo, tipo_m, tipo_p = get_morfologia(galName)
        ALL_morfType_GAL__g[iGal] = tipo
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
        K.EL._forceMask = None
        
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
        
        # a fake all galaxy mass in stars per zone for all galaxies, creating a stamp for each zone
        aux = np.ones_like(K.Mcor__z) * ALL_Mcor_GAL__g[iGal]
        _ALL_Mcor_GAL_zones__g.append(aux)

        ALL_McorSD_GAL__rg[:, iGal] = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)

        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        ALL_at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()
        aux = np.ones_like(K.Mcor__z) * ALL_at_flux_GAL__g[iGal] 
        _ALL_at_flux_GAL_zones__g.append(aux)
        
        for iT, tSF in enumerate(tSF__T):
            #--------------------------------------------------------------------------
            # Define mask to pick only populations younger than the input tSF in the computation of SFR & SFRSD.
            flag__t = K.ageBase <= tSF
        
            # SRFSD "raw" image
            # Note that we are NOT dezonifying SFR__z, since it will be compared to the un-dezonifiable tauV!
            aux = K.Mini__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0) / tSF
            SFR__z = np.ma.masked_array(aux)
            SFRSD__z = SFR__z / K.zoneArea_pc2
            SFRSD_kpc__z = SFRSD__z * 1e6
            
            aux = K.integrated_Mini__tZ[flag__t, :].sum() / tSF
            integrated_SFR = np.ma.masked_array(aux)
            integrated_SFRSD = integrated_SFR / K.zoneArea_pc2.sum()
            integrated_SFRSD_kpc = integrated_SFRSD * 1e6
            
            tau_V__z = np.ma.masked_array(K.tau_V__z)
                            
            if mask_xOk:
                # Compute xOk "raw" image
                x__tZz = K.popx / K.popx.sum(axis = 1).sum(axis = 0)
                xOk__z = x__tZz[flag__t, :, :].sum(axis = 1).sum(axis = 0)
                
                maskNotOk__z = (xOk__z < xOkMin) | (tau_V__z < tauVOkMin) 
                 
                if BPTLowS06:
                    maskNotOk__z |= ~maskBPT

                tau_V__z[maskNotOk__z] = np.ma.masked
                SFR__z[maskNotOk__z] = np.ma.masked
                SFRSD__z[maskNotOk__z] = np.ma.masked
                SFRSD_kpc__z[maskNotOk__z] = np.ma.masked

            weiRadProf = True
                            
            aux = calc_alogZ_Stuff(K, flag__t, Rbin__r, weiRadProf, xOkMin = xOkMin)
            ALL_alogZ_mass_GAL__Tg[iT, iGal] = aux[0]
            ALL_alogZ_flux_GAL__Tg[iT, iGal] = aux[1]
            ALL_isOkFrac_GAL__Tg[iT, iGal] = aux[2]
            ALL_alogZ_mass__Trg[iT, :, iGal] = aux[3]
            ALL_alogZ_flux__Trg[iT, :, iGal] = aux[4]
            ALL_integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL_integrated_SFRSD__Tg[iT, iGal] = integrated_SFRSD
            ALL_integrated_SFRSD_kpc__Tg[iT, iGal] = integrated_SFRSD_kpc

            aSFRSD__yx = K.zoneToYX(SFR__z, extensive = True)
            aSFRSD__r = K.radialProfile(aSFRSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            aSFRSD_kpc__r = K.radialProfile(aSFRSD__yx * 1e6, Rbin__r, rad_scale = K.HLR_pix)

            ALL_aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
            ALL_aSFRSD_kpc__Trg[iT, :, iGal] = aSFRSD_kpc__r

            _ALL_tau_V__Tg[iT].append(tau_V__z.data)
            _ALL_tau_V_mask__Tg[iT].append(tau_V__z.mask)
            
            tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
            ALL_tau_V__Trg[iT, :, iGal] = K.radialProfile(tau_V__yx, Rbin__r, rad_scale = K.HLR_pix)
            
            _ALL_SFR__Tg[iT].append(SFR__z)
            _ALL_SFRSD__Tg[iT].append(SFRSD__z)
            _ALL_SFRSD_kpc__Tg[iT].append(SFRSD_kpc__z)
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
        SFRSD_Ha_kpc__z = SFRSD_Ha__z * 1e6
        
        integrated_SFR_Ha = 3.13 * integrated_L_int_Ha / (1.e8)
        integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()
        integrated_SFRSD_Ha_kpc = integrated_SFRSD_Ha * 1e6
        
        SFRSD_Ha__yx = K.zoneToYX(SFR_Ha__z, extensive = True)
        aSFRSD_Ha__r = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        aSFRSD_Ha_kpc__r = K.radialProfile(SFRSD_Ha__yx * 1e6, Rbin__r, rad_scale = K.HLR_pix)

        ALL_aSFRSD_Ha__rg[:, iGal] = aSFRSD_Ha__r
        ALL_aSFRSD_Ha_kpc__rg[:, iGal] = aSFRSD_Ha_kpc__r
        
        ALL_integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha
        ALL_integrated_SFRSD_Ha__g[iGal] = integrated_SFRSD_Ha
        ALL_integrated_SFRSD_Ha_kpc__g[iGal] = integrated_SFRSD_Ha_kpc
        
        _ALL_SFR_Ha__g.append(SFR_Ha__z)
        _ALL_SFRSD_Ha__g.append(SFRSD_Ha__z)
        _ALL_SFRSD_Ha_kpc__g.append(SFRSD_Ha__z * 1e6)
        ##########################
        ##########################
        ##########################
        
        if nebular_plot:
            f, axArr = plt.subplots(4, 4)
            f.set_size_inches(24, 20)
            
            for ax in f.axes:
                ax.set_axis_off()
            
            ax = axArr[0, 0]
            ax.set_axis_on()
            galaxyImgFile = imgDir + K.califaID + '.jpg'
            galimg = plt.imread(galaxyImgFile)
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.imshow(galimg)
            
            ax = axArr[0, 1]
            ax.set_axis_on()
            ax.plot(RbinCenter__r, tau_V_neb__r, 'o-k')
            #ax.tick_params(axis='x', pad=30)
            ax.set_title(r'$\tau_V^{neb}(R)$')
            
            ax = axArr[0, 2]
            ax.set_axis_on()
            L_int_Ha__yx = K.zoneToYX(L_int_Ha__z, extensive = True)
            L_int_Ha__r = K.radialProfile(L_int_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
            ax.plot(RbinCenter__r, L_int_Ha__r, 'o-k')
            ax.set_title(r'$L_{H\alpha}^{int}(R)$')
        
            ax = axArr[0, 3]
            ax.set_axis_on()
            #Lobn__yx            = K.zoneToYX(K.Lobn__z, extensive = True)
            #logZNeb__r          = radialProfileWeighted(logZ_neb__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
            logZ_neb__yx = K.zoneToYX(K.EL.logZ_neb_S06__z, extensive = False)
            logZNeb__r = K.radialProfile(logZ_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
            ax.plot(RbinCenter__r, logZNeb__r, 'o-k')
            #ax.set_title(r'$\langle \log\ Z_{neb}\rangle_L (HLR)$')
            ax.set_title(r'$\log\ Z_{neb}(R)$')
        
            ax = axArr[1, 1]
            ax.set_axis_on()
            ax.set_title(r'$\tau_V^{neb}$')
            sigma = tau_V_neb__yx.std()
            mean = tau_V_neb__yx.mean()
            vmax = mean + 2. * sigma
            vmin = mean - 2. * sigma
            im = ax.imshow(tau_V_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[2, 1]
            ax.set_axis_on()
            ax.set_title(r'$\epsilon (\tau_V^{neb})$')
            sigma = tau_V_neb_err__yx.std()
            mean = tau_V_neb_err__yx.mean()
            vmax = mean + 2. * sigma
            im = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[3, 1]
            ax.set_axis_on()
            ax.set_title(r'$F^{H\alpha}_{H\beta} / \epsilon(F^{H\alpha}_{H\beta})$')
            HaHb__yx = K.zoneToYX(K.EL.HaHb__z, extensive = True)
            err_HaHb__yx = K.zoneToYX(K.EL.HaHb_err__z, extensive = True)
            signalToNoise = np.abs(HaHb__yx) / np.abs(err_HaHb__yx) 
            sigma = signalToNoise.std()
            mean = signalToNoise.mean()
            vmax = mean + 2. * sigma
            im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
        
            ax = axArr[1, 2]
            ax.set_axis_on()
            ax.set_title(r'$L_{H\alpha}^{int}$')
            sigma = L_int_Ha__yx.std()
            mean = L_int_Ha__yx.mean()
            vmax = mean + sigma
            vmin = mean - sigma
            im = ax.imshow(L_int_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[2, 2]
            ax.set_axis_on()
            ax.set_title(r'$\epsilon (L_{H\alpha}^{int})$')
            L_int_Ha_err__yx = K.zoneToYX(L_int_Ha_err__z, extensive = True)
            sigma = L_int_Ha_err__yx.std()
            mean = L_int_Ha_err__yx.mean()
            vmax = mean + sigma
            im = ax.imshow(L_int_Ha_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[3, 2]
            ax.set_axis_on()
            ax.set_title(r'$L_{H\alpha}^{int} / \epsilon(L_{H\alpha}^{int})$')
            signalToNoise = np.abs(L_int_Ha__yx) / np.abs(L_int_Ha_err__yx) 
            sigma = signalToNoise.std()
            mean = signalToNoise.mean()
            vmax = mean + sigma
            im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
        
            ax = axArr[1, 3]
            ax.set_axis_on()
            ax.set_title(r'$\log\ Z_{neb}$')
            sigma = logZ_neb__yx.std()
            mean = logZ_neb__yx.mean()
            vmax = mean + sigma
            vmin = mean - sigma
            im = ax.imshow(logZ_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[2, 3]
            ax.set_axis_on()
            ax.set_title(r'$\epsilon (log\ Z_{neb})$')
            logZ_neb_err__yx = K.zoneToYX(K.EL.logZ_neb_S06_err__z, extensive = False)
            sigma = logZ_neb_err__yx.std()
            mean = logZ_neb_err__yx.mean()
            vmax = mean + sigma
            im = ax.imshow(logZ_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
            ax = axArr[3, 3]
            ax.set_axis_on()
            ax.set_title(r'$\log\ Z_{neb} / \epsilon(log\ Z_{neb})$')
            signalToNoise = np.abs(logZ_neb__yx) / np.abs(logZ_neb_err__yx) 
            sigma = signalToNoise.std()
            mean = signalToNoise.mean()
            vmax = mean + sigma
            im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
            f.colorbar(ax = ax, mappable = im, use_gridspec = False)
            
        #    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
            f.savefig(K.califaID + '_' + versionSuffix + '_nebular.png')
            plt.close(f)

        K.close()

    print 'Total of %d galaxies (%d zones): %s zones for tau_V' % (ALL_N_gals, ALL_N_zones, ALL_N_zones_tau_V)
    
    ALL_dist_zone__g = np.ma.masked_array(np.hstack(np.asarray(_ALL_dist_zone__g)))

    aux = np.hstack(_ALL_tau_V_neb__g)
    auxMask = np.hstack(_ALL_tau_V_neb_mask__g)
    ALL_tau_V_neb__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_tau_V_neb_err__g = np.ma.masked_array(np.hstack(_ALL_tau_V_neb_err__g), mask = auxMask)

    aux = np.hstack(_ALL_L_int_Ha__g)
    auxMask = np.hstack(_ALL_L_int_Ha_mask__g)
    ALL_L_int_Ha__g = np.ma.masked_array(aux, mask = auxMask)
    ALL_F_obs_Ha__g = np.ma.masked_array(np.hstack(_ALL_F_obs_Ha__g), mask = auxMask)
    ALL_SFR_Ha__g = np.ma.masked_array(np.hstack(_ALL_SFR_Ha__g), mask = auxMask)
    ALL_SFRSD_Ha__g = np.ma.masked_array(np.hstack(_ALL_SFRSD_Ha__g), mask = auxMask)
    ALL_SFRSD_Ha_kpc__g = np.ma.masked_array(np.hstack(_ALL_SFRSD_Ha_kpc__g), mask = auxMask)
    
    ALL_Mcor__g = np.ma.masked_array(np.hstack(_ALL_Mcor__g))
    ALL_McorSD__g = np.ma.masked_array(np.hstack(_ALL_McorSD__g))
    ALL_Mcor_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_Mcor_GAL_zones__g))
    ALL_McorSD_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_McorSD_GAL_zones__g))
    ALL_morfType_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_morfType_GAL_zones__g))
    ALL_at_flux_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_at_flux_GAL_zones__g))

    ALL_Mr_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_Mr_GAL_zones__g))
    ALL_ur_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_ur_GAL_zones__g))
    ALL_califaID_GAL_zones__g = np.ma.masked_array(np.hstack(_ALL_califaID_GAL_zones__g))

    ALL_tau_V__Tg = []
    ALL_SFR__Tg = []
    ALL_SFRSD__Tg = []
    ALL_SFRSD_kpc__Tg = []
    correl_SFR__T = np.ones_like(tSF__T)
    correl_SFRSD__T = np.ones_like(tSF__T)
    correl_aSFRSD__rT = np.ones((1 + len(RRange), tSF__T.shape[0]))

    for iT, tSF in enumerate(tSF__T):
        aux = np.hstack(_ALL_tau_V__Tg[iT])
        auxMask = np.hstack(_ALL_tau_V_mask__Tg[iT])
        ALL_tau_V__Tg.append(np.ma.masked_array(aux, mask = auxMask))
        ALL_SFR__Tg.append(np.ma.masked_array(np.hstack(_ALL_SFR__Tg[iT]), mask = auxMask))
        ALL_SFRSD__Tg.append(np.ma.masked_array(np.hstack(_ALL_SFRSD__Tg[iT]), mask = auxMask))
        ALL_SFRSD_kpc__Tg.append(np.ma.masked_array(np.hstack(_ALL_SFRSD_kpc__Tg[iT]), mask = auxMask))

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

    if hdf5:
        import h5py
        
        if BPTLowS06:       
            filename = 'SFR_BelowS06.h5'
        else:
            filename = 'SFR.h5'
            
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
            '/masked/data/aSFRSD_Ha_kpc__rg' : ALL_aSFRSD_Ha_kpc__rg.data,
            '/masked/data/aSFRSD__Trg' : ALL_aSFRSD__Trg.data,
            '/masked/data/aSFRSD_kpc__Trg' : ALL_aSFRSD_kpc__Trg.data,
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
            '/masked/data/SFRSD_Ha_kpc__g' : ALL_SFRSD_Ha_kpc__g.data,
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
            '/masked/data/integrated_SFRSD_kpc__Tg' : ALL_integrated_SFRSD_kpc__Tg.data,
            '/masked/data/integrated_SFR_Ha__g' : ALL_integrated_SFR_Ha__g.data,
            '/masked/data/integrated_SFRSD_Ha__g' : ALL_integrated_SFRSD_Ha__g.data,
            '/masked/data/integrated_SFRSD_Ha_kpc__g' : ALL_integrated_SFRSD_Ha_kpc__g.data,
            '/masked/data/dist_zone__g' : ALL_dist_zone__g.data,
            '/masked/data/N_zones__g' : ALL_N_zones__g.data, 
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
            '/masked/mask/aSFRSD_Ha_kpc__rg' : ALL_aSFRSD_Ha_kpc__rg.mask,
            '/masked/mask/aSFRSD__Trg' : ALL_aSFRSD__Trg.mask,
            '/masked/mask/aSFRSD_kpc__Trg' : ALL_aSFRSD_kpc__Trg.mask,
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
            '/masked/mask/SFRSD_Ha_kpc__g' : ALL_SFRSD_Ha_kpc__g.mask,
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
            '/masked/mask/integrated_SFRSD_kpc__Tg' : ALL_integrated_SFRSD_kpc__Tg.mask,
            '/masked/mask/integrated_SFR_Ha__g' : ALL_integrated_SFR_Ha__g.mask,
            '/masked/mask/integrated_SFRSD_Ha__g' : ALL_integrated_SFRSD_Ha__g.mask,
            '/masked/mask/integrated_SFRSD_Ha_kpc__g' : ALL_integrated_SFRSD_Ha_kpc__g.mask,
            '/masked/mask/dist_zone__g' : ALL_dist_zone__g.mask,
            '/masked/mask/N_zones__g' : ALL_N_zones__g.mask,
            '/data/correl_SFR__T' : correl_SFR__T,
            '/data/correl_SFRSD__T' : correl_SFRSD__T,
            '/data/correl_aSFRSD__rT' : correl_aSFRSD__rT,
            '/data/zones_tau_V' : ALL_N_zones_tau_V,
            '/data/gals' : ALL_N_gals,
            '/data/zones' : ALL_N_zones,
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
            D['/masked/data/SFRSD_kpc__Tg/%d' % iT] = ALL_SFRSD_kpc__Tg[iT].data
            D['/masked/mask/tau_V__Tg/%d' % iT] = ALL_tau_V__Tg[iT].mask
            D['/masked/mask/SFR__Tg/%d' % iT] = ALL_SFR__Tg[iT].mask
            D['/masked/mask/SFRSD__Tg/%d' % iT] = ALL_SFRSD__Tg[iT].mask
            D['/masked/mask/SFRSD_kpc__Tg/%d' % iT] = ALL_SFRSD_kpc__Tg[iT].mask

        for k in D.keys():
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
                
        h5.close()
