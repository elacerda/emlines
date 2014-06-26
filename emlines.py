#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
from rgbread import read_rgb_fits
import numpy as np
from pycasso import fitsQ3DataCube
import pyfits
import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
import os

plot = True
#plot = False
debug = False
#debug = True

mpl.rcParams['font.size'] = 16
mpl.rcParams['font.family'] = 'sans-serif'
    
galaxiesListFile    = '/Users/lacerda/CALIFA/listOf300GalPrefixes.txt'
baseCode            = 'Bgsd6e'
versionSuffix      = 'px1_q043.d14a'
#versionSuffix       = 'v20_q043.d14a'
superFitsDir        = '/Users/lacerda/CALIFA/gal_fits/' + versionSuffix + '/'
emLinesFitsDir      = '/Users/lacerda/CALIFA/rgb-gas/' + versionSuffix + '/'

ZSol = 0.019
LSol = 3.826e33 # erg/s
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
    listOfPrefixes = ['K0001\n']
    
N_gals = len(listOfPrefixes)

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
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / LSol

def calc_Lobs(f_obs__lz, distance_Mpc):
    '''
    Calculate luminosity using 
    L_\lambda = 4 \pi d^2 F_\lambda
    Distance in Mpc due to flux_to_LSol() uses this unit.
    
    Lacerda@Saco - 25/Jun/2014
    '''
    solidAngle = 4. * np.pi * distance_Mpc
    
    Lobs__lz        = np.ma.zeros(f_obs__lz.shape)
    err_Lobs__lz    = np.ma.zeros(f_obs__lz.shape)
    
    for line in range(f_obs__lz.shape[0]):
        Lobs__lz[line]      = flux_to_LSol(f_obs__lz[line, :], K.distance_Mpc)
        err_Lobs__lz[line]  = flux_to_LSol(err_f_obs__lz[line, :], K.distance_Mpc)
        
    return Lobs__lz, err_Lobs__lz

def calc_tauVNeb(K, f_obs__lz, err_f_obs__lz, lines):
    '''
    Calculate Balmer optical depth (tau_V).
    
    Lacerda@Saco - 25/Jun/2014
    '''
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    f_int_HaHb = 2.86
    f_obs_HaHb__z = f_obs__lz[i_Ha, :] / f_obs__lz[i_Hb, :]
    q = qCCM['4861'] - qCCM['6563']
    #LintHaHb = 2.86
    #LobsHaHb = Lobs['6563'] / Lobs['4861'] # OBS: LHa / LHb = fluxHa / fluxHb 
    
    #tauVNeb__z = np.ma.log(LobsHaHb / LintHaHb) / q
    lnHaHb = np.ma.log(f_obs_HaHb__z / f_int_HaHb)
    tauVNeb__z = lnHaHb / q 
    
    #err_tauVNeb__z = np.sqrt((err_Lobs['6563'] / Lobs['6563']) ** 2.0 + (err_Lobs['4861'] / Lobs['4861']) ** 2.0) / (qCCM['4861'] - qCCM['6563'])
    a = err_f_obs__lz[i_Ha, :] / f_obs__lz[i_Ha, :]
    b = err_f_obs__lz[i_Hb, :] / f_obs__lz[i_Hb, :]
    err_tauVNeb__z = np.sqrt(a ** 2.0 + b ** 2.0) / q
    
    a                 = err_f_obs__lz[i_Ha, :]
    b                 = (f_obs__lz[i_Ha, :] / f_obs__lz[i_Hb, :]) * err_f_obs__lz[i_Hb, :]
    err_f_obs_HaHb__z  = np.ma.sqrt(a ** 2. + b ** 2.) / f_obs__lz[i_Hb, :]
        
    return tauVNeb__z, err_tauVNeb__z, f_obs_HaHb__z, err_f_obs_HaHb__z

def calc_Lint_Ha(Lobs__lz, err_Lobs__lz, tauVNeb__z,lines):
    i_Ha = lines.index('6563')
    i_Hb = lines.index('4861')
    
    q = qCCM['6563'] / (qCCM['4861'] - qCCM['6563'])
    
    eHa = np.ma.exp(qCCM['6563'] * tauVNeb__z)
    LobsHaHb = Lobs__lz[i_Ha, :] / Lobs__lz[i_Hb, :]

    Lint_Ha__z = Lobs__lz[i_Ha, :] * eHa
    
    a = err_Lobs__lz[i_Ha, :]
    b = q * LobsHaHb * err_Lobs__lz[i_Hb, :]
    err_Lint_Ha__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
    
    return Lint_Ha__z, err_Lint_Ha__z

def calc_logZNeb(K, f_obs__lz, err_f_obs__lz, lines):
    '''
    Calculates Z_Neb using Asari et al (2007)
    Z_Neb = log((O/H) / (O/H)_Sol) = - 0.14 - 0.25 log([OIII]5007 / [NII]6583)
    
    Lacerda@Saco - 25/Jun/2014
    '''
    i_O3 = lines.index('5007')
    i_N2 = lines.index('6583')

    O3N2__z     = f_obs__lz[i_O3, :] / f_obs__lz[i_N2, :]
    logZNeb__z  = - 0.14 - (0.25 * np.log10(O3N2__z))
    
    a               = err_f_obs__lz[i_O3, :] / f_obs__lz[i_O3, :]
    b               = err_f_obs__lz[i_N2, :] / f_obs__lz[i_N2, :]
    err_logZNeb__z  = 0.25 * np.sqrt(a ** 2. + b ** 2.)

    a           = err_f_obs__lz[i_O3, :]
    b           = (f_obs__lz[i_O3, :] / f_obs__lz[i_N2, :]) * err_f_obs__lz[i_N2, :]
    err_O3N2__z = np.ma.sqrt(a ** 2. + b ** 2.) / f_obs__lz[i_N2, :]
    
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

if __name__ == '__main__':
    for iGal in np.arange(N_gals):
        galName         = listOfPrefixes[iGal][:-1]

        CALIFASuffix    = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.fits'
        CALIFAFitsFile  = superFitsDir + galName + CALIFASuffix
        emLinesSuffix   = '_synthesis_eBR_' + versionSuffix + '512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
        emLinesFitsFile = emLinesFitsDir + galName + emLinesSuffix
        
        if not os.path.isfile(CALIFAFitsFile):
            continue
        
        K = fitsQ3DataCube(CALIFAFitsFile)
        
        # read FITSFILE containing galaxy emission lines measured by R.G.B.
        # read_rgb_fits returns False if emLinesFitsFile does not exists.
        read = read_rgb_fits(emLinesFitsFile, read_lines)
        
        if read:
            print '>>> Doing' , iGal , galName , '|  Nzones=' , K.N_zone
            
            f_obs__lz       = read[0]
            err_f_obs__lz   = read[1]
            fwhm__lz        = read[2]
            err_fwhm__lz    = read[3]
            ew__lz          = read[4]
            lines           = read[5]

            #################################################
            ##################### MASK ######################
            #################################################            
            # minimum value of f_lz / err_f_lz
            minSNR = 3.
            
            i_Hb = lines.index('4861')
            i_O3 = lines.index('5007')
            i_Ha = lines.index('6563')
            i_N2 = lines.index('6583')
            HbOk = (f_obs__lz[i_Hb, :] / err_f_obs__lz[i_Hb, :]) >= minSNR
            O3Ok = (f_obs__lz[i_O3, :] / err_f_obs__lz[i_O3, :]) >= minSNR
            HaOk = (f_obs__lz[i_Ha, :] / err_f_obs__lz[i_Ha, :]) >= minSNR
            N2Ok = (f_obs__lz[i_N2, :] / err_f_obs__lz[i_N2, :]) >= minSNR
            maskOk = HbOk & O3Ok & HaOk & N2Ok
            
            f_obs__lz[:, ~maskOk] = np.ma.masked
            err_f_obs__lz[:, ~maskOk] = np.ma.masked
            ##################################################
            
            aux                 = calc_Lobs(f_obs__lz, K.distance_Mpc)
            Lobs__lz            = aux[0]
            err_Lobs__lz        = aux[1]
            
            aux                 = calc_tauVNeb(K, f_obs__lz, err_f_obs__lz, lines)
            tauVNeb__z          = aux[0] 
            err_tauVNeb__z      = aux[1] 
            f_obs_HaHb__z       = aux[2] 
            err_f_obs_HaHb__z   = aux[3]
            tauVNeb__yx         = K.zoneToYX(tauVNeb__z, extensive = False)
            err_tauVNeb__yx     = K.zoneToYX(err_tauVNeb__z, extensive = False)
            f_obs_HaHb__yx      = K.zoneToYX(f_obs_HaHb__z, extensive = True)
            err_f_obs_HaHb__yx  = K.zoneToYX(err_f_obs_HaHb__z, extensive = True)
        
            aux                 = tauV_to_AV(tauVNeb__z, err_tauVNeb__z)
            AVNeb__z            = aux[0]
            err_AVNeb__z        = aux[1]         
            AVNeb__yx           = K.zoneToYX(AVNeb__z, extensive = False)
            err_AVNeb__yx       = K.zoneToYX(err_AVNeb__z, extensive = False)
        
            aux                 = calc_Lint_Ha(Lobs__lz, err_Lobs__lz, tauVNeb__z, lines)
            Lint_Ha__z          = aux[0]
            err_Lint_Ha__z      = aux[1] 
            Lint_Ha__yx         = K.zoneToYX(Lint_Ha__z, extensive = True)
            err_Lint_Ha__yx     = K.zoneToYX(err_Lint_Ha__z, extensive = True)
            
            aux                 = calc_logZNeb(K, f_obs__lz, err_f_obs__lz, lines)
            logZNeb__z          = aux[0]
            err_logZNeb__z      = aux[1]
            O3N2__z             = aux[2]
            err_O3N2__z         = aux[3]
            logZNeb__yx         = K.zoneToYX(logZNeb__z, extensive = False)
            err_logZNeb__yx     = K.zoneToYX(err_logZNeb__z, extensive = False)
            O3N2__yx            = K.zoneToYX(O3N2__z, extensive = True)
            err_O3N2__yx        = K.zoneToYX(err_O3N2__z, extensive = True)
            
            if plot:
                tauVNeb_implot(tauVNeb__yx, err_tauVNeb__yx, f_obs_HaHb__yx, err_f_obs_HaHb__yx, K.califaID + '_' + versionSuffix + '-tauVNeb.png')
                Lint_Ha_implot(Lint_Ha__yx, err_Lint_Ha__yx, K.califaID + '_' + versionSuffix + '-Lint.png')
            