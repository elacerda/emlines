#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
from rgbread import read_rgb_fits
import numpy as np
from pycasso import fitsQ3DataCube
import pyfits

ZSol = 0.019
LSol = 3.826e33 # erg/s
qCCM = {
    'Hb' : 1.16427,
    'O3' : 1.12022,
    'Ha' : 0.81775,
    'N2' : 0.81466,
}

gal = 'K0277'
baseCode = 'Bgsd6e'

superFitsDir = '/Users/lacerda/CALIFA/gal_fits/px1_q043.d14a/'
CALIFASuffix = '_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.' + baseCode + '.fits'
CALIFAFitsFile = superFitsDir + gal + CALIFASuffix

emLinesFitsDir = '/Users/lacerda/CALIFA/rgb-gas/'
emLinesSuffix = '_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.' + baseCode + '.EML.MC100.fits'
emLinesFitsFile = emLinesFitsDir + gal + emLinesSuffix

#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
def Mpc_to_cm(dist):
    # google: 1 Mpc = 3.08567758e24 cm
    return dist * 3.08567758e24 
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
def flux_to_LSol(flux, distance):
    return 4. * np.pi * Mpc_to_cm(distance) ** 2.0 * flux / LSol
#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#####################################################################################
#####################################################################################
#####################################################################################
if __name__ == '__main__':
    K = fitsQ3DataCube(CALIFAFitsFile)
    
    # read FITSFILE containing galaxy emission lines measured by R.G.B.
    flux__z, err_flux__z, fwhm__z, err_fwhm__z, ew__z = read_rgb_fits(emLinesFitsFile)
    
    lines = {
        'Hb' : '4861',
        'O3' : '5007',
        'Ha' : '6563',
        'N2' : '6583',
    }
    
    solidAngle = 4. * np.pi * K.distance_Mpc
    
    Lobs = {}
    err_Lobs = {}
    
    for line in lines.keys():
        Lobs[line] = flux_to_LSol(flux__z[lines[line]], K.distance_Mpc)
        err_Lobs[line] = flux_to_LSol(err_flux__z[lines[line]], K.distance_Mpc)
    
    LintHaHb = 2.86
    LobsHaHb = Lobs['Ha'] / Lobs['Hb']
    
    tauVNeb__z = np.ma.log(LobsHaHb / LintHaHb) / (qCCM['Hb'] - qCCM['Ha'])
    tauVNeb__yx = K.zoneToYX(tauVNeb__z, extensive = False)
    
    e = np.ma.exp(qCCM['Ha'] * tauVNeb__z)
    q = qCCM['Ha'] / (qCCM['Hb'] - qCCM['Ha'])
    
    Lint_Ha__z = Lobs['Ha'] * e
    Lint_Ha__yx = K.zoneToYX(Lint_Ha__z, extensive = True)
     
    err_Lint_Ha__z = e * np.sqrt(err_Lobs['Ha'] ** 2.0 * q ** 2.0 * LobsHaHb ** 2.0 * err_Lobs['Hb'] ** 2.0)
    err_Lint_Ha__yx = K.zoneToYX(err_Lint_Ha__z, extensive = True)
    
    err_tauVNeb__z = np.sqrt((err_Lobs['Ha'] / Lobs['Ha']) ** 2.0 + (err_Lobs['Hb'] / Lobs['Hb']) ** 2.0) / (qCCM['Hb'] - qCCM['Ha'])
    err_tauVNeb__yx = K.zoneToYX(err_tauVNeb__z, extensive = False)
