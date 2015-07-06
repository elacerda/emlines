import sys
import numpy as np
import argparse as ap
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
import CALIFAUtils as C
from CALIFAUtils.lines import Lines

def parser_args():
    paths = C.paths
    default_args = {
        'debug' : False,
        'hdf5' : None,
        'mintauvneb' : 0.05,
        'maxtauvneberr' : 999.,
        'gals_filename' : paths.califa_work_dir + 'listv20_q050.d15a.txt',
        'rgbcuts' : False,
        'gasprop' : False,
        'v_run' :-1,
    }
    
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--gasprop', '-G',
                        action = 'store_true',
                        default = default_args['gasprop'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
    parser.add_argument('--v_run',
                        metavar = 'INT',
                        type = int,
                        default = default_args['v_run'])
    parser.add_argument('--gals_filename', '-L',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['gals_filename'])
    parser.add_argument('--rgbcuts',
                        action = 'store_true',
                        default = default_args['rgbcuts'])
    parser.add_argument('--mintauvneb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauvneb'])
    parser.add_argument('--maxtauvneberr',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['maxtauvneberr'])

    return parser.parse_args()

def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return 0, False
    if EL == True and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return 1, False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return 2, False
    if GP is True and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return 2, False
    # Problem in FITS file
    return 0, True

def O_O3N2_M13(Ha_obs, Hb_obs, O3_obs, N2_obs, mask_zones, tau_V = None, correct = False):
    Hb = np.ma.masked_array(Hb_obs, mask = mask_zones)
    O3 = np.ma.masked_array(O3_obs, mask = mask_zones)
    Ha = np.ma.masked_array(Ha_obs, mask = mask_zones)
    N2 = np.ma.masked_array(N2_obs, mask = mask_zones)
    if correct is True:
        tau_V_m = np.ma.masked_array(tau_V, mask = mask_zones)
        from pystarlight.util import redenninglaws
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        Hb *= np.ma.exp(q[0] * tau_V_m) 
        O3 *= np.ma.exp(q[1] * tau_V_m) 
        Ha *= np.ma.exp(q[2] * tau_V_m) 
        N2 *= np.ma.exp(q[3] * tau_V_m)
    return np.ma.log10(O3 * Ha / (N2 * Hb)) 

def O3N2(Ha_obs, Hb_obs, O3_obs, N2_obs, mask_zones, tau_V = None, correct = False):
    Hb = np.ma.masked_array(Hb_obs, mask = mask_zones)
    O3 = np.ma.masked_array(O3_obs, mask = mask_zones)
    Ha = np.ma.masked_array(Ha_obs, mask = mask_zones)
    N2 = np.ma.masked_array(N2_obs, mask = mask_zones)
    if correct is True:
        tau_V_m = np.ma.masked_array(tau_V, mask = mask_zones)
        from pystarlight.util import redenninglaws
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        Hb *= np.ma.exp(q[0] * tau_V_m) 
        O3 *= np.ma.exp(q[1] * tau_V_m) 
        Ha *= np.ma.exp(q[2] * tau_V_m) 
        N2 *= np.ma.exp(q[3] * tau_V_m)
    O3Hb = np.ma.log10(O3/Hb)
    N2Ha = np.ma.log10(N2/Ha)
    O3N2 = np.ma.log10(O3 * Ha / (N2 * Hb))
    return O3Hb, N2Ha, O3N2

if __name__ == '__main__':
    # Parse arguments 
    args = parser_args()
    C.debug_var(args.debug, args = args)    
    
    Zsun = 0.019
    
    # Reading galaxies file,
    gals, _ = C.sort_gals(args.gals_filename)
    N_gals = len(gals)
    maxGals = None
    if args.debug:
        maxGals = 10
        if N_gals > maxGals:
            N_gals = maxGals
            
    _tau_V_neb__g = []
    _tau_V_neb_RGB__g = []
    _O_O3N2_M13__g = []
    _O_O3N2_c__g = [] # M13 computed by myself
    _O_O3N2_nc__g = [] # M13 computed by myself without extinction correction
    _O_O3N2_RGB_c__g = [] # M13 computed by myself with extinction correction by RGB
    _O3N2_c__g = [] # M13 computed by myself
    _O3N2_nc__g = [] # M13 computed by myself without extinction correction
    _O3N2_RGB_c__g = [] # M13 computed by myself with extinction correction by RGB
    _maskOkTauVNeb__g = []
    _maskOkTauVNebRGB__g = []
    _maskOkLines__g = []
    _O_O3N2_M13_mask__g = []
    _N_zones__g = []
    N_zones__G = np.ma.empty((N_gals))
    califaIDs__G = np.ma.empty((N_gals), dtype = '|S5')     

    for iGal, K in C.loop_cubes(gals.tolist(), imax = maxGals, EL = True, GP = args.gasprop, v_run = args.v_run):        
    #for iGal in xrange(len(gals)):
        califaID = gals[iGal] 
        
        sit, verify = verify_files(K, califaID, EL = True, GP = args.gasprop)
        
        if verify is not True:
            N_zones__G[iGal] = np.ma.masked
            califaIDs__G[iGal] = np.ma.masked
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue
        
        califaIDs__G[iGal] = califaID 
        N_zones__G[iGal] = K.N_zone
        _N_zones__g.append(np.arange(K.N_zone))
        
        ####################################################
        ######## EmLines ###################################
        ####################################################
        Hb_central_wl = '4861'
        O3_central_wl = '5007'
        Ha_central_wl = '6563'
        N2_central_wl = '6583'
        
        lines_central_wl = [ 
            Hb_central_wl,
            O3_central_wl,
            Ha_central_wl,
            N2_central_wl,
        ]
         
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)
        
        ###### MASK EmLines ######
        maskOkLine = {}
        if args.rgbcuts is True:
            for l in lines_central_wl:
                if args.gasprop is True:
                    pos, sigma, snr = K.GP._dlcons[l]['pos'], K.GP._dlcons[l]['sigma'], K.GP._dlcons[l]['SN']
                else:
                    pos, sigma, snr = 3.0, 3.0, 3.0
                C.debug_var(args.debug, pref = l, pos = pos, sigma = sigma, snr = snr)
                maskOkLine[l] = K.EL._setMaskLineFluxNeg(l) 
                maskOkLine[l] |= K.EL._setMaskLineDisplacement(l, pos)
                maskOkLine[l] |= K.EL._setMaskLineSigma(l, sigma)
                maskOkLine[l] |= K.EL._setMaskLineSNR(l, snr)
                maskOkLine[l] = ~(maskOkLine[l])
                #C.debug_var(args.debug, maskOkLine = maskOkLine[l])
        else:
            for l in lines_central_wl:
                C.debug_var(args.debug, l = l)
                minSNR = 3
                maskOkLine[l] = K.EL._setMaskLineFluxNeg(l)
                maskOkLine[l] |= K.EL._setMaskLineSNR(l, minSNR)
                maskOkLine[l] = ~(maskOkLine[l])
                #C.debug_var(args.debug, maskOkLine = maskOkLine[l])
            
        maskOkLines__z = np.bitwise_and(maskOkLine[Hb_central_wl], maskOkLine[O3_central_wl])
        maskOkLines__z = np.bitwise_and(maskOkLines__z, maskOkLine[Ha_central_wl])
        maskOkLines__z = np.bitwise_and(maskOkLines__z, maskOkLine[N2_central_wl])  

        ########## tau_V #########
        maskOkTauVNeb__z = np.ones((K.N_zone), dtype = np.bool)
        if args.mintauvneb >= 0:
            maskOkTauVNeb__z = (K.EL.tau_V_neb__z >= args.mintauvneb) & (K.EL.tau_V_neb_err__z <= args.maxtauvneberr)
        maskOkTauVNeb__z &= maskOkLine[Ha_central_wl]
        maskOkTauVNeb__z &= maskOkLine[Hb_central_wl]
        
        tau_V_neb_RGB__z = K.GP.AVtoTau(K.GP.REDDENING.AV)
        tau_V_neb_RGB_err__z = K.GP.AVtoTau(K.GP.REDDENING.e_AV)
        maskOkTauVNebRGB__z = (tau_V_neb_RGB__z >= 0.)

        O3N2_nc__z = O_O3N2_M13(K.EL.Ha_obs__z, K.EL.Hb_obs__z,
                                  K.EL.O3_obs__z, K.EL.N2_obs__z,
                                  mask_zones = ~(maskOkLines__z))
        O3N2_c__z = O_O3N2_M13(K.EL.Ha_obs__z, K.EL.Hb_obs__z,
                                 K.EL.O3_obs__z, K.EL.N2_obs__z,
                                 mask_zones = ~(maskOkTauVNeb__z & maskOkLines__z),
                                 tau_V = K.EL.tau_V_neb__z, correct = True)
        O3N2_c_RGB__z = O_O3N2_M13(K.EL.Ha_obs__z, K.EL.Hb_obs__z,
                                     K.EL.O3_obs__z, K.EL.N2_obs__z,
                                     mask_zones = ~(maskOkTauVNebRGB__z & maskOkLines__z),
                                     tau_V = tau_V_neb_RGB__z, correct = True)

        O_O3N2_nc__z = 8.533 - 0.214 * O3N2_nc__z
        O_O3N2_c__z = 8.533 - 0.214 * O3N2_c__z
        O_O3N2_c_RGB__z = 8.533 - 0.214 * O3N2_c_RGB__z
        
        _tau_V_neb__g.append(K.EL.tau_V_neb__z)
        _tau_V_neb_RGB__g.append(tau_V_neb_RGB__z)
        _O3N2_c__g.append(O3N2_c__z)
        _O3N2_nc__g.append(O3N2_nc__z)
        _O3N2_RGB_c__g.append(O3N2_c_RGB__z)
        _O_O3N2_c__g.append(O_O3N2_c__z)
        _O_O3N2_nc__g.append(O_O3N2_nc__z)
        _O_O3N2_RGB_c__g.append(O_O3N2_c_RGB__z)
        _O_O3N2_M13__g.append(K.GP.EMPAB.O_O3N2_M13)
        _maskOkTauVNeb__g.append(maskOkTauVNeb__z)
        _maskOkTauVNebRGB__g.append(maskOkTauVNebRGB__z)
        _maskOkLines__g.append(maskOkLines__z)
        _O_O3N2_M13_mask__g.append(np.isnan(K.GP.EMPAB.O_O3N2_M13))

    N_zones__g = np.hstack(_N_zones__g)    
    tau_V_neb__g = np.hstack(_tau_V_neb__g)
    tau_V_neb_RGB__g = np.hstack(_tau_V_neb_RGB__g)
    O3N2_c__g = np.hstack(_O3N2_c__g)
    O3N2_nc__g = np.hstack(_O3N2_nc__g)
    O3N2_RGB_c__g = np.hstack(_O3N2_RGB_c__g)
    O_O3N2_c__g = np.hstack(_O_O3N2_c__g)
    O_O3N2_nc__g = np.hstack(_O_O3N2_nc__g)
    O_O3N2_RGB_c__g = np.hstack(_O_O3N2_RGB_c__g)
    O_O3N2_M13__g = np.hstack(_O_O3N2_M13__g)
    maskOkTauVNeb__g = np.hstack(_maskOkTauVNeb__g)
    maskOkTauVNebRGB__g = np.hstack(_maskOkTauVNebRGB__g)
    maskOkLines__g = np.hstack(_maskOkLines__g)
    O_O3N2_M13_mask__g = np.hstack(_O_O3N2_M13_mask__g)
