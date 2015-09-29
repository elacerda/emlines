#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import time
import sys
import argparse as ap
import numpy as np
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
import CALIFAUtils as C
from CALIFAUtils.lines import Lines
from CALIFAUtils.objects import ALLGals
from CALIFAUtils.scripts import calc_SFR
from CALIFAUtils.scripts import calc_xY
from CALIFAUtils.scripts import calc_alogZ_Stuff
from CALIFAUtils.scripts import radialProfileWeighted

def parser_args():
    paths = C.paths
    default_args = {
        'debug' : False,
        'underS06' : False,
        'whanSF' : False,
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
        'maxiTSF' :-1,
        'gals_filename' : paths.califa_work_dir + 'listv20_q050.d15a.txt',
        'rgbcuts' : False,
        'filter_residual' : False,
        'gasprop' : False,
        'v_run' : -1,
    }
    
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--gasprop', '-G',
                        action = 'store_true',
                        default = default_args['gasprop'])
    parser.add_argument('--spiral', '-S',
                        action = 'store_true',
                        default = default_args['spiral'])
    parser.add_argument('--filter_residual', '-R',
                        action = 'store_true',
                        default = default_args['filter_residual'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default_args['weiradprof'])
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
    parser.add_argument('--underS06',
                        action = 'store_true',
                        default = default_args['underS06'])
    parser.add_argument('--whanSF',
                        action = 'store_true',
                        default = default_args['whanSF'])
    parser.add_argument('--rgbcuts',
                        action = 'store_true',
                        default = default_args['rgbcuts'])
    parser.add_argument('--minpopx',
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minpopx'])
    parser.add_argument('--mintauv',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauv'])
    parser.add_argument('--mintauvneb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauvneb'])
    parser.add_argument('--maxtauvneberr',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['maxtauvneberr'])
    parser.add_argument('--rbinini',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinini'])
    parser.add_argument('--rbinfin',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinfin'])
    parser.add_argument('--rbinstep',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['rbinstep'])
    parser.add_argument('--maxiTSF',
                        metavar = 'iT',
                        type = int,
                        default = default_args['maxiTSF'])

    return parser.parse_args()

def print_args(args):
    # dummy function to print a dictionary.
    for k, v in args.__dict__.iteritems():
        print k, v 

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

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parser_args()
    print_args(args)    
    
    Zsun = 0.019

    # Creating radial bins.
    Rbin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins = len(RbinCenter__r)
    RColor = [ 'r', 'y', 'b', 'k' ]
    RRange = [  .5, 1., 1.5, 2.  ]
    
    # Reading galaxies file,
    gals, _ = C.sort_gals(gals = args.gals_filename, order = 1)
    N_gals = len(gals)
    maxGals = None
    if args.debug:
        maxGals = 10
        if N_gals > maxGals:
            N_gals = maxGals

    # SFR-time-scale array (index __T)
    base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # if args.maxiTSF > 0:
    #     tSF__T = np.asarray(base.ageBase[0:args.maxiTSF + 1])
    # else:
    #     tSF__T = np.asarray(base.ageBase)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    #tSF__T = np.array([0.032 , 0.3 , 1.5, 14.2]) * 1.e9
    tSF__T = np.array([1, 2.5, 3.2, 10]) * 1.e7
    #tSF__T = np.logspace(np.log10(base.ageBase[1]), np.log10(base.ageBase[23]), 100)
    N_T = len(tSF__T)

    # Z-time-scale array (index __U).
    tZ__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
    N_U = len(tZ__U)

    ALL = ALLGals(N_gals, NRbins, N_T, N_U)

    # automatically read PyCASSO and EmLines data cubes.
    for iGal, K in C.loop_cubes(gals.tolist(), imax = maxGals, EL = True, GP = args.gasprop, v_run = args.v_run):        
    #for iGal in xrange(len(gals)):
        t_init_gal = time.clock()
        califaID = gals[iGal] 
        #K =  C.read_one_cube(califaID, EL = True, GP = args.gasprop, v_run = args.v_run)
        
        # Saving for later :D        
        ALL.califaIDs[iGal] = califaID
        
        sit, verify = verify_files(K, califaID, EL = True, GP = args.gasprop)
        
        if verify is not True:
            ALL.mask_gal(iGal)
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        tipos, tipo, tipo_m, tipo_p = C.get_morfologia(califaID)
        
        # Only spiral
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # [
        # ('E0', 0),
        # ('E1', 1),
        # ('E2', 2),
        # ('E3', 3),
        # ('E4', 4),
        # ('E5', 5),
        # ('E6', 6),
        # ('E7', 7),
        # ('S0', 8),
        # ('S0a', 8.5),
        # ('Sa', 9),
        # ('Sab', 9.5),
        # ('Sb', 10),
        # ('Sbc', 10.5),
        # ('Sc', 11),
        # ('Scd', 11.5),
        # ('Sd', 12),
        # ('Sdm', 12.5),
        # ('Sm', 13),
        # ('Ir', 14),
        # ]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        if args.spiral and (tipo <= 8.5 or tipo >= 12.5):  # between Sa ... Sd
            ALL.mask_gal(iGal)
            if args.gasprop is True:
                K.GP.close()
            K.EL.close()
            K.close()
            print '<<< %s galaxy: is not a spiral (type: %f)' % (califaID, tipo) 
            continue

        N_zone = K.N_zone

        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        # Saving for later :D
        ALL.N_zones__g[iGal] = N_zone
        ALL.ba_PyCASSO_GAL__g[iGal] = ba
        
        print '>>> Doing' , iGal , califaID , 'hubtyp=', tipo, '|  Nzones=' , N_zone
        
        # zone distance in HLR
        #zoneDistHLR = np.sqrt((K.zonePos['x'] - K.x0) ** 2. + (K.zonePos['y'] - K.y0) ** 2.) / K.HLR_pix
        zone_distance_HLR = K.zoneDistance_HLR
        zone_area_pc2 = K.zoneArea_pc2
        ba_GAL = np.float(K.masterListData['ba'])
        Mr_GAL = np.float(K.masterListData['Mr'])
        ur_GAL = np.float(K.masterListData['u-r'])
        
        # Saving for later :D
        ALL._zone_dist_HLR__g.append(zone_distance_HLR)
        ALL._zone_area_pc2__g.append(zone_area_pc2)
        ALL.ba_GAL__g[iGal] = ba_GAL
        ALL.Mr_GAL__g[iGal] = Mr_GAL
        ALL.ur_GAL__g[iGal] = ur_GAL
        ALL.morfType_GAL__g[iGal] = tipo
        
        ####################################################
        ####### STARLIGHT ##################################
        ####################################################
        AVtotauV = 1. / (np.log10(np.exp(1)) / 0.4)
        ALL.integrated_tau_V__g[iGal] = K.integrated_keywords['A_V'] * AVtotauV 
        
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        Mcor__z = K.Mcor__z
        McorSD__z = K.Mcor__z / K.zoneArea_pc2
        Mcor_GAL = K.Mcor_tot.sum()
        McorSD_GAL = K.McorSD__yx.mean()
        McorSD__r = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale=K.HLR_pix)
        at_flux__z = K.at_flux__z
        at_mass__z = K.at_mass__z
        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        at_flux_GAL = numerator__z.sum() / denominator__z.sum()
        at_flux__r = K.radialProfile(K.at_flux__yx, Rbin__r, rad_scale=K.HLR_pix)
        at_mass__r = K.radialProfile(K.at_mass__yx, Rbin__r, rad_scale=K.HLR_pix)
        
        # Saving for later :D
        ALL._Mcor__g.append(Mcor__z)
        ALL._McorSD__g.append(McorSD__z)
        ALL.Mcor_GAL__g[iGal] = Mcor_GAL
        ALL.McorSD_GAL__g[iGal] = McorSD_GAL
        ALL.McorSD__rg[:, iGal] = McorSD__r
        ALL.at_flux_GAL__g[iGal] = at_flux_GAL
        ALL._at_flux__g.append(at_flux__z)
        ALL._at_mass__g.append(at_mass__z)
        ALL.at_flux__rg[:, iGal] = at_flux__r
        ALL.at_mass__rg[:, iGal] = at_mass__r
        
        # Composition by StarForming time scale
        for iT, tSF in enumerate(tSF__T):
            Mcor__z = np.ma.masked_array(K.Mcor__z)
            McorSD__z = np.ma.masked_array(K.Mcor__z / K.zoneArea_pc2)
            tau_V__z = np.ma.masked_array(K.tau_V__z)
            at_flux__z = np.ma.masked_array(K.at_flux__z)
            at_mass__z = np.ma.masked_array(K.at_mass__z)
            
            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            aux = calc_SFR(K, tSF)
            SFR__z = np.ma.masked_array(aux[0])
            SFRSD__z = np.ma.masked_array(aux[1])
            maskNotOk__z = np.zeros((K.N_zone), dtype=np.bool)
            
            if args.filter_residual is True:
                maskNotOk__z |= ~(K.filterResidual(w2 = 4600))
                C.debug_var(args.debug, pref = '    >>>',
                    tSF = '%.2fMyr' % (tSF / 1e6),
                    NnotOkResidual = (~(K.filterResidual())).sum()
                )

            if args.minpopx >= 0.:
                # Compute xOk "raw" image
                maskNotOk__z |= (x_Y__z < args.minpopx) 
                maskNotOk__z |= (tau_V__z < args.mintauv)
                 
                C.debug_var(args.debug, pref = '    >>>',
                    tSF = '%.2fMyr' % (tSF / 1e6),
                    NnotOkxY = (x_Y__z < args.minpopx).sum(),
                    NnotOktauV = (tau_V__z < args.mintauv).sum(),
                    NnotOk = maskNotOk__z.sum(),
                )

            tau_V__z[maskNotOk__z] = np.ma.masked
            SFR__z[maskNotOk__z] = np.ma.masked
            SFRSD__z[maskNotOk__z] = np.ma.masked
            Mcor__z[maskNotOk__z] = np.ma.masked
            McorSD__z[maskNotOk__z] = np.ma.masked
            at_flux__z[maskNotOk__z] = np.ma.masked
            at_mass__z[maskNotOk__z] = np.ma.masked
            integrated_SFR = SFR__z.sum()

            # Radial Profiles:
            x_Y__r = K.zoneToRad(x_Y__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            McorSD__r = K.zoneToRad(Mcor__z, Rbin__r, rad_scale=K.HLR_pix)
            aSFRSD__r = K.zoneToRad(SFRSD__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            tau_V__r = K.zoneToRad(tau_V__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            at_flux__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            at_mass__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale=K.HLR_pix, extensive=False, surface_density=False)
            at_flux_dezon__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale=K.HLR_pix)
            at_mass_dezon__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale=K.HLR_pix)
            Lobn__yx = K.zoneToYX(K.Lobn__z, extensive = False, surface_density = False)
            at_flux__yx = K.zoneToYX(at_flux__z, extensive = False, surface_density = False)
            McorSD__yx = K.zoneToYX(Mcor__z)
            at_mass__yx = K.zoneToYX(at_mass__z, extensive = False, surface_density = False)
            at_flux_wei__r = radialProfileWeighted(at_flux__yx, Lobn__yx, r_func=K.radialProfile, bin_r=Rbin__r, rad_scale=K.HLR_pix)
            at_mass_wei__r = radialProfileWeighted(at_mass__yx, McorSD__yx, r_func=K.radialProfile, bin_r=Rbin__r, rad_scale=K.HLR_pix)

            # Saving for later :D                
            ALL._x_Y__Tg[iT].append(x_Y__z)
            ALL._tau_V__Tg[iT].append(tau_V__z.data)
            ALL._tau_V_mask__Tg[iT].append(tau_V__z.mask)
            ALL._Mcor__Tg[iT].append(Mcor__z)
            ALL._McorSD__Tg[iT].append(McorSD__z)
            ALL._SFR__Tg[iT].append(SFR__z.data)
            ALL._SFR_mask__Tg[iT].append(SFR__z.mask)
            ALL._SFRSD__Tg[iT].append(SFRSD__z.data)
            ALL._SFRSD_mask__Tg[iT].append(SFRSD__z.mask)
            ALL._at_flux__Tg[iT].append(at_flux__z)
            ALL._at_mass__Tg[iT].append(at_mass__z)
            ALL.integrated_x_Y__Tg[:, iGal] = integrated_x_Y
            ALL.integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL.integrated_SFRSD__Tg[iT, iGal] = integrated_SFR / K.zoneArea_pc2.sum()
            ALL.x_Y__Trg[iT, :, iGal] = x_Y__r
            ALL.McorSD__Trg[iT, :, iGal] = McorSD__r
            ALL.aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
            ALL.tau_V__Trg[iT, :, iGal] = tau_V__r
            ALL.at_flux__Trg[iT, :, iGal] = at_flux__r
            ALL.at_mass__Trg[iT, :, iGal] = at_mass__r
            ALL.at_flux_dezon__Trg[iT, :, iGal] = at_flux_dezon__r
            ALL.at_mass_dezon__Trg[iT, :, iGal] = at_mass_dezon__r
            ALL.at_flux_wei__Trg[iT, :, iGal] = at_flux_wei__r
            ALL.at_mass_wei__Trg[iT, :, iGal] = at_mass_wei__r
            
        for iU, tZ in enumerate(tZ__U):
            aux = calc_alogZ_Stuff(K, tZ, args.minpopx, Rbin__r)
            alogZ_mass__z = aux[0]
            alogZ_flux__z = aux[1]
            alogZ_mass_GAL = aux[2]
            alogZ_flux_GAL = aux[3]
            alogZ_mass__r = aux[4]
            alogZ_flux__r = aux[5]
            alogZ_mass_wei__r = aux[6]
            alogZ_flux_wei__r = aux[7]
            isOkFrac_GAL = aux[8]
            
            # Saving for later :D
            ALL._alogZ_mass__Ug[iU].append(alogZ_mass__z.data)
            ALL._alogZ_mass_mask__Ug[iU].append(alogZ_mass__z.mask)
            ALL._alogZ_flux__Ug[iU].append(alogZ_flux__z.data)
            ALL._alogZ_flux_mask__Ug[iU].append(alogZ_flux__z.mask)
            ALL.alogZ_mass_GAL__Ug[iU, iGal] = alogZ_mass_GAL
            ALL.alogZ_flux_GAL__Ug[iU, iGal] = alogZ_flux_GAL
            ALL.alogZ_mass__Urg[iU, :, iGal] = alogZ_mass__r
            ALL.alogZ_flux__Urg[iU, :, iGal] = alogZ_flux__r
            ALL.alogZ_mass_wei__Urg[iU, :, iGal] = alogZ_mass_wei__r
            ALL.alogZ_flux_wei__Urg[iU, :, iGal] = alogZ_flux_wei__r
        ####################################################
        ####################################################
        ####################################################    
        
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
                    pos = K.GP._dlcons[l]['pos']
                    sigma = K.GP._dlcons[l]['sigma']
                    snr = K.GP._dlcons[l]['SN']
                    if snr < 3.0: snr = 3.0
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

        maskOkBPT__z = np.ones((K.N_zone), dtype = np.bool)
        maskOkwhan__z = np.ones((K.N_zone), dtype = np.bool)
        if args.underS06:
            L = Lines()
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
            maskOkBPT__z = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        if args.whanSF:
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            WHa = K.EL.EW[i_Ha, :]
            maskOkwhan__z = np.bitwise_and(np.greater(WHa, 3.), np.less_equal(N2Ha, -0.4))
        ##########################

        ########## tau_V #########
        maskOkTauVNeb__z = np.ones((K.N_zone), dtype = np.bool)
        if args.mintauvneb >= 0:
            maskOkTauVNeb__z = (K.EL.tau_V_neb__z >= args.mintauvneb) & (K.EL.tau_V_neb_err__z <= args.maxtauvneberr)
        maskOkTauVNeb__z &= maskOkLine[Ha_central_wl]
        maskOkTauVNeb__z &= maskOkLine[Hb_central_wl]
        maskOkTauVNeb__z &= maskOkBPT__z
        maskOkTauVNeb__z &= maskOkwhan__z
        
        C.debug_var(args.debug, 
            NOkHb = maskOkLine[Hb_central_wl].sum(),
            NOkO3 = maskOkLine[O3_central_wl].sum(),
            NOkHa = maskOkLine[Ha_central_wl].sum(),
            NOkN2 = maskOkLine[N2_central_wl].sum(),
            NOkBPT = maskOkBPT__z.sum(),
            NOkminTauVNeb = (K.EL.tau_V_neb__z >= args.mintauvneb).sum(),
            NOkmaxTauVNebErr = (K.EL.tau_V_neb_err__z <= args.maxtauvneberr).sum(),
            NOkTauVNeb = maskOkTauVNeb__z.sum(),
            NOkHaHb = (maskOkLine[Ha_central_wl] & maskOkLine[Hb_central_wl]).sum(),
            NOkO3N2 = (maskOkLine[O3_central_wl] & maskOkLine[N2_central_wl]).sum(),
            NOkLines = (maskOkLines__z).sum(),
        )
        
        tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = ~maskOkTauVNeb__z)
        tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = ~maskOkTauVNeb__z)
        tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
        tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
        tau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
        tau_V_neb_err__r = K.radialProfile(tau_V_neb_err__yx, Rbin__r, rad_scale = K.HLR_pix)

        # Saving for later :D
        ALL._tau_V_neb__g.append(tau_V_neb__z.data)
        ALL._tau_V_neb_err__g.append(tau_V_neb_err__z.data)
        ALL._tau_V_neb_mask__g.append(tau_V_neb__z.mask)
        ALL.tau_V_neb__rg[:, iGal] = tau_V_neb__r
        ALL.tau_V_neb_err__rg[:, iGal] = tau_V_neb_err__r
        ALL.integrated_tau_V_neb__g[iGal] = K.EL.integrated_tau_V_neb
        ALL.integrated_tau_V_neb_err__g[iGal] = K.EL.integrated_tau_V_neb_err
        ##########################

        ######### Z_neb ##########
        maskOkZNeb__z = (maskOkTauVNeb__z & maskOkLines__z & maskOkBPT__z)
        logZ_neb_S06__z = np.ma.masked_array(K.EL.logZ_neb_S06__z, mask = ~maskOkZNeb__z)
        logZ_neb_S06_err__z = np.ma.masked_array(K.EL.logZ_neb_S06_err__z, mask = ~maskOkZNeb__z)
        logZ_neb_S06__yx = K.zoneToYX(logZ_neb_S06__z, extensive = False)
        logZ_neb_S06_err__yx = K.zoneToYX(logZ_neb_S06_err__z, extensive = False)
        logZ_neb_S06__r = K.radialProfile(logZ_neb_S06__yx, Rbin__r, rad_scale = K.HLR_pix)
        logZ_neb_S06_err__r = K.radialProfile(logZ_neb_S06_err__yx, Rbin__r, rad_scale = K.HLR_pix)

        # Saving for later :D
        ALL._logZ_neb_S06__g.append(logZ_neb_S06__z.data)
        ALL._logZ_neb_S06_mask__g.append(logZ_neb_S06__z.mask)
        ALL._logZ_neb_S06_err__g.append(logZ_neb_S06_err__z.data)
        ALL.logZ_neb_S06__rg[:, iGal] = logZ_neb_S06__r
        ALL.logZ_neb_S06_err__rg[:, iGal] = logZ_neb_S06_err__r
        ALL.integrated_logZ_neb_S06__g[iGal] = K.EL.integrated_logZ_neb_S06
        ALL.integrated_logZ_neb_S06_err__g[iGal] = K.EL.integrated_logZ_neb_S06_err
        ##########################

        ########### EW ###########
        EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        EW_Hb__z = np.ma.masked_array(K.EL.EW[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
        EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
        EW_Ha__r = K.radialProfile(EW_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        EW_Hb__r = K.radialProfile(EW_Hb__yx, Rbin__r, rad_scale = K.HLR_pix)
        
        # Saving for later :D        
        ALL._EW_Ha__g.append(EW_Ha__z.data)
        ALL._EW_Ha_mask__g.append(EW_Ha__z.mask)
        ALL._EW_Hb__g.append(EW_Hb__z.data)
        ALL._EW_Hb_mask__g.append(EW_Hb__z.mask)
        ALL.EW_Ha__rg[:, iGal] = EW_Ha__r
        ALL.EW_Hb__rg[:, iGal] = EW_Hb__r
        ALL.integrated_EW_Ha__g[iGal] = K.EL.integrated_EW[i_Ha]
        ALL.integrated_EW_Hb__g[iGal] = K.EL.integrated_EW[i_Hb]
        ##########################

        #### intrinsic Ha Lum ####
        F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun        
        L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = ~(maskOkLine[Ha_central_wl]))
        L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = ~(maskOkLine[Hb_central_wl]))
        L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z
        # L_int_Ha__Lz intrinsic Ha luminosity 
        eHa = np.ma.exp(K.EL._qCCM['6563'] * tau_V_neb__z)
        # For the zones where I don't have values for tau_V_neb I don't correct the Lum_Ha
        L_int_Ha__z = np.where(maskOkTauVNeb__z, L_obs_Ha__z * eHa, L_obs_Ha__z)
        L_int_Ha__z = np.ma.masked_array(L_int_Ha__z, mask = ~maskOkTauVNeb__z)
        integrated_eHa = np.ma.exp(K.EL._qCCM['6563'] * K.EL.integrated_tau_V_neb)
        integrated_L_obs__L = K.EL._F_to_L(K.EL.integrated_flux) / L_sun
        integrated_L_int_Ha = integrated_L_obs__L[i_Ha] * integrated_eHa
        # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
        q = K.EL._qCCM['6563'] / (K.EL._qCCM['4861'] - K.EL._qCCM['6563'])
        a = L_obs_Ha_err__z
        b = q * L_obs_HaHb__z * L_obs_Hb_err__z
        L_int_Ha_err__z = np.where(maskOkTauVNeb__z == True, L_obs_Ha_err__z, eHa * np.sqrt(a ** 2.0 + b ** 2.0))
        L_int_Ha_err__z = np.ma.masked_array(L_int_Ha_err__z, mask = ~maskOkTauVNeb__z)
        
        # Saving for later :D
        ALL._F_obs_Ha__g.append(F_obs_Ha__z.data)
        ALL._F_obs_Ha_mask__g.append(F_obs_Ha__z.mask)
        ALL._L_obs_Ha__g.append(L_obs_Ha__z.data)
        ALL._L_obs_Ha_mask__g.append(L_obs_Ha__z.mask)
        ALL._L_obs_Ha_err__g.append(L_obs_Ha_err__z.data)
        ALL._L_int_Ha__g.append(L_int_Ha__z.data)
        ALL._L_int_Ha_mask__g.append(L_int_Ha__z.mask)
        ALL._L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        ALL.integrated_L_obs_Ha__g[iGal] = integrated_L_obs__L[i_Ha]
        ALL.integrated_L_int_Ha__g[iGal] = integrated_L_int_Ha
        ##########################

        #### SFR and SigmaSFR ####
        # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter
        SFR_Ha__z = np.ma.masked_array(3.13 * L_int_Ha__z.data / (1.e8), mask = L_int_Ha__z.mask)
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        integrated_SFR_Ha = 3.13 * integrated_L_int_Ha / (1.e8)
        integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()
        SFRSD_Ha__yx = K.zoneToYX(SFRSD_Ha__z, extensive = False)
        aSFRSD_Ha__r = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
        SFRSD_Ha_masked__yx = K.zoneToYX(np.ma.masked_array(SFRSD_Ha__z, mask = ~maskOkTauVNeb__z), extensive = False)
        aSFRSD_Ha_masked__r = K.radialProfile(SFRSD_Ha_masked__yx, Rbin__r, rad_scale = K.HLR_pix)

        # Saving for later :D
        ALL.aSFRSD_Ha__rg[:, iGal] = aSFRSD_Ha__r
        ALL.aSFRSD_Ha_masked__rg[:, iGal] = aSFRSD_Ha_masked__r
        ALL.integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha
        ALL.integrated_SFRSD_Ha__g[iGal] = integrated_SFRSD_Ha
        ALL._SFR_Ha__g.append(SFR_Ha__z.data)
        ALL._SFR_Ha_mask__g.append(SFR_Ha__z.mask)
        ALL._SFRSD_Ha__g.append(SFRSD_Ha__z.data)
        ALL._SFRSD_Ha_mask__g.append(SFRSD_Ha__z.mask)
        ####################################################
        ####################################################
        ####################################################

        if args.gasprop:
            ####################################################
            # GasProp Ruben ####################################
            ####################################################
            # Values in GasProp could be NaN. This values will be masked at 
            # ALL.stiack_zones_data() with mask = np.isnan(vect)
            chb_in__z = K.GP.REDDENING.chb_in
            c_Ha_Hb__z = K.GP.REDDENING.c_Ha_Hb
            # O_HIICHIM may have zeros.
            O_HIICHIM__z = np.where(K.GP.EMPAB.O_HIICHIM == 0., np.nan, K.GP.EMPAB.O_HIICHIM)
            O_O3N2_M13__z = K.GP.EMPAB.O_O3N2_M13
            O_O3N2_PP04__z = K.GP.EMPAB.O_O3N2_PP04
            O_direct_O_23__z = K.GP.ELEMAB.O_direct_O_23
            _O_HIICHIM__z = np.ma.masked_array(O_HIICHIM__z, mask = np.isnan(O_HIICHIM__z)) 
            _O_O3N2_M13__z = np.ma.masked_array(O_O3N2_M13__z, mask = np.isnan(O_O3N2_M13__z))
            _O_O3N2_PP04__z = np.ma.masked_array(O_O3N2_PP04__z, mask = np.isnan(O_O3N2_PP04__z))
            _O_direct_O_23__z = np.ma.masked_array(O_direct_O_23__z, mask = np.isnan(O_direct_O_23__z))
            O_HIICHIM__yx = K.zoneToYX(_O_HIICHIM__z, extensive = False)
            O_HIICHIM__r = K.radialProfile(O_HIICHIM__yx, Rbin__r, rad_scale = K.HLR_pix)
            O_O3N2_M13__yx = K.zoneToYX(_O_O3N2_M13__z, extensive = False)
            O_O3N2_M13__r = K.radialProfile(O_O3N2_M13__yx, Rbin__r, rad_scale = K.HLR_pix)
            O_O3N2_PP04__yx = K.zoneToYX(_O_O3N2_PP04__z, extensive = False)
            O_O3N2_PP04__r = K.radialProfile(O_O3N2_PP04__yx, Rbin__r, rad_scale = K.HLR_pix)
            O_direct_O_23__yx = K.zoneToYX(_O_direct_O_23__z, extensive = False)
            O_direct_O_23__r = K.radialProfile(O_direct_O_23__yx, Rbin__r, rad_scale = K.HLR_pix)
            integrated_chb_in = K.GP.REDDENING.integrated_chb_in
            integrated_c_Ha_Hb = K.GP.REDDENING.integrated_c_Ha_Hb
            integrated_O_HIICHIM = K.GP.EMPAB.integrated_O_HIICHIM
            integrated_O_O3N2_M13 = K.GP.EMPAB.integrated_O_O3N2_M13
            integrated_O_O3N2_PP04 = K.GP.EMPAB.integrated_O_O3N2_PP04
            integrated_O_direct_O_23 = K.GP.ELEMAB.integrated_O_direct_O_23
            
            # Saving for later :D
            ALL._chb_in__g.append(chb_in__z)
            ALL._c_Ha_Hb__g.append(c_Ha_Hb__z)
            ALL._O_HIICHIM__g.append(O_HIICHIM__z)
            ALL._O_O3N2_M13__g.append(O_O3N2_M13__z)
            ALL._O_O3N2_PP04__g.append(O_O3N2_PP04__z)
            ALL._O_direct_O_23__g.append(O_direct_O_23__z)
            ALL.integrated_chb_in__g[iGal] = integrated_chb_in
            ALL.integrated_c_Ha_Hb__g[iGal] = integrated_c_Ha_Hb
            ALL.integrated_O_HIICHIM__g[iGal] = integrated_O_HIICHIM
            ALL.integrated_O_O3N2_M13__g[iGal] = integrated_O_O3N2_M13
            ALL.integrated_O_O3N2_PP04__g[iGal] = integrated_O_O3N2_PP04
            ALL.integrated_O_direct_O_23__g[iGal] = integrated_O_direct_O_23
            ALL.O_HIICHIM__rg[:, iGal] = O_HIICHIM__r
            ALL.O_O3N2_M13__rg[:, iGal] = O_O3N2_M13__r
            ALL.O_O3N2_PP04__rg[:, iGal] = O_O3N2_PP04__r 
            ALL.O_direct_O_23__rg[:, iGal] = O_direct_O_23__r

        if args.gasprop is True:
            K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)

    t_init_stack = time.clock()        
    ALL.stack_zones_data()
    ALL.integrated_mask_nan()
    print 'time stack: %.2f' % (time.clock() - t_init_stack)
    
    if args.hdf5 != None:
        t_init_hdf5 = time.clock()        
        import h5py
        filename = args.hdf5
        h5 = h5py.File(filename, 'w')
        D = ALL.create_dict_h5()
        D['data/RbinIni'] = args.rbinini
        D['data/RbinFin'] = args.rbinfin
        D['data/RbinStep'] = args.rbinstep
        D['data/Rbin__r'] = Rbin__r
        D['data/RbinCenter__r'] = RbinCenter__r
        D['data/NRbins'] = NRbins
        D['data/RColor'] = RColor
        D['data/RRange'] = RRange
        D['data/tSF__T'] = tSF__T
        D['data/tZ__U'] = tZ__U
        D['data/xOkMin'] = args.minpopx
        D['data/tauVOkMin'] = args.mintauv
        D['data/tauVNebOkMin'] = args.mintauvneb
        D['data/tauVNebErrMax'] = args.maxtauvneberr
        for k in D.keys():
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
        h5.close()
        print 'time hdf5: %.2f' % (time.clock() - t_init_hdf5)
        
    print 'total time: %.2f' % (time.clock() - t_init_prog)

    