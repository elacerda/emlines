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
from CALIFAUtils.lines import Lines
from CALIFAUtils.globals import califa_work_dir
from CALIFAUtils.scripts import get_morfologia
from CALIFAUtils.scripts import ALLGals
from CALIFAUtils.scripts import calc_SFR
from CALIFAUtils.scripts import calc_xY
from CALIFAUtils.scripts import calc_alogZ_Stuff
from CALIFAUtils.scripts import radialProfileWeighted
from CALIFAUtils.scripts import loop_cubes
from CALIFAUtils.scripts import sort_gals

def parser_args():
    default_args = {
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
        'maxiTSF' :-1,
        'gals_filename' : califa_work_dir + 'listOf300GalPrefixes.txt'
    }
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--spiral', '-S',
                        action = 'store_true',
                        default = default_args['spiral'])
    parser.add_argument('--underS06',
                        action = 'store_true',
                        default = default_args['underS06'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default_args['weiradprof'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
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
    parser.add_argument('--gals_filename', '-L',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['gals_filename'])
    return parser.parse_args()

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    t_init_prog = time.clock()

    args = parser_args()
    print_args(args)    
    
    imgDir = califa_work_dir + 'images/'
    
    Zsun = 0.019

    Rbin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins = len(RbinCenter__r)
    RColor = [ 'r', 'y', 'b', 'k' ]
    RRange = [  .5, 1., 1.5, 2.  ]
    
    gals, _ = sort_gals(args.gals_filename)
    maxGals = None
    if args.debug:
        maxGals = 10
    N_gals = len(gals)

    # SFR-time-scale array (index __T)
    base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
    if args.maxiTSF > 0:
        tSF__T = np.asarray(base.ageBase[0:args.maxiTSF + 1])
    else:
        tSF__T = np.asarray(base.ageBase)
    N_T = len(tSF__T)

    tZ__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
    N_U = len(tZ__U)

    ALL = ALLGals(N_gals, NRbins, N_T, N_U)
    
    # automatically read PyCASSO and EmLines data cubes.
    for iGal, K in loop_cubes(gals.tolist(), imax = maxGals, EL = True):
        t_init_gal = time.clock()
        califaID = gals[iGal] 

        ALL.califaIDs[iGal] = califaID
        
        if K == None:
            ALL.mask_gal(iGal)
            print '<<< %s galaxy: miss files' % califaID 
            continue

        ALL.N_zones__g[iGal] = K.N_zone
        
        if K.EL == None:
            ALL.mask_gal(iGal)
            print '<<< %s galaxy: miss EmLines files' % califaID 
            continue

        # Problem in FITS file
        if K.EL.flux[0, :].sum() == 0.:
            ALL.mask_gal(iGal)
            print '<<< %s EmLines FITS problem' % califaID
            continue
        
        tipos, tipo, tipo_m, tipo_p = get_morfologia(califaID)
        
        # Only spiral
        if args.spiral and tipo <= 8: 
            ALL.mask_gal(iGal)
            print '<<< %s galaxy: is not a spiral (type: %d)' % (califaID, tipo) 
            continue
        
        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        print '>>> Doing' , iGal , califaID , 'hubtyp=', tipo, '|  Nzones=' , K.N_zone
        
        # zone distance in HLR
        zoneDistHLR = np.sqrt((K.zonePos['x'] - K.x0) ** 2. + (K.zonePos['y'] - K.y0) ** 2.) / K.HLR_pix
        ALL._dist_zone__g.append(zoneDistHLR)
        
        ALL.ba_PyCASSO_GAL__g[iGal] = ba
        ALL.ba_GAL__g[iGal] = np.float(K.masterListData['ba'])
        ALL.Mr_GAL__g[iGal] = np.float(K.masterListData['Mr'])
        ALL.ur_GAL__g[iGal] = np.float(K.masterListData['u-r'])
        ALL.morfType_GAL__g[iGal] = tipo
        
        ####################################################
        ####### STARLIGHT ##################################
        ####################################################
        AVtotauV = 1. / (np.log10(np.exp(1)) / 0.4)
        ALL.integrated_tau_V__g[iGal] = K.integrated_keywords['A_V'] * AVtotauV 
        
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        ALL._Mcor__g.append(K.Mcor__z)
        ALL._McorSD__g.append(K.Mcor__z / K.zoneArea_pc2)
        ALL.Mcor_GAL__g[iGal] = K.Mcor_tot.sum()
        ALL.McorSD_GAL__g[iGal] = K.McorSD__yx.mean()
        
        ALL.McorSD__rg[:, iGal] = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)

        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        ALL.at_flux_GAL__g[iGal] = numerator__z.sum() / denominator__z.sum()
        
        ALL._at_flux__g.append(K.at_flux__z)
        ALL._at_mass__g.append(K.at_flux__z)
        
        ALL.at_flux__rg[:, iGal] = K.radialProfile(K.at_flux__yx, Rbin__r, rad_scale = K.HLR_pix)
        ALL.at_mass__rg[:, iGal] = K.radialProfile(K.at_mass__yx, Rbin__r, rad_scale = K.HLR_pix)
        
        for iT, tSF in enumerate(tSF__T):
            Mcor__z = np.ma.masked_array(K.Mcor__z)
            McorSD__z = np.ma.masked_array(K.Mcor__z / K.zoneArea_pc2)
            tau_V__z = np.ma.masked_array(K.tau_V__z)
            at_flux__z = np.ma.masked_array(K.at_flux__z)
            at_mass__z = np.ma.masked_array(K.at_mass__z)
            
            x_Y__z = calc_xY(K, tSF)
            aux = calc_SFR(K, tSF)
            SFR__z = np.ma.masked_array(aux[0])
            SFRSD__z = np.ma.masked_array(aux[1])

            if args.minpopx >= 0.:
                # Compute xOk "raw" image
                maskNotOk__z = (x_Y__z < args.minpopx) | (tau_V__z < args.mintauv) 
                 
                tau_V__z[maskNotOk__z] = np.ma.masked
                SFR__z[maskNotOk__z] = np.ma.masked
                SFRSD__z[maskNotOk__z] = np.ma.masked
                Mcor__z[maskNotOk__z] = np.ma.masked
                McorSD__z[maskNotOk__z] = np.ma.masked
                at_flux__z[maskNotOk__z] = np.ma.masked
                at_mass__z[maskNotOk__z] = np.ma.masked
                
            tau_V__yx = K.zoneToYX(tau_V__z, extensive = False, surface_density = False)
            aSFRSD__yx = K.zoneToYX(SFRSD__z, extensive = False, surface_density = False)
            McorSD__yx = K.zoneToYX(Mcor__z, extensive = True)
            at_flux__yx = K.zoneToYX(at_flux__z, extensive = False, surface_density = False)
            at_mass__yx = K.zoneToYX(at_mass__z, extensive = False, surface_density = False)
            at_flux_dezon__yx = K.zoneToYX(at_flux__z, extensive = True)
            at_mass_dezon__yx = K.zoneToYX(at_mass__z, extensive = True)
            
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

            integrated_SFR = SFR__z.sum()
            ALL.integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL.integrated_SFRSD__Tg[iT, iGal] = integrated_SFR / K.zoneArea_pc2.sum()

            ALL.McorSD__Trg[iT, :, iGal] = K.radialProfile(McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.aSFRSD__Trg[iT, :, iGal] = K.radialProfile(aSFRSD__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.tau_V__Trg[iT, :, iGal] = K.radialProfile(tau_V__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.at_flux__Trg[iT, :, iGal] = K.radialProfile(at_flux__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.at_mass__Trg[iT, :, iGal] = K.radialProfile(at_mass__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.at_flux_dezon__Trg[iT, :, iGal] = K.radialProfile(at_flux_dezon__yx, Rbin__r, rad_scale = K.HLR_pix)
            ALL.at_mass_dezon__Trg[iT, :, iGal] = K.radialProfile(at_mass_dezon__yx, Rbin__r, rad_scale = K.HLR_pix)
            Lobn__yx = K.zoneToYX(K.Lobn__z, extensive = True)
            ALL.at_flux_wei__Trg[iT, :, iGal] = radialProfileWeighted(at_flux__yx, Lobn__yx, r_func = K.radialProfile, bin_r = Rbin__r, rad_scale = K.HLR_pix)
            ALL.at_mass_wei__Trg[iT, :, iGal] = radialProfileWeighted(at_mass__yx, McorSD__yx, r_func = K.radialProfile, bin_r = Rbin__r, rad_scale = K.HLR_pix)
            
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

        ###### MASK EmLines ######
        # minimum value of f_lz / err_f_lz
        minSNR = 3.
        
        i_Hb = K.EL.lines.index('4861')
        i_O3 = K.EL.lines.index('5007')
        i_Ha = K.EL.lines.index('6563')
        i_N2 = K.EL.lines.index('6583')
        Ha = K.EL.flux[i_Ha, :]
        Hb = K.EL.flux[i_Hb, :]
        O3 = K.EL.flux[i_O3, :]
        N2 = K.EL.flux[i_N2, :]
        eHa = K.EL.eflux[i_Ha, :]
        eHb = K.EL.eflux[i_Hb, :]
        eO3 = K.EL.eflux[i_O3, :]
        eN2 = K.EL.eflux[i_N2, :]
        
        HbOk = np.array((Hb / eHb) >= minSNR, dtype = np.bool)
        O3Ok = np.array((O3 / eO3) >= minSNR, dtype = np.bool)
        HaOk = np.array((Ha / eHa) >= minSNR, dtype = np.bool)
        N2Ok = np.array((N2 / eN2) >= minSNR, dtype = np.bool)
        
        maskLinesSNR__z = HbOk & O3Ok & HaOk & N2Ok
        maskOkFlux__z = (Ha >= 0) & (Hb >= 0) & (O3 >= 0) & (N2 >= 0)
        maskOkBPT__z = np.ones((K.N_zone), dtype = np.bool)
        
        if args.underS06:
            L = Lines()
            N2Ha = np.ma.log10(K.EL.N2_obs__z / K.EL.Ha_obs__z)
            O3Hb = np.ma.log10(K.EL.O3_obs__z / K.EL.Hb_obs__z)
            maskOkBPT__z = L.maskBelowlinebpt('S06', N2Ha, O3Hb)
        ##########################

        ########## tau_V #########
        maskOkTauVNeb__z = np.ones((K.N_zone), dtype = np.bool)
        
        if args.mintauvneb >= 0:
            maskOkTauVNeb__z = (K.EL.tau_V_neb__z >= args.mintauvneb) & (K.EL.tau_V_neb_err__z <= args.maxtauvneberr)

        maskOkNeb__z = (maskOkTauVNeb__z & maskLinesSNR__z & maskOkFlux__z & maskOkBPT__z)
        
        tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = ~maskOkNeb__z)
        tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = ~maskOkNeb__z)
        tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
        tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z, extensive = False)
        tau_V_neb__r = K.radialProfile(tau_V_neb__yx, Rbin__r, rad_scale = K.HLR_pix)

        ALL._tau_V_neb__g.append(tau_V_neb__z.data)
        ALL._tau_V_neb_err__g.append(tau_V_neb_err__z.data)
        ALL._tau_V_neb_mask__g.append(tau_V_neb__z.mask)
        
        ALL.tau_V_neb__rg[:, iGal] = tau_V_neb__r

        ALL.integrated_tau_V_neb__g[iGal] = K.EL.integrated_tau_V_neb
        ##########################

        ######### Z_neb ##########
        logZ_neb_S06__z = np.ma.masked_array(K.EL.logZ_neb_S06__z, mask = ~maskOkNeb__z)
        logZ_neb_S06__yx = K.zoneToYX(logZ_neb_S06__z, extensive = False)
        logZ_neb_S06__r = K.radialProfile(logZ_neb_S06__yx, Rbin__r, rad_scale = K.HLR_pix)
        logZ_neb_S06_err__z = np.ma.masked_array(K.EL.logZ_neb_S06_err__z, mask = ~maskOkNeb__z)
        logZ_neb_S06_err__yx = K.zoneToYX(logZ_neb_S06_err__z, extensive = False)

        ALL._logZ_neb_S06__g.append(logZ_neb_S06__z.data)
        ALL._logZ_neb_S06_mask__g.append(logZ_neb_S06__z.mask)
        ALL._logZ_neb_S06_err__g.append(logZ_neb_S06_err__z.data)
        
        ALL.logZ_neb_S06__rg[:, iGal] = logZ_neb_S06__r
        
        ALL.integrated_logZ_neb_S06__g[iGal] = K.EL.integrated_logZ_neb_S06
        ##########################

        ########### EW ###########
        EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = ~maskOkNeb__z)
        EW_Hb__z = np.ma.masked_array(K.EL.EW[i_Hb, :], mask = ~maskOkNeb__z)
        
        ALL._EW_Ha__g.append(EW_Ha__z)
        ALL._EW_Hb__g.append(EW_Hb__z)
        ##########################

        #### intrinsic Ha Lum ####
        F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = ~maskOkNeb__z)
        
        ALL._F_obs_Ha__g.append(F_obs_Ha__z)
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun
        
        # L_int_Ha__Lz intrinsic Ha luminosity 
        q = K.EL._qCCM['6563'] / (K.EL._qCCM['4861'] - K.EL._qCCM['6563'])
        eHa = np.ma.exp(K.EL._qCCM['6563'] * tau_V_neb__z)
        
        L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = ~maskOkNeb__z)
        L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = ~maskOkNeb__z)
        L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = ~maskOkNeb__z)
        L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = ~maskOkNeb__z)
        
        L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z 
        L_int_Ha__z = L_obs_Ha__z * eHa
        
        integrated_eHa = np.ma.exp(K.EL._qCCM['6563'] * K.EL.integrated_tau_V_neb)
        integrated_L_obs_Ha__L = K.EL._F_to_L(K.EL.integrated_flux) / L_sun
        integrated_L_int_Ha = integrated_L_obs_Ha__L[i_Ha] * integrated_eHa
        
        # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
        a = L_obs_Ha_err__z
        b = q * L_obs_HaHb__z * L_obs_Hb_err__z
        L_int_Ha_err__z = eHa * np.sqrt(a ** 2.0 + b ** 2.0)
        
        ALL._L_int_Ha__g.append(L_int_Ha__z.data)
        ALL._L_int_Ha_mask__g.append(L_int_Ha__z.mask)
        ALL._L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        
        ALL.integrated_L_int_Ha__g[iGal] = integrated_L_int_Ha
        ##########################

        #### SFR and SigmaSFR ####
        # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter        
        SFR_Ha__z = 3.13 * L_int_Ha__z / (1.e8)
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        
        integrated_SFR_Ha = 3.13 * integrated_L_int_Ha / (1.e8)
        integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()
        
        SFRSD_Ha__yx = K.zoneToYX(SFRSD_Ha__z, extensive = False)
        aSFRSD_Ha__r = K.radialProfile(SFRSD_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)

        ALL.aSFRSD_Ha__rg[:, iGal] = aSFRSD_Ha__r
        
        ALL.integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha
        ALL.integrated_SFRSD_Ha__g[iGal] = integrated_SFRSD_Ha
        
        ALL._SFR_Ha__g.append(SFR_Ha__z)
        ALL._SFRSD_Ha__g.append(SFRSD_Ha__z)
        ####################################################
        ####################################################
        ####################################################
        
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)
        K.close()

    t_init_stack = time.clock()        
    ALL.stack_zones_data()
    print 'time stack: %.2f' % (time.clock() - t_init_stack)
    
    if args.hdf5 != None:
        t_init_hdf5 = time.clock()        
        import h5py
        filename = args.hdf5
        h5 = h5py.File(filename, 'w')
        D = ALL.create_dict_h5()
        D['/data/RbinIni'] = args.rbinini
        D['/data/RbinFin'] = args.rbinfin
        D['/data/RbinStep'] = args.rbinstep
        D['/data/Rbin__r'] = Rbin__r
        D['/data/RbinCenter__r'] = RbinCenter__r
        D['/data/NRbins'] = NRbins
        D['/data/RColor'] = RColor
        D['/data/RRange'] = RRange
        D['/data/tSF__T'] = tSF__T
        D['/data/tZ__U'] = tZ__U
        D['/data/xOkMin'] = args.minpopx
        D['/data/tauVOkMin'] = args.mintauv
        D['/data/tauVNebOkMin'] = args.mintauvneb
        D['/data/tauVNebErrMax'] = args.maxtauvneberr
        for k in D.keys():
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
        h5.close()
        print 'time hdf5: %.2f' % (time.clock() - t_init_hdf5)
        
    print 'total time: %.2f' % (time.clock() - t_init_prog)
