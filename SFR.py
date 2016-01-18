#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import sys
import time
import numpy as np
import argparse as ap
import CALIFAUtils as C
from CALIFAUtils.objects import ALLGals
from CALIFAUtils.scripts import calc_xY
from CALIFAUtils.scripts import calc_SFR
from pystarlight.util import redenninglaws
from pystarlight.util.constants import L_sun
from pystarlight.util.base import StarlightBase
from CALIFAUtils.scripts import calc_alogZ_Stuff
from CALIFAUtils.scripts import create_masks_gal 

def parser_args(args_str):
    paths = C.paths
    default_args = {
        'hdf5' : None,
        'debug' : False,
        'whanSF' : None,
        'spiral' : False,
        'rgbcuts' : False,
        'gasprop' : False,
        'underS06' : False,
        'weiradprof' : False,
        'nolinecuts' : False,
        'filter_residual' : False,
        'rbinini' : 0.0,
        'rbinfin' : 2.0,
        'rbinstep' : 0.1,
        'minpopx' : np.finfo(np.float_).min,
        'mintauv' : np.finfo(np.float_).min,
        'minEWHb' : np.finfo(np.float_).min,
        'mintauvneb' : np.finfo(np.float_).min,
        'maxtauvneberr' : np.finfo(np.float_).max,
        'v_run' : -1,
        'minSNR' : 3,
        'minSNRHb' : 3,
        'gals_filename' : paths.califa_work_dir + 'listv20_q050.d15a.txt',
    }
    
    parser = ap.ArgumentParser(description = '%s' % args_str)
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--gasprop', '-G',
                        action = 'store_true',
                        default = default_args['gasprop'])
    parser.add_argument('--nolinecuts' ,
                        action = 'store_true',
                        default = default_args['nolinecuts'])
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
                        metavar = 'INT',
                        type = int,
                        #action = 'store_true',
                        default = default_args['whanSF'])
    parser.add_argument('--rgbcuts',
                        action = 'store_true',
                        default = default_args['rgbcuts'])
    parser.add_argument('--minpopx',
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minpopx'])
    parser.add_argument('--minEWHb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minEWHb'])    
    parser.add_argument('--minSNR',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minSNR'])    
    parser.add_argument('--minSNRHb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minSNRHb'])    
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
    
    args = parser.parse_args()
    
    if args.nolinecuts:
        args.rgbcuts = False

    return args

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
    args = parser_args(sys.argv[0])
    C.debug_var(True, args = args.__dict__)    
    
    Zsun = 0.019

    # Creating radial bins.
    Rbin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins = len(RbinCenter__r)
    RColor = [ 'r', 'y', 'b', 'k' ]
    RRange = [  .5, 1., 1.5, 2.  ]
    Rbin_oneHLR = [1. - args.rbinstep, 1. + args.rbinstep]
    
    # Reading galaxies file,
    # gals, _ = C.sort_gals(gals = args.gals_filename, order = 1)
    gals = np.array(['K0195', 'K0201', 'K0272', 'K0136', 'K0197', 'K0017', 'K0708',
                 'K0090', 'K0631', 'K0128', 'K0780', 'K0067', 'K0160', 'K0923',
                 'K0035', 'K0098', 'K0859', 'K0816', 'K0815', 'K0092', 'K0806',
                 'K0602', 'K0018', 'K0210', 'K0782', 'K0051', 'K0279', 'K0903',
                 'K0127', 'K0864', 'K0846', 'K0900', 'K0341', 'K0004', 'K0044',
                 'K0829', 'K0171', 'K0589', 'K0911', 'K0588', 'K0814', 'K0781',
                 'K0835', 'K0881', 'K0112', 'K0387', 'K0845', 'K0076', 'K0612',
                 'K0888', 'K0068', 'K0851', 'K0832', 'K0893', 'K0318', 'K0101',
                 'K0840', 'K0870', 'K0705', 'K0138', 'K0139', 'K0633', 'K0037',
                 'K0391', 'K0046', 'K0121', 'K0744', 'K0162', 'K0080', 'K0339',
                 'K0912', 'K0047', 'K0860', 'K0908', 'K0093', 'K0844', 'K0787',
                 'K0103', 'K0703', 'K0055', 'K0085', 'K0059', 'K0072', 'K0118',
                 'K0173', 'K0134', 'K0916', 'K0919', 'K0050', 'K0281', 'K0917',
                 'K0822', 'K0865', 'K0087', 'K0875', 'K0170', 'K0826', 'K0063',
                 'K0562', 'K0096', 'K0479', 'K0099', 'K0613', 'K0778', 'K0189',
                 'K0874', 'K0607', 'K0119', 'K0883', 'K0858', 'K0592',             
    #gals = np.array([
                'K0673',
                'K0811', 'K0274', 'K0100', 'K0091', 'K0867', 'K0061', 'K0105',
                'K0024', 'K0502', 'K0653', 'K0077', 'K0818', 'K0131', 'K0872',
                'K0057', 'K0838', 'K0111', 'K0186', 'K0174', 'K0634', 'K0933',
                'K0036', 'K0049', 'K0936', 'K0163', 'K0032', 'K0220', 'K0314',
                'K0889', 'K0132', 'K0026', 'K0169', 'K0029', 'K0078', 'K0075',
                'K0066', 'K0319', 'K0135', 'K0902', 'K0083', 'K0020', 'K0386',
                'K0853', 'K0007', 'K0863', 'K0809', 'K0123', 'K0156', 'K0925',
                'K0833', 'K0850', 'K0894', 'K0168', 'K0194', 'K0038', 'K0651',
                'K0219', 'K0791', 'K0381', 'K0062', 'K0886', 'K0364', 'K0624',
                'K0113', 'K0177', 'K0740', 'K0663', 'K0311', 'K0932', 'K0019',
                'K0783', 'K0910', 'K0797', 'K0664', 'K0185', 'K0518', 'K0924',
                'K0414', 'K0013', 'K0804', 'K0326', 'K0615', 'K0672', 'K0837',
                'K0914', 'K0856', 'K0848', 'K0569', 'K0676', 'K0915', 'K0065',
                'K0842', 'K0146', 'K0176', 'K0115', 'K0043', 'K0001', 'K0021',
                'K0824', 'K0192', 'K0023', 'K0097', 'K0890', 'K0188', 'K0868',
                'K0086', 'K0180', 'K0153', 'K0054', 'K0133', 'K0191', 'K0107',
                'K0307', 'K0610', 'K0025', 'K0164', 'K0073', 'K0807', 'K0854',
                'K0834', 'K0821', 'K0871', 'K0151', 'K0070', 'K0010', 'K0830',
                'K0102', 'K0278', 'K0869', 'K0280', 'K0009', 'K0684', 'K0665',
                'K0126', 'K0873', 'K0927', 'K0774', 'K0931', 'K0789', 'K0297',
                'K0934', 'K0383', 'K0825', 'K0630', 'K0715', 'K0152', 'K0768',
                'K0907', 'K0823', 'K0041', 'K0388', 'K0005', 'K0769', 'K0208',
                'K0831', 'K0748', 'K0196', 'K0130', 'K0659', 'K0611', 'K0798',
                'K0880', 'K0184', 'K0879', 'K0820', 'K0754', 'K0309', 'K0861',
                'K0137', 'K0436', 'K0116', 'K0110', 'K0088', 'K0500', 'K0275',
                'K0124', 'K0181', 'K0476', 'K0779', 'K0580', 'K0437', 'K0183',
                'K0901', 'K0042', 'K0094', 'K0147', 'K0108', 'K0052', 'K0898',
                'K0203', 'K0813', 'K0190', 'K0876', 'K0277', 'K0810', 'K0515',
                'K0008', 'K0489', 'K0109', 'K0608', 'K0929', 'K0714', 'K0904',
                'K0896', 'K0149', 'K0122', 'K0857', 'K0165', 'K0002', 'K0178',
                'K0028', 'K0140', 'K0764', 'K0887', 'K0849', 'K0141', 'K0652',
                'K0895', 'K0885', 'K0878', 'K0711', 'K0361', 'K0843', 'K0758',
                'K0713', 'K0069', 'K0852', 'K0089', 'K0697', 'K0827', 'K0095',
                'K0158', 'K0143', 'K0805', 'K0117', 'K0040', 'K0486', 'K0921',
                'K0355', 'K0841', 'K0817', 'K0260', 'K0166', 'K0603', 'K0884',
                'K0148', 'K0039', 'K0930', 'K0707', 'K0125', 'K0935', 'K0016',
                'K0891', 'K0828', 'K0609', 'K0548', 'K0129', 'K0836', 'K0012',
                'K0581', 'K0232', 'K0866', 'K0003', 'K0909', 'K0205', 'K0144',
                'K0906', 'K0159', 'K0030', 'K0187', 'K0157', 'K0053', 'K0775',
                'K0273', 'K0045', 'K0084', 'K0614', 'K0031', 'K0071', 'K0204',
                'K0081', 'K0033', 'K0034', 'K0862', 'K0231', 'K0226', 'K0182',
                'K0027', 'K0312', 'K0060', 'K0606', 'K0528', 'K0150', 'K0306',
                'K0657', 'K0161', 'K0353', 'K0749', 'K0209', 'K0058', 'K0014',
                'K0179', 'K0475', 'K0937'], 
                dtype='|S5')  
      
    N_gals = len(gals)
    maxGals = None
    if args.debug:
        maxGals = 10
        if N_gals > maxGals:
            N_gals = maxGals

    # SFR-time-scale array (index __T)
    base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
    #tSF__T = np.array([0.032 , 0.3 , 1.5, 14.2]) * 1.e9
    tSF__T = np.asarray(base.ageBase)
    #tSF__T = np.array([1, 3.2, 10, 100]) * 1e7
    N_T = len(tSF__T)

    # Z-time-scale array (index __U).
    tZ__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
    N_U = len(tZ__U)

    ALL = ALLGals(N_gals, NRbins, N_T, N_U)
    
    gals_not_masked = 0

    # automatically read PyCASSO and EmLines data cubes.
    for iGal, K in C.loop_cubes(gals.tolist(), imax = maxGals, 
                                EL = True, GP = args.gasprop, 
                                v_run = args.v_run, debug = 'True'):        
        t_init_gal = time.clock()
        califaID = gals[iGal] 
        
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
        my_type = C.my_morf(tipos)
        
        # Only spiral
        if args.spiral and my_type < 0:  # between Sa ... Sd
            ALL.mask_gal(iGal)
            if args.gasprop is True:
                K.GP.close()
            K.EL.close()
            K.close()
            print '<<< %s galaxy: is not a spiral (type: %f (%s))' % (califaID, my_type, tipos) 
            continue

        N_zone = K.N_zone

        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        # Saving for later :D
        ALL.N_zones__g[iGal] = N_zone
        ALL.ba_PyCASSO_GAL__g[iGal] = ba
        
        print '>>> Doing %d %s %s (%d)' % (iGal, califaID, tipos, my_type)
        
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
        ALL.morfType_GAL__g[iGal] = my_type

        # masks: more info. in create_masks() 
        mask__Tz, mask_syn__Tz, mask_eml__z, mask_popx__Tz, \
        mask_tau_V__z, mask_residual__z, mask_tau_V_neb__z, \
        mask_tau_V_neb_err__z, mask_EW_Hb__z, mask_whan__z, \
        mask_bpt__z, mask_lines_dict__Lz = create_masks_gal(K = K, tSF__T = tSF__T, args = args)
        
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
        McorSD__r = K.radialProfile(K.McorSD__yx, Rbin__r, rad_scale = K.HLR_pix)
        McorSD_oneHLR = K.radialProfile(K.McorSD__yx, Rbin_oneHLR, rad_scale = K.HLR_pix)
        at_flux__z = K.at_flux__z
        at_mass__z = K.at_mass__z
        HLR_pix = K.HLR_pix
        HMR_pix = K.getHalfRadius(K.McorSD__yx)
        parsecPerPixel = K.parsecPerPixel
        
        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        at_flux_GAL = numerator__z.sum() / denominator__z.sum()
        at_flux__r = K.radialProfile(K.at_flux__yx, Rbin__r, rad_scale = K.HLR_pix)
        at_mass__r = K.radialProfile(K.at_mass__yx, Rbin__r, rad_scale = K.HLR_pix)
        at_flux_oneHLR = K.radialProfile(K.at_flux__yx, Rbin_oneHLR, rad_scale = K.HLR_pix)
        at_mass_oneHLR = K.radialProfile(K.at_mass__yx, Rbin_oneHLR, rad_scale = K.HLR_pix)
        
        # Saving for later :D
        ALL._Mcor__g.append(Mcor__z)
        ALL._McorSD__g.append(McorSD__z)
        ALL.Mcor_GAL__g[iGal] = Mcor_GAL
        ALL.McorSD_GAL__g[iGal] = McorSD_GAL
        ALL.McorSD__rg[:, iGal] = McorSD__r
        ALL.McorSD_oneHLR__g[iGal] = McorSD_oneHLR
        ALL.HLR_pix_GAL__g[iGal] = HLR_pix
        ALL.HMR_pix_GAL__g[iGal] = HMR_pix
        ALL.parsecPerPixel__g[iGal] = parsecPerPixel
        ALL.at_flux_GAL__g[iGal] = at_flux_GAL
        ALL._at_flux__g.append(at_flux__z)
        ALL._at_mass__g.append(at_mass__z)
        ALL.at_flux__rg[:, iGal] = at_flux__r
        ALL.at_mass__rg[:, iGal] = at_mass__r
        ALL.at_flux_oneHLR__g[iGal] = at_flux_oneHLR
        ALL.at_mass_oneHLR__g[iGal] = at_mass_oneHLR
        
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

            tau_V__z[mask_syn__Tz[iT]] = np.ma.masked
            SFR__z[mask_syn__Tz[iT]] = np.ma.masked
            SFRSD__z[mask_syn__Tz[iT]] = np.ma.masked
            Mcor__z[mask_syn__Tz[iT]] = np.ma.masked
            McorSD__z[mask_syn__Tz[iT]] = np.ma.masked
            at_flux__z[mask_syn__Tz[iT]] = np.ma.masked
            at_mass__z[mask_syn__Tz[iT]] = np.ma.masked
            integrated_SFR = SFR__z.sum()

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
            ALL.integrated_x_Y__Tg[iT, iGal] = integrated_x_Y
            ALL.integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL.integrated_SFRSD__Tg[iT, iGal] = integrated_SFR / K.zoneArea_pc2.sum()
                        
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
            alogZ_mass_oneHLR = aux[9]
            alogZ_flux_oneHLR = aux[10]
            
            # Saving for later :D
            ALL._alogZ_mass__Ug[iU].append(alogZ_mass__z.data)
            ALL._alogZ_mass_mask__Ug[iU].append(alogZ_mass__z.mask)
            ALL._alogZ_flux__Ug[iU].append(alogZ_flux__z.data)
            ALL._alogZ_flux_mask__Ug[iU].append(alogZ_flux__z.mask)
            ALL.alogZ_mass_GAL__Ug[iU, iGal] = alogZ_mass_GAL
            ALL.alogZ_flux_GAL__Ug[iU, iGal] = alogZ_flux_GAL
            ALL.alogZ_mass__Urg[iU, :, iGal] = alogZ_mass__r
            ALL.alogZ_flux__Urg[iU, :, iGal] = alogZ_flux__r
            ALL.alogZ_mass_oneHLR__Ug[iU, iGal] = alogZ_mass_oneHLR
            ALL.alogZ_flux_oneHLR__Ug[iU, iGal] = alogZ_flux_oneHLR
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
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)

        ########## tau_V #########
        mask_tau_V_neb_aux__z = np.zeros((K.N_zone), dtype = np.bool_)
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_lines_dict__Lz[Ha_central_wl])
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_lines_dict__Lz[Hb_central_wl])
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_tau_V_neb__z)
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_tau_V_neb_err__z)
        
        print '# N_mask_tau_V_neb_aux (Ha+Hb+tauVNeb+tauVNebErr): ', mask_tau_V_neb_aux__z.astype(int).sum()

        tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = mask_tau_V_neb_aux__z)
        tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = mask_tau_V_neb_aux__z)

        # Saving for later :D
        ALL._tau_V_neb__g.append(tau_V_neb__z.data)
        ALL._tau_V_neb_err__g.append(tau_V_neb_err__z.data)
        ALL._tau_V_neb_mask__g.append(tau_V_neb__z.mask)
        ALL.integrated_tau_V_neb__g[iGal] = K.EL.integrated_tau_V_neb
        ALL.integrated_tau_V_neb_err__g[iGal] = K.EL.integrated_tau_V_neb_err
        ##########################

        ########### EW ###########
        EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        EW_Hb__z = np.ma.masked_array(K.EL.EW[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        baseline_Hb__z = K.EL.baseline[i_Hb]
        baseline_Ha__z = K.EL.baseline[i_Ha]
        integrated_baseline_Ha = np.ma.masked_array(baseline_Ha__z, mask = EW_Ha__z.mask).sum()
        integrated_baseline_Hb = np.ma.masked_array(baseline_Hb__z, mask = EW_Hb__z.mask).sum()

        # Saving for later :D        
        ALL._EW_Ha__g.append(EW_Ha__z.data)
        ALL._EW_Ha_mask__g.append(EW_Ha__z.mask)
        ALL._EW_Hb__g.append(EW_Hb__z.data)
        ALL._EW_Hb_mask__g.append(EW_Hb__z.mask)
        ALL._baseline_Hb__g.append(K.EL.baseline[i_Hb])
        ALL._baseline_Ha__g.append(K.EL.baseline[i_Ha])
        ALL.integrated_baseline_Ha__g[iGal] = integrated_baseline_Ha
        ALL.integrated_baseline_Hb__g[iGal] = integrated_baseline_Hb
        ##########################

        #### intrinsic Ha Lum ####
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun        
        L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z
        # L_int_Ha__Lz intrinsic Ha luminosity 
        eHa = np.ma.exp(q[2] * tau_V_neb__z)
        # For the zones where I don't have values for tau_V_neb I don't correct the Lum_Ha
        L_int_Ha__z = np.where(~mask_tau_V_neb_aux__z, L_obs_Ha__z * eHa, L_obs_Ha__z)
        L_int_Ha__z = np.ma.masked_array(L_int_Ha__z, mask = mask_tau_V_neb_aux__z)
        integrated_eHa = np.ma.exp(q[2] * K.EL.integrated_tau_V_neb)
        integrated_L_obs_Ha = K.EL._F_to_L(F_obs_Ha__z.sum()) / L_sun
        integrated_L_int_Ha = integrated_L_obs_Ha * integrated_eHa
        # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
        qq = q[2] / (q[0] - q[2])
        a = L_obs_Ha_err__z
        b = qq * L_obs_HaHb__z * L_obs_Hb_err__z
        L_int_Ha_err__z = np.where(~mask_tau_V_neb_aux__z, L_obs_Ha_err__z, eHa * np.sqrt(a ** 2.0 + b ** 2.0))
        L_int_Ha_err__z = np.ma.masked_array(L_int_Ha_err__z, mask = mask_tau_V_neb_aux__z)
        
        # Saving for later :D
        ALL._F_obs_Ha__g.append(F_obs_Ha__z.data)
        ALL._F_obs_Ha_mask__g.append(F_obs_Ha__z.mask)
        ALL._L_obs_Ha__g.append(L_obs_Ha__z.data)
        ALL._L_obs_Ha_mask__g.append(L_obs_Ha__z.mask)
        ALL._L_obs_Ha_err__g.append(L_obs_Ha_err__z.data)
        ALL._L_int_Ha__g.append(L_int_Ha__z.data)
        ALL._L_int_Ha_mask__g.append(L_int_Ha__z.mask)
        ALL._L_int_Ha_err__g.append(L_int_Ha_err__z.data)
        ALL.integrated_L_obs_Ha__g[iGal] = integrated_L_obs_Ha
        ALL.integrated_L_int_Ha__g[iGal] = integrated_L_int_Ha
        ##########################
        
        ###### OTH BPT LINES #####
        F_obs_Hb__z = np.ma.masked_array(K.EL.Hb_obs__z, mask = mask_lines_dict__Lz[Hb_central_wl])
        F_obs_O3__z = np.ma.masked_array(K.EL.O3_obs__z, mask = mask_lines_dict__Lz[O3_central_wl])
        F_obs_N2__z = np.ma.masked_array(K.EL.N2_obs__z, mask = mask_lines_dict__Lz[N2_central_wl])
        eHb = np.ma.exp(q[0] * tau_V_neb__z)
        eO3 = np.ma.exp(q[1] * tau_V_neb__z)
        eN2 = np.ma.exp(q[3] * tau_V_neb__z)
        F_int_Hb__z = np.where(~mask_tau_V_neb_aux__z, F_obs_Hb__z * eHb, F_obs_Hb__z)
        F_int_O3__z = np.where(~mask_tau_V_neb_aux__z, F_obs_O3__z * eO3, F_obs_O3__z)
        F_int_Ha__z = np.where(~mask_tau_V_neb_aux__z, F_obs_Ha__z * eHa, F_obs_Ha__z)
        F_int_N2__z = np.where(~mask_tau_V_neb_aux__z, F_obs_N2__z * eN2, F_obs_N2__z)
        baseline_O3__z = K.EL.baseline[i_O3]
        baseline_N2__z = K.EL.baseline[i_N2]
        maskNone = np.zeros((K.N_zone), dtype = np.bool_)
        eF_obs_Hb__z = np.ma.masked_array(K.EL.eflux[i_Hb], mask = maskNone)
        eF_obs_O3__z = np.ma.masked_array(K.EL.eflux[i_O3], mask = maskNone)
        eF_obs_Ha__z = np.ma.masked_array(K.EL.eflux[i_Ha], mask = maskNone)
        eF_obs_N2__z = np.ma.masked_array(K.EL.eflux[i_N2], mask = maskNone)
        integrated_F_obs_Hb = F_obs_Hb__z.sum()
        integrated_F_obs_O3 = F_obs_O3__z.sum()
        integrated_F_obs_Ha = F_obs_Ha__z.sum()
        integrated_F_obs_N2 = F_obs_N2__z.sum()
        integrated_F_int_Hb = integrated_F_obs_Hb * np.ma.exp(q[0] * K.EL.integrated_tau_V_neb)
        integrated_F_int_O3 = integrated_F_obs_O3 * np.ma.exp(q[1] * K.EL.integrated_tau_V_neb)
        integrated_F_int_Ha = integrated_F_obs_Ha * integrated_eHa
        integrated_F_int_N2 = integrated_F_obs_N2 * np.ma.exp(q[3] * K.EL.integrated_tau_V_neb)
        integrated_baseline_O3 = baseline_O3__z.sum()
        integrated_baseline_N2 = baseline_N2__z.sum()
        
        #integrated_tau_V
        qq = (1./(q[0] - q[2]))
        integrated_eflux_Ha = np.ma.masked_array(K.EL.eflux[i_Ha], mask = mask_tau_V_neb_aux__z).sum()
        integrated_eflux_Hb = np.ma.masked_array(K.EL.eflux[i_Hb], mask = mask_tau_V_neb_aux__z).sum()
        a = integrated_eflux_Ha / integrated_F_obs_Ha
        b = integrated_eflux_Hb / integrated_F_obs_Hb
        
        # Saving for later :D
        ALL._F_obs_Hb__g.append(F_obs_Hb__z.data)
        ALL._F_obs_O3__g.append(F_obs_O3__z.data)
        ALL._F_obs_N2__g.append(F_obs_N2__z.data)
        ALL._F_obs_Hb_mask__g.append(F_obs_Hb__z.mask)
        ALL._F_obs_O3_mask__g.append(F_obs_O3__z.mask)
        ALL._F_obs_N2_mask__g.append(F_obs_N2__z.mask)
        ALL._eF_obs_Hb__g.append(eF_obs_Hb__z.data)
        ALL._eF_obs_O3__g.append(eF_obs_O3__z.data)
        ALL._eF_obs_Ha__g.append(eF_obs_Ha__z.data)
        ALL._eF_obs_N2__g.append(eF_obs_N2__z.data)
        ALL._F_int_Ha__g.append(F_int_Ha__z)
        ALL._F_int_Hb__g.append(F_int_Hb__z)
        ALL._F_int_O3__g.append(F_int_O3__z)
        ALL._F_int_N2__g.append(F_int_N2__z)
        ALL._baseline_O3__g.append(K.EL.baseline[i_O3])
        ALL._baseline_N2__g.append(K.EL.baseline[i_N2])
        ALL.integrated_F_obs_Ha__g[iGal] = integrated_F_obs_Ha
        ALL.integrated_F_obs_Hb__g[iGal] = integrated_F_obs_Hb
        ALL.integrated_F_obs_O3__g[iGal] = integrated_F_obs_O3
        ALL.integrated_F_obs_N2__g[iGal] = integrated_F_obs_N2
        ALL.integrated_F_int_Ha__g[iGal] = integrated_F_int_Ha
        ALL.integrated_F_int_Hb__g[iGal] = integrated_F_int_Hb
        ALL.integrated_F_int_O3__g[iGal] = integrated_F_int_O3
        ALL.integrated_F_int_N2__g[iGal] = integrated_F_int_N2
        ALL.integrated_tau_V_neb_resolved__g[iGal] = qq * (np.ma.log(integrated_F_obs_Ha/integrated_F_obs_Hb) - np.log(2.86))
        ALL.integrated_tau_V_neb_err_resolved__g[iGal] = qq * (a**2. + b**2.)**0.5
        ALL.integrated_EW_Ha__g[iGal] = F_obs_Ha__z.sum()/baseline_Ha__z.sum()
        ALL.integrated_EW_Hb__g[iGal] = F_obs_Hb__z.sum()/baseline_Hb__z.sum()
        ##########################

        #### SFR and SigmaSFR ####
        # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter
        SFR_Ha__z = np.ma.masked_array(3.13 * L_int_Ha__z.data / (1.e8), mask = L_int_Ha__z.mask)
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        integrated_SFR_Ha = 3.13 * integrated_L_int_Ha / (1.e8)
        integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()

        # Saving for later :D
        ALL.integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha
        ALL.integrated_SFRSD_Ha__g[iGal] = integrated_SFRSD_Ha
        ALL._SFR_Ha__g.append(SFR_Ha__z.data)
        ALL._SFR_Ha_mask__g.append(SFR_Ha__z.mask)
        ALL._SFRSD_Ha__g.append(SFRSD_Ha__z.data)
        ALL._SFRSD_Ha_mask__g.append(SFRSD_Ha__z.mask)
        ####################################################
        ####################################################
        ####################################################
        
        # M13 Zneb calib.
        mask_Zneb_aux__z = np.zeros((K.N_zone), dtype = np.bool_)
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_tau_V_neb_aux__z)
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_lines_dict__Lz[O3_central_wl])
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_lines_dict__Lz[N2_central_wl])
        print '# N_mask_Zneb_aux__z (Ha+Hb+tauVNeb+tauVNebErr+O3+N2): ', mask_Zneb_aux__z.astype(int).sum()
        logO3N2__z = K.EL.Zneb_M13__z.copy()
        logO3N2__z[mask_Zneb_aux__z] = np.ma.masked
        ALL._logO3N2_M13__g.append(logO3N2__z.data)
        ALL._logO3N2_M13_mask__g.append(logO3N2__z.mask)
        O3Hb = integrated_F_int_O3/integrated_F_int_Hb
        N2Ha = integrated_F_int_N2/integrated_F_int_Ha
        ALL.integrated_logO3N2_M13__g[iGal] = 8.533 - 0.214 * (np.log10(O3Hb) - np.log10(N2Ha))
        
        minzones = 5

        # SYN Radial Profiles:
        for iT, tSF in enumerate(tSF__T):
            #iG = gals_not_masked - 1
            iG = -1
            ########### Now I have to define the masks
            mask_radial = mask__Tz[iT]
            print '# N_mask_radialProfiles (%.3f Myrs): %d' % ((tSF / 1e6), mask_radial.astype(int).sum())
                        
            ##########################################
            if (mask_radial.astype(int).sum() < (K.N_zone - minzones)):
                C.debug_var(args.debug, 
                            gal = califaID, 
                            iT = iT, 
                            tSF = '%.3f Myr' % (tSF / 1e6), 
                            radial_profiles = 'computing...',
                )
                
                x_Y__z = np.ma.masked_array(ALL._x_Y__Tg[iT][iG], mask = mask_radial)
                Mcor__z = np.ma.masked_array(ALL._Mcor__Tg[iT][iG], mask = mask_radial)
                McorSD__z = np.ma.masked_array(ALL._McorSD__Tg[iT][iG], mask = mask_radial)
                SFR__z = np.ma.masked_array(ALL._SFR__Tg[iT][iG], mask = mask_radial)
                SFRSD__z = np.ma.masked_array(ALL._SFRSD__Tg[iT][iG], mask = mask_radial)
                tau_V__z = np.ma.masked_array(ALL._tau_V__Tg[iT][iG], mask = mask_radial)
                at_flux__z = np.ma.masked_array(ALL._at_flux__Tg[iT][iG], mask = mask_radial)
                at_mass__z = np.ma.masked_array(ALL._at_mass__Tg[iT][iG], mask = mask_radial)
                
                #Mcor__yx = K.zoneToYX(Mcor__z, extensive = False, surface_density = False)
                #v__r, npts = K.radialProfile(Mcor__yx, Rbin__r, return_npts = True, rad_scale = K.HLR_pix, mode = 'sum', mask = mask__yx)
                #McorSD__r = v__r / (npts * K.parsecPerPixel**2.)
                #SFR__yx = K.zoneToYX(SFR__z, extensive = False, surface_density = False)
                #v__r, npts = K.radialProfile(SFR__yx, Rbin__r, return_npts = True, rad_scale = K.HLR_pix, mode = 'sum', mask = mask__yx)
                #aSFRSD__r = v__r / (npts * K.parsecPerPixel**2.)
                x_Y__r = K.zoneToRad(x_Y__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                McorSD__r = K.zoneToRad(McorSD__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                aSFRSD__r = K.zoneToRad(SFRSD__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                tau_V__r = K.zoneToRad(tau_V__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_flux__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_mass__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_flux_dezon__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale = K.HLR_pix)
                at_mass_dezon__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale = K.HLR_pix)
                
                x_Y_oneHLR = K.zoneToRad(x_Y__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                McorSD_oneHLR = K.zoneToRad(McorSD__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
                aSFRSD_oneHLR = K.zoneToRad(SFRSD__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                tau_V_oneHLR = K.zoneToRad(tau_V__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_flux_oneHLR = K.zoneToRad(at_flux__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_mass_oneHLR = K.zoneToRad(at_mass__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_flux_dezon_oneHLR = K.zoneToRad(at_flux__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
                at_mass_dezon_oneHLR = K.zoneToRad(at_mass__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
    
                Lobn__yx = K.zoneToYX(np.ma.masked_array(K.Lobn__z, mask = mask_radial), extensive = False)
                at_flux__yx = K.zoneToYX(at_flux__z, extensive = False)
                w__r = K.radialProfile(Lobn__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                v_w__r = K.radialProfile(at_flux__yx * Lobn__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                at_flux_wei__r = v_w__r / w__r
                
                Mcor__yx = K.zoneToYX(Mcor__z, extensive = False)
                at_mass__yx = K.zoneToYX(at_mass__z, extensive = False)
                w__r = K.radialProfile(Mcor__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                v_w__r = K.radialProfile(at_mass__yx * Mcor__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                at_mass_wei__r = v_w__r / w__r
        
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
                ALL.x_Y_oneHLR__Tg[iT, iGal] = x_Y_oneHLR
                ALL.McorSD_oneHLR__Tg[iT, iGal] = McorSD_oneHLR
                ALL.aSFRSD_oneHLR__Tg[iT, iGal] = aSFRSD_oneHLR
                ALL.tau_V_oneHLR__Tg[iT, iGal] = tau_V_oneHLR
                ALL.at_flux_oneHLR__Tg[iT, iGal] = at_flux_oneHLR
                ALL.at_mass_oneHLR__Tg[iT, iGal] = at_mass_oneHLR
                ALL.at_flux_dezon_oneHLR__Tg[iT, iGal] = at_flux_dezon_oneHLR
                ALL.at_mass_dezon_oneHLR__Tg[iT, iGal] = at_mass_dezon_oneHLR
                
                #NEB
                SFRSD_Ha__z = np.ma.masked_array(ALL._SFRSD_Ha__g[iG], mask = mask_radial)
                tau_V_neb__z = np.ma.masked_array(ALL._tau_V_neb__g[iG], mask = mask_radial)
                logO3N2_M13__z = np.ma.masked_array(ALL._logO3N2_M13__g[iG], mask = mask_radial)
                EW_Ha__z = np.ma.masked_array(ALL._EW_Ha__g[iG], mask = mask_radial)
                EW_Hb__z = np.ma.masked_array(ALL._EW_Hb__g[iG], mask = mask_radial)
                F_obs_Hb__z = np.ma.masked_array(ALL._F_obs_Hb__g[iG], mask = mask_radial)
                F_obs_O3__z = np.ma.masked_array(ALL._F_obs_O3__g[iG], mask = mask_radial)
                F_obs_Ha__z = np.ma.masked_array(ALL._F_obs_Ha__g[iG], mask = mask_radial)
                F_obs_N2__z = np.ma.masked_array(ALL._F_obs_N2__g[iG], mask = mask_radial)
                F_int_Hb__z = np.ma.masked_array(ALL._F_int_Hb__g[iG], mask = mask_radial)
                F_int_O3__z = np.ma.masked_array(ALL._F_int_O3__g[iG], mask = mask_radial)
                F_int_Ha__z = np.ma.masked_array(ALL._F_int_Ha__g[iG], mask = mask_radial)
                F_int_N2__z = np.ma.masked_array(ALL._F_int_N2__g[iG], mask = mask_radial)
                baseline_Ha__z = np.ma.masked_array(ALL._baseline_Ha__g[iG], mask = mask_radial)
                baseline_Hb__z = np.ma.masked_array(ALL._baseline_Hb__g[iG], mask = mask_radial)
                
                tau_V_neb__r = K.zoneToRad(tau_V_neb__z, Rbin__r,  rad_scale = K.HLR_pix, extensive = False)
                tau_V_neb_oneHLR = K.zoneToRad(tau_V_neb__z, Rbin_oneHLR,  rad_scale = K.HLR_pix, extensive = False)
                logO3N2_M13__r = K.zoneToRad(logO3N2_M13__z, Rbin__r,  rad_scale = K.HLR_pix, extensive = False)
                logO3N2_M13_oneHLR = K.zoneToRad(logO3N2_M13__z, Rbin_oneHLR,  rad_scale = K.HLR_pix, extensive = False)
                aSFRSD_Ha__r = K.zoneToRad(SFRSD_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                aSFRSD_Ha_oneHLR = K.zoneToRad(SFRSD_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)

                EW_Ha__r = K.zoneToRad(EW_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                EW_Hb__r = K.zoneToRad(EW_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                EW_Ha_oneHLR = K.zoneToRad(EW_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                EW_Hb_oneHLR = K.zoneToRad(EW_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
    
                F_Ha__yx = K.zoneToYX(F_obs_Ha__z, extensive = False)
                baseline_Ha__yx = K.zoneToYX(baseline_Ha__z, extensive = False)
                v__r = K.radialProfile(F_Ha__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                w__r = K.radialProfile(baseline_Ha__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                EW_Ha_wei__r = v__r / w__r
    
                F_Hb__yx = K.zoneToYX(F_obs_Hb__z, extensive = False)
                baseline_Hb__yx = K.zoneToYX(baseline_Hb__z, extensive = False)
                v__r = K.radialProfile(F_Hb__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                w__r = K.radialProfile(baseline_Hb__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                EW_Hb_wei__r = v__r / w__r
    
                F_obs_Ha__r = K.zoneToRad(F_obs_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Hb__r = K.zoneToRad(F_obs_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_O3__r = K.zoneToRad(F_obs_O3__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_N2__r = K.zoneToRad(F_obs_N2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_Ha__r = K.zoneToRad(F_int_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_Hb__r = K.zoneToRad(F_int_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_O3__r = K.zoneToRad(F_int_O3__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_N2__r = K.zoneToRad(F_int_N2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Ha_oneHLR = K.zoneToRad(F_obs_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Hb_oneHLR = K.zoneToRad(F_obs_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_O3_oneHLR = K.zoneToRad(F_obs_O3__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_N2_oneHLR = K.zoneToRad(F_obs_N2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_Ha_oneHLR = K.zoneToRad(F_int_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_Hb_oneHLR = K.zoneToRad(F_int_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_O3_oneHLR = K.zoneToRad(F_int_O3__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_N2_oneHLR = K.zoneToRad(F_int_N2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                
                ALL.tau_V_neb__Trg[iT, :, iGal] = tau_V_neb__r
                ALL.logO3N2_M13__Trg[iT, :, iGal] = logO3N2_M13__r
                ALL.EW_Ha__Trg[iT, :, iGal] = EW_Ha__r
                ALL.EW_Hb__Trg[iT, :, iGal] = EW_Hb__r
                ALL.EW_Ha_wei__Trg[iT, :, iGal] = EW_Ha_wei__r
                ALL.EW_Hb_wei__Trg[iT, :, iGal] = EW_Hb_wei__r
                ALL.F_obs_Ha__Trg[iT, :, iGal] = F_obs_Ha__r
                ALL.F_obs_Hb__Trg[iT, :, iGal] = F_obs_Hb__r
                ALL.F_obs_O3__Trg[iT, :, iGal] = F_obs_O3__r
                ALL.F_obs_N2__Trg[iT, :, iGal] = F_obs_N2__r
                ALL.F_int_Ha__Trg[iT, :, iGal] = F_int_Ha__r
                ALL.F_int_Hb__Trg[iT, :, iGal] = F_int_Hb__r
                ALL.F_int_O3__Trg[iT, :, iGal] = F_int_O3__r
                ALL.F_int_N2__Trg[iT, :, iGal] = F_int_N2__r
                ALL.aSFRSD_Ha__Trg[iT, :, iGal] = aSFRSD_Ha__r
                ALL.tau_V_neb_oneHLR__Tg[iT, iGal] = tau_V_neb_oneHLR
                ALL.logO3N2_M13_oneHLR__Tg[iT, iGal] = logO3N2_M13_oneHLR
                ALL.EW_Ha_oneHLR__Tg[iT, iGal] = EW_Ha_oneHLR
                ALL.EW_Hb_oneHLR__Tg[iT, iGal] = EW_Hb_oneHLR
                ALL.F_obs_Ha_oneHLR__Tg[iT, iGal] = F_obs_Ha_oneHLR
                ALL.F_obs_Hb_oneHLR__Tg[iT, iGal] = F_obs_Hb_oneHLR
                ALL.F_obs_O3_oneHLR__Tg[iT, iGal] = F_obs_O3_oneHLR
                ALL.F_obs_N2_oneHLR__Tg[iT, iGal] = F_obs_N2_oneHLR
                ALL.F_int_Ha_oneHLR__Tg[iT, iGal] = F_int_Ha_oneHLR
                ALL.F_int_Hb_oneHLR__Tg[iT, iGal] = F_int_Hb_oneHLR
                ALL.F_int_O3_oneHLR__Tg[iT, iGal] = F_int_O3_oneHLR
                ALL.F_int_N2_oneHLR__Tg[iT, iGal] = F_int_N2_oneHLR
                ALL.aSFRSD_Ha_oneHLR__Tg[iT, iGal] = aSFRSD_Ha_oneHLR
        if args.gasprop:
            K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)

    t_init_stack = time.clock()        
    ALL.stack_zones_data()
    ALL.integrated_mask()
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

    
