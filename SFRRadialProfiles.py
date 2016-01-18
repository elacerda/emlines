#!/usr/bin/python
#
# Lacerda@Saco - 9/Jan/2016
#
import sys
import time
import numpy as np
import argparse as ap
import CALIFAUtils as C

class tupperware(object): 
    def __init__(self):
        pass
    
    def __getattr__(self, attr):
        r = self.__dict__.get(attr, None)
        return r
    
    def mask_gal(self, iGal):
        for v in self.__dict__.keys():
            if isinstance(self.__dict__[v], np.ma.core.MaskedArray):
                self.__dict__[v][..., iGal] = np.ma.masked

    def new_empty(self, attr, **kw):
        shape = kw.get('shape', 1)
        dtype = kw.get('dtype', float)
        setattr(self, attr, np.ma.masked_all(shape, dtype = dtype)) 

def mask_zones_iT(iT, H, bamin, gals_slice):
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # mask__g = np.bitwise_or(H.SFRSD_Ha__g.mask, H.SFRSD__Tg[iT].mask)
    # mask__g = np.bitwise_or(mask__g, H.tau_V__Tg[iT].mask)
    # mask__g = np.bitwise_or(mask__g, H.tau_V_neb__g.mask)
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD_Ha__g * 1e6).mask, np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.my_O_O3N2_M13__g.mask)
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), bamin))
    mask__g = np.bitwise_or(mask__g, ~gals_slice)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    return mask__g

def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return False
    if EL == True and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return False
    if GP is True and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return False
    # Problem in FITS file
    return True    

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default_args = {
        'debug' : False,
        'hdf5' : None,
        'slice_gals' : None,
        'bamin' : 0,
        'rbinini' : 0.0,
        'rbinfin' : 2.0,
        'rbinstep' : 0.1,
        'weiradprof' : False,
        'output' : 'SFRRadialProfiles',
        'v_run' : 'v20_q050.d15a',
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
    parser.add_argument('--v_run',
                        metavar = 'RUNCODE',
                        type = str,
                        default = default_args['v_run'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['slice_gals'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default_args['output'])
    parser.add_argument('--weiradprof', '-W',
                        action = 'store_true',
                        default = default_args['weiradprof'])
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
    parser.add_argument('--bamin',
                        metavar = 'FLOAT',
                        type = float,
                        default = default_args['bamin'])


    return parser.parse_args()

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE


Zsun = 0.019

if __name__ == '__main__':
    t_init_gal = time.clock()

    args = parser_args()
    debug = args.debug
    h5fname = args.hdf5
    
    # Creating radial bins.
    Rbin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    RbinCenter__r = (Rbin__r[:-1] + Rbin__r[1:]) / 2.0
    NRbins = len(RbinCenter__r)
    RColor = [ 'r', 'y', 'b', 'k' ]
    RRange = [  .5, 1., 1.5, 2.  ]
    Rbin_oneHLR = 1. + np.asarray([ -1.5, -0.5, 0.5, 1.5 ]) * args.rbinstep
    
    H = C.H5SFRData(args.hdf5)
    tSF__T = H.tSF__T[10:12]
    tZ__U = H.tZ__U
    N_g = H.N_gals
    N_T = 2
    N_U = H.N_U
    shape2d = (N_T, N_g)
    shape3dT = (N_T, NRbins, N_g)
    shape3dU = (N_U
                , NRbins, N_g)
    Dvars = dict(
        #aZ_mass__Urg = dict(shape = shape3dU, dtype = float),
        #L_int_Ha__Trg = dict(shape = shape3dT, dtype = float),
        #at_flux_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        #aZ_mass_oneHLR__Ug = dict(shape = shape2d, dtype = float),
        #L_int_Ha_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        aSFRSD__Trg = dict(shape = shape3dT, dtype = float),
        McorSD__Trg = dict(shape = shape3dT, dtype = float),
        tau_V__Trg = dict(shape = shape3dT, dtype = float),
        EW_Ha__Trg = dict(shape = shape3dT, dtype = float),
        EW_Hb__Trg = dict(shape = shape3dT, dtype = float),
        aSFRSD_Ha__Trg = dict(shape = shape3dT, dtype = float),
        tau_V_neb__Trg = dict(shape = shape3dT, dtype = float),
        O_O3N2_M13__Trg = dict(shape = shape3dT, dtype = float),
        aSFRSD_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        McorSD_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        tau_V_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        EW_Ha_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        EW_Hb_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        aSFRSD_Ha_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        tau_V_neb_oneHLR__Tg = dict(shape = shape2d, dtype = float),
        O_O3N2_M13_oneHLR__Tg = dict(shape = shape2d, dtype = float),
    )

    # Creating variables    
    radial_profiles = tupperware()
    for k, v in Dvars.iteritems():
        radial_profiles.new_empty(k, **dict(shape = v['shape'], dtype = v['dtype']))
    
    fnamesuffix = '.pdf'
    
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1]
        fnamesuffix = '_%s%s' % (gals_txt, fnamesuffix)
        #integrated
        l_gals, _ = C.sort_gals(args.slice_gals)
        l_gals = sorted(l_gals)
        gals_slice__integr = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
        for g in l_gals:
            i = H.califaIDs_all.tolist().index(g)
            gals_slice__integr[i] = True

    mask__Tg = np.zeros((N_T, H.N_zones__g.sum()), dtype = bool)
    for iT in xrange(N_T):
        mask__Tg[iT] = mask_zones_iT(iT, H, args.bamin, gals_slice__g)
            
    for iGal, califaID in enumerate(H.califaIDs):
        K = C.read_one_cube(califaID, EL = True, GP = True, debug = args.debug, v_run = args.v_run)
        tipos, tipo, tipo_m, tipo_p = C.get_morfologia(califaID)
        my_type = C.my_morf(tipos)
        N_zone = K.N_zone
        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        minzones = 5
        
        for iT, tSF in enumerate(tSF__T):
            mask__z, where_slice = H.get_prop_gal(mask__Tg[iT], califaID, return_slice = True)
            N_mask = mask__z.astype(int).sum() 
            C.debug_var(args.debug, N_mask = N_mask, N_zone = K.N_zone)
            if (N_mask < (K.N_zone - minzones)):
                mask__yx = K.zoneToYX(mask__z, extensive = False, surface_density = False)
                SFR__z = H.SFR__Tg[iT][where_slice]
                zoneArea_pc2__z = H.zone_area_pc2__g[where_slice]
                SFR_Ha__z = H.SFR_Ha__g[where_slice]
                Mcor__z = H.Mcor__Tg[iT][where_slice]
                at_flux__z = H.at_flux__Tg[iT][where_slice]
                tau_V__z = H.tau_V__Tg[iT][where_slice]
                tau_V_neb__z = H.tau_V_neb__g[where_slice]
                EW_Ha__z = H.EW_Ha__g[where_slice]
                EW_Hb__z = H.EW_Hb__g[where_slice]
                F_int_Hb__z = H.F_int_Hb__g[where_slice]
                F_int_O3__z = H.F_int_O3__g[where_slice]
                F_int_Ha__z = H.F_int_Ha__g[where_slice]
                F_int_N2__z = H.F_int_N2__g[where_slice]
                O3__z = F_int_O3__z / F_int_Hb__z
                N2__z = F_int_N2__z / F_int_Ha__z
                logO3N2__z = np.ma.log10(O3__z) - np.ma.log10(N2__z)
                O_O3N2_M13__z = 8.533 - 0.214 * logO3N2__z
                #print califaID, tau_V__z.max(), tau_V__z.min(), tau_V_neb__z.max(), tau_V_neb__z.min()
                radial_profiles.aSFRSD__Trg[iT, :, iGal] = K.zoneToRad(SFR__z/zoneArea_pc2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact') 
                radial_profiles.McorSD__Trg[iT, :, iGal] = K.zoneToRad(Mcor__z/zoneArea_pc2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.tau_V__Trg[iT, :, iGal] = K.zoneToRad(tau_V__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.EW_Ha__Trg[iT, :, iGal] = K.zoneToRad(EW_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.EW_Hb__Trg[iT, :, iGal] = K.zoneToRad(EW_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.aSFRSD_Ha__Trg[iT, :, iGal] = K.zoneToRad(SFR_Ha__z/zoneArea_pc2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.tau_V_neb__Trg[iT, :, iGal] = K.zoneToRad(tau_V_neb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.O_O3N2_M13__Trg[iT, :, iGal] = K.zoneToRad(O_O3N2_M13__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)#, mode = 'mean_exact')
                radial_profiles.aSFRSD_oneHLR__Tg[iT, iGal] = K.zoneToRad(SFR__z/zoneArea_pc2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.McorSD_oneHLR__Tg[iT, iGal] = K.zoneToRad(Mcor__z/zoneArea_pc2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.tau_V_oneHLR__Tg[iT, iGal] = K.zoneToRad(tau_V__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.EW_Ha_oneHLR__Tg[iT, iGal] = K.zoneToRad(EW_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.EW_Hb_oneHLR__Tg[iT, iGal] = K.zoneToRad(EW_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.aSFRSD_Ha_oneHLR__Tg[iT, iGal] = K.zoneToRad(SFR_Ha__z/zoneArea_pc2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.tau_V_neb_oneHLR__Tg[iT, iGal] = K.zoneToRad(tau_V_neb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
                radial_profiles.O_O3N2_M13_oneHLR__Tg[iT, iGal] = K.zoneToRad(O_O3N2_M13__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False, surface_density = False, mask = mask__yx)[1]#, mode = 'mean_exact')[1]
        K.GP.close()
        K.EL.close()
        K.close()
        del K
    
    if args.output != None:
        import h5py
        filename = args.output
        h5 = h5py.File(filename, 'w')
        D = {}
        for k in radial_profiles.__dict__.keys():
            #print k, k[-3:], k[-2:]
            if (k[-3:] == 'Trg') | (k[-2:] == 'Tg'):
                #print 'in'
                tmp_data = {'masked/data/%s' % k : (getattr(radial_profiles, k)).data }
                tmp_mask = {'masked/mask/%s' % k : (getattr(radial_profiles, k)).mask }
                D.update(tmp_data)
                D.update(tmp_mask)
        D['data/RbinIni'] = args.rbinini
        D['data/RbinFin'] = args.rbinfin
        D['data/RbinStep'] = args.rbinstep
        D['data/Rbin__r'] = Rbin__r
        D['data/RbinCenter__r'] = RbinCenter__r
        D['data/NRbins'] = NRbins
        D['data/RRange'] = RRange
        D['data/tSF__T'] = tSF__T
        D['data/tZ__U'] = tZ__U
        D['data/xOkMin'] = H.xOkMin
        D['data/tauVOkMin'] = H.tauVOkMin
        D['data/tauVNebOkMin'] = H.tauVNebOkMin
        D['data/tauVNebErrMax'] = H.tauVNebErrMax
        for k in D.keys():
            print k
            try:
                h5.create_dataset(k, data = D[k], compression = 'gzip', compression_opts = 4)
            except TypeError:
                h5.create_dataset(k, data = D[k])
        h5.close()
