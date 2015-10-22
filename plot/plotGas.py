#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
from CALIFAUtils.scripts import get_CALIFAID_by_NEDName
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plot_text_ax
#from CALIFAUtils.plots import plot_zbins
from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

RNuc = 0.5

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'scatter' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--scatter', '-s',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--itSF', '-T',
                        help = 'SF age index.',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    parser.add_argument('--maskradius', '-R',
                        help = 'initial RDisc value in HLR',
                        metavar = 'NUM',
                        type = float,
                        default = default['maskradius'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default['output'])

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    H = C.H5SFRData(args.hdf5)
    iT = args.itSF
    iU = -1

    minR = 0
    
    if args.maskradius is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
        maskRadiusOk__g = (H.zone_dist_HLR__g >= args.maskradius) & (H.zone_dist_HLR__g <= H.Rbin__r[-1]) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= args.maskradius)).T
        
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
    
    Area_GAL__g = (H.Mcor_GAL__g / H.McorSD_GAL__g)
    dustdim = 0.2 # md / rhod
    
    ######################
    # SK Law 
    ######################
    SK_zero = 1.6e-4
    SK_slope = 1.4
    aux = 1e6 ** (1. / SK_slope)
    SK_SigmaGas__g = aux * (H.SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__g = aux * (H.SFRSD_Ha__g / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas__rg = aux * (H.aSFRSD__Trg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__rg = aux * (H.aSFRSD_Ha__rg / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas = aux * (H.integrated_SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas_Ha = aux * (H.integrated_SFRSD_Ha__g / SK_zero) ** (1. / SK_slope) 
    SK_DGR__g = dustdim * H.tau_V__Tg[iT] / SK_SigmaGas__g
    SK_DGR_Ha__g = dustdim * H.tau_V_neb__g / SK_SigmaGas_Ha__g
    SK_DGR__rg = dustdim * H.tau_V__Trg[iT] / SK_SigmaGas__rg
    SK_DGR_Ha__rg = dustdim * H.tau_V_neb__rg / SK_SigmaGas_Ha__rg
    SK_integrated_DGR = dustdim * H.integrated_tau_V__g / SK_integrated_SigmaGas
    SK_integrated_DGR_Ha = dustdim * H.integrated_tau_V_neb__g / SK_integrated_SigmaGas_Ha
    SK_GSR__g = SK_SigmaGas__g / H.McorSD__Tg[iT]
    SK_GSR_Ha__g = SK_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    SK_GSR__rg = SK_SigmaGas__rg / H.McorSD__Trg[iT]
    SK_GSR_Ha__rg = SK_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    SK_integrated_GSR = SK_integrated_SigmaGas / H.McorSD_GAL__g
    SK_integrated_GSR_Ha = SK_integrated_SigmaGas_Ha / H.McorSD_GAL__g
    SK_f_gas__g = 1. / (1. + 1. / SK_GSR__g)
    SK_f_gas_Ha__g = 1. / (1. + 1. / SK_GSR_Ha__g)
    SK_f_gas__rg = 1. / (1. + 1. / SK_GSR__rg)
    SK_f_gas_Ha__rg = 1. / (1. + 1. / SK_GSR_Ha__rg)
    SK_integrated_f_gas = 1. / (1. + 1. / SK_integrated_GSR)
    SK_integrated_f_gas_Ha = 1. / (1. + 1. / SK_integrated_GSR_Ha)
    ######################
    
    ######################
    # Remy-Ruyer
    ######################
    RR_DGR = 10. ** (-2.21) # 0.006166
    RR_SigmaGas__g = dustdim * H.tau_V__Tg[iT] / RR_DGR
    RR_SigmaGas_Ha__g = dustdim * H.tau_V_neb__g / RR_DGR
    RR_SigmaGas__rg = dustdim * H.tau_V__Trg[iT] / RR_DGR
    RR_SigmaGas_Ha__rg = dustdim * H.tau_V_neb__rg / RR_DGR
    RR_integrated_SigmaGas = dustdim * H.integrated_tau_V__g / RR_DGR
    RR_integrated_SigmaGas_Ha = dustdim * H.integrated_tau_V_neb__g / RR_DGR
    RR_GSR__g = RR_SigmaGas__g / H.McorSD__Tg[iT]
    RR_GSR_Ha__g = RR_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    RR_GSR__rg = RR_SigmaGas__rg / H.McorSD__Trg[iT]
    RR_GSR_Ha__rg = RR_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    RR_integrated_GSR = RR_integrated_SigmaGas / H.McorSD_GAL__g
    RR_integrated_GSR_Ha = RR_integrated_SigmaGas_Ha / H.McorSD_GAL__g 
    RR_f_gas__g = 1. / (1. + 1. / RR_GSR__g)
    RR_f_gas_Ha__g = 1. / (1. + 1. / RR_GSR_Ha__g)
    RR_f_gas__rg = 1. / (1. + 1. / RR_GSR__rg)
    RR_f_gas_Ha__rg = 1. / (1. + 1. / RR_GSR_Ha__rg)
    RR_integrated_f_gas = 1. / (1. + 1. / RR_integrated_GSR)
    RR_integrated_f_gas_Ha = 1. / (1. + 1. / RR_integrated_GSR_Ha)
    ######################
    
    ######################
    #Brinchmann
    ######################
    DGR_conv_lim_sup = 1.1e-2
    DGR_conv_lim_inf = 5.3e-3
    DGR_interval = np.array([DGR_conv_lim_inf, DGR_conv_lim_sup])
    DGR_cte = DGR_interval.mean()
    OHSunBrinch_inv = 1 / (10.**(8.82 - 12))
    ######################
    # from OLS Bisector
    p_ols = np.array([1.87, -6.98])
    BR_OHBrinch_ols__g = np.ma.masked_all((H.O_O3N2_M13__g.shape))
    BR_OHBrinch_ols__g[~H.O_O3N2_M13__g.mask] = np.polyval(p_ols, H.O_O3N2_M13__g.compressed()) 
    BR_OHBrinch_ols__rg = np.ma.masked_all((H.O_O3N2_M13__rg.shape))
    BR_OHBrinch_ols__rg[~H.O_O3N2_M13__rg.mask] = np.polyval(p_ols, H.O_O3N2_M13__rg.compressed()) 
    BR_DGR_up_ols__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_ols__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_ols__g = DGR_cte * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_ols__rg = DGR_cte * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_up_ols__g
    BR_SigmaGas_up_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_up_ols__rg
    BR_SigmaGas_Ha_up_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_up_ols__g
    BR_SigmaGas_Ha_up_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_up_ols__rg
    BR_SigmaGas_down_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_down_ols__g
    BR_SigmaGas_down_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_down_ols__rg
    BR_SigmaGas_Ha_down_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_down_ols__g
    BR_SigmaGas_Ha_down_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_down_ols__rg
    BR_SigmaGas_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_ols__g
    BR_SigmaGas_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_ols__rg
    BR_SigmaGas_Ha_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_ols__g
    BR_SigmaGas_Ha_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_ols__rg
    BR_GSR_up_ols__g = BR_SigmaGas_up_ols__g / H.McorSD__Tg[iT]
    BR_GSR_up_ols__rg = BR_SigmaGas_up_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_up_ols__g = BR_SigmaGas_Ha_up_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_up_ols__rg = BR_SigmaGas_Ha_up_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_down_ols__g = BR_SigmaGas_down_ols__g / H.McorSD__Tg[iT]
    BR_GSR_down_ols__rg = BR_SigmaGas_down_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_down_ols__g = BR_SigmaGas_Ha_down_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_down_ols__rg = BR_SigmaGas_Ha_down_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_ols__g = BR_SigmaGas_ols__g / H.McorSD__Tg[iT]
    BR_GSR_ols__rg = BR_SigmaGas_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_ols__g = BR_SigmaGas_Ha_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_ols__rg = BR_SigmaGas_Ha_ols__rg / H.McorSD__Trg[iT]
    BR_f_gas_up_ols__g = 1. / (1. + 1. / BR_GSR_up_ols__g)
    BR_f_gas_up_ols__rg = 1. / (1. + 1. / BR_GSR_up_ols__rg)
    BR_f_gas_Ha_up_ols__g = 1. / (1. + 1. / BR_GSR_Ha_up_ols__g)
    BR_f_gas_Ha_up_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_up_ols__rg)
    BR_f_gas_down_ols__g = 1. / (1. + 1. / BR_GSR_down_ols__g)
    BR_f_gas_down_ols__rg = 1. / (1. + 1. / BR_GSR_down_ols__rg)
    BR_f_gas_Ha_down_ols__g = 1. / (1. + 1. / BR_GSR_Ha_down_ols__g)
    BR_f_gas_Ha_down_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_down_ols__rg)
    BR_f_gas_ols__g = 1. / (1. + 1. / BR_GSR_ols__g)
    BR_f_gas_ols__rg = 1. / (1. + 1. / BR_GSR_ols__rg)
    BR_f_gas_Ha_ols__g = 1. / (1. + 1. / BR_GSR_Ha_ols__g)
    BR_f_gas_Ha_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_ols__rg)
    #### integrated ####
    aux = np.polyval(p_ols, H.integrated_O_O3N2_M13__g)
    BR_integrated_OHBrinch_ols = np.ma.masked_array(aux, mask = (np.isnan(aux) | H.integrated_O_O3N2_M13__g.mask))
    BR_integrated_DGR_ols = DGR_cte * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_integrated_DGR_up_ols = DGR_conv_lim_sup * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_integrated_DGR_down_ols = DGR_conv_lim_inf * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_integrated_SigmaGas_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_ols
    BR_integrated_SigmaGas_up_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_up_ols
    BR_integrated_SigmaGas_down_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_down_ols
    BR_integrated_SigmaGas_Ha_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_ols
    BR_integrated_SigmaGas_Ha_up_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_up_ols
    BR_integrated_SigmaGas_Ha_down_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_down_ols
    BR_integrated_GSR_ols = BR_integrated_SigmaGas_ols / H.McorSD_GAL__g
    BR_integrated_GSR_up_ols = BR_integrated_SigmaGas_up_ols / H.McorSD_GAL__g
    BR_integrated_GSR_down_ols = BR_integrated_SigmaGas_down_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_ols = BR_integrated_SigmaGas_Ha_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_up_ols = BR_integrated_SigmaGas_Ha_up_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_down_ols = BR_integrated_SigmaGas_Ha_down_ols / H.McorSD_GAL__g
    BR_integrated_f_gas_ols = 1. / (1. + 1. / BR_integrated_GSR_ols)
    BR_integrated_f_gas_up_ols = 1. / (1. + 1. / BR_integrated_GSR_up_ols)
    BR_integrated_f_gas_down_ols = 1. / (1. + 1. / BR_integrated_GSR_down_ols)
    BR_integrated_f_gas_Ha_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_ols)
    BR_integrated_f_gas_Ha_up_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_up_ols)
    BR_integrated_f_gas_Ha_down_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_down_ols)
    
    ######################
    # from cubic polynomial fit
    #p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    ######################
    
    ####################################
    #### Mass from CALIFA Datafiles ####
    ####################################
    ####################################
    '''
    FILE: M_H1_CALIFA.csv - Miguel A. Perez
    delimiter = ','
    comment = '#'
    columns:
        1 - CALIFA No
        2 - NED Name
        3 - Distance
        4 - RA(J2000.0)
        5 - DEC(J2000.0)
        6 - Sobs
        7 - Scor
        8 - Sabs
        9 - e_Sabs,
        10 - M(HI)abs
        11 - e_M(HI)abs
        
    FILE: CALIFA_HI_angel.dat - Angel R. Lopez-Sanchez
    delimiter = ','
    comments = ';'
    columns:
        1 - NED Name
        2 - redshift
        3 - log(Ms)
        4 - e_log(Ms)
        5 - 12+log(O/H)
        6 - e_12+log(O/H)
        7 - log(SFR)
        8 - e_log(SFR)
        9 - log(Mg = 1.32 MHI)
        10 - e_log(Mg = 1.32 MHI)
    '''
    dirs = C.CALIFAPaths()
    file_miguel = '%sM_HI_CALIFA.csv' % dirs.califa_work_dir
    dtype_miguel = np.dtype([('califaID', '|S5'), ('M_HI', np.float)])
    read_miguel = np.loadtxt(file_miguel, delimiter = ',', usecols = [0, 9], dtype = dtype_miguel)
    map_miguel = {}
    for i, g in enumerate(read_miguel['califaID']):
        map_miguel[g] = i
    aux = set(H.califaIDs.tolist())
    gals_miguel_intersect = sorted([g for g in map_miguel.keys() if g in aux])
    gals_miguel_slice = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
    integrated_M_HI_miguel__g = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    for g in gals_miguel_intersect:
        i = H.califaIDs_all.tolist().index(g)
        i_r = map_miguel[g]
        #print g, i, i_r
        gals_miguel_slice[i] = True
        integrated_M_HI_miguel__g[i] = read_miguel['M_HI'][i_r]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # file_angel = '%sCALIFA_HI_angel.dat' % dirs.califa_work_dir
    # dtype_angel = np.dtype([('NEDName', '|S15'), ('logMg', np.float)])
    # read_angel = np.loadtxt(file_angel, delimiter = ',', comments=';', usecols = [0, 8], dtype = dtype_angel)
    # map_angel = {}
    # for i, g in enumerate(read_angel['NEDName']):
    #     k = get_CALIFAID_by_NEDName(g)
    #     if k: map_angel[k[0]] = i
    # gals_angel_intersect = sorted(np.intersect1d(map_angel.keys(), H.califaIDs).tolist())
    # gals_angel_slice = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
    # integrated_log_M_gas_angel__g = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    # for g in gals_angel_intersect:
    #     i = H.califaIDs_all.tolist().index(g)
    #     i_r = map_angel[g]
    #     gals_angel_slice[i] = True
    #     integrated_log_M_gas_angel__g[i] = read_angel['logMg'][i_r]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    integrated_f_gas_miguel = 1. / (1. + 1 / ((integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_miguel = 1./(1. + 1/((1.32 * integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_angel = 1./(1. + 1/((10. ** integrated_log_M_gas_angel__g) / H.Mcor_GAL__g))

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # logMHII_NEDName_sanitycheck = {
    #     'NGC5682' : { 'MH2' : 6.72, 'MHI' : 9.17 },
    #     'UGC07012': { 'MH2' : 7.63, 'MHI' : 9.85 },
    #     'NGC3381' : { 'MH2' : 7.27, 'MHI' : 9.42 },
    #     'NGC6168' : { 'MH2' : 8.04, 'MHI' : 9.45 },
    #     'NGC4961' : { 'MH2' : 7.91, 'MHI' : 9.77 },
    #     'NGC5205' : { 'MH2' : 7.35, 'MHI' : 9.50 },
    #     'NGC2347' : { 'MH2' : 8.96, 'MHI' : 10.13 },
    #     'NGC6060' : { 'MH2' : 9.10, 'MHI' : 10.40 },
    #     'NGC2410' : { 'MH2' : 8.91, 'MHI' : None },
    #     'NGC2639' : { 'MH2' : 8.85, 'MHI' : None },
    #     'NGC5614' : { 'MH2' : 9.43, 'MHI' : 9.49 },
    #     'NGC1167' : { 'MH2' : None, 'MHI' : 9.98 },
    # }
    # 
    # M_gas_sanitycheck = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    # f_gas_sanitycheck = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    # 
    # for k, v in logMHII_NEDName_sanitycheck.iteritems():
    #     g = get_CALIFAID_by_NEDName(k)
    #     if g:
    #         logMHII = v['MH2']
    #         logMHI = v['MHI']
    #         if logMHI != None and logMHII != None:
    #             M_gas = 10.**logMHI + 10.**logMHII
    #             i = H.califaIDs_all.tolist().index(g[0])
    #             f_gas = 1. / (1. + 1./(M_gas/H.Mcor_GAL__g[i]))
    #             #print g, k, '%.2f' % np.log10(H.Mcor_GAL__g[i])
    #             M_gas_sanitycheck[i] = M_gas
    #             f_gas_sanitycheck[i] = f_gas
    #         
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    Area_GAL__g = (H.Mcor_GAL__g / H.McorSD_GAL__g)
    #SigmaGas_sanitycheck = M_gas_sanitycheck / Area_GAL__g

    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # print '#CALIFAID,EDGE,ANGEL,MIGUEL,BR,SK'
    # for i, g in enumerate(H.califaIDs_all):
    #     if not H.califaIDs_all.mask[i]:
    #         print '%s & %.4f & %.4f & %.4f & %.4f & %.4f \\' % (g, 
    #                                                             f_gas_sanitycheck[i], 
    #                                                             integrated_f_gas_angel[i], 
    #                                                             integrated_f_gas_miguel[i],
    #                                                             BR_integrated_f_gas_Ha_ols[i],
    #                                                             SK_integrated_f_gas_Ha[i])
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    
    default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)
    
    BR_f_gas_at_oneHLR__g = BR_f_gas_ols__rg[9, :] + BR_f_gas_ols__rg[10, :]
    BR_f_gas_at_oneHLR__g /= 2. 
    BR_f_gas_Ha_at_oneHLR__g = BR_f_gas_Ha_ols__rg[9, :] + BR_f_gas_Ha_ols__rg[10, :]
    BR_f_gas_Ha_at_oneHLR__g /= 2.
    SK_f_gas_at_oneHLR__g = SK_f_gas__rg[9, :] + SK_f_gas__rg[10, :]
    SK_f_gas_at_oneHLR__g /= 2.
    SK_f_gas_Ha_at_oneHLR__g = SK_f_gas_Ha__rg[9, :] + SK_f_gas_Ha__rg[10, :]
    SK_f_gas_Ha_at_oneHLR__g /= 2.

    if args.dryrun is True:
        sys.exit('dryrun')

 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 #    if args.debug is True:
 #        import pandas as pd
 #        import seaborn as sns
 #        from scipy.stats import spearmanr
 #        #sns.set(style = 'white')
 #        sns.set_context('talk', font_scale = 1.)
 #         
 #        # Generate a random correlated bivariate dataset
 #        xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Ug[-1], y = np.ma.log10(BR_f_gas_Ha_ols__g))
 #        x1 = pd.Series(xm.compressed(), name = r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]')
 #        x2 = pd.Series(ym.compressed(), name = r'$\log\ f_{gas}$')
 #         
 #        # Show the joint distribution using kernel density estimation
 #        NRows = 1
 #        NCols = 4
 #        page_size_inches = (NCols * 3.5, NRows * 4)
 #        grid_shape = (NRows, NCols)
 # 
 #        f = plt.figure()
 #        f.set_size_inches(page_size_inches)
 #        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
 #        ax.set_axis_on()
 #        sns.kdeplot(x1, x2, cmap = 'Blues', shade = True, shade_lowest = False)
 #        sns.rugplot(x1, color = 'k', ax = ax)
 #        #sns.rugplot(x2, vertical = True, color = 'k', ax = ax_joint)
 #        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 #        # g = sns.jointplot(
 #        #     x1,
 #        #     x2,
 #        #     space = 0,
 #        #     stat_func = spearmanr,
 #        #     kind = 'kde',
 #        #     joint_kws = dict(s = 5, edgecolor = 'none'),
 #        #     annot_kws = dict(stat = r'$R_s$', fontsize = 10),
 #        #     marginal_kws = dict(color = ".5"),
 #        # )
 #        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 #        rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, smooth = True, sigma = 1.2, overlap = 0.4)
 #        ax.plot(rs.xS, rs.yS, '.-', c = 'w')
 #        #g = g.plot_joint(sns.kdeplot, n_levels=10, cmap = 'Blues_r')
 #        # f = plt.gcf()
 #        # f.set_dpi(100)
 #        # f.set_size_inches(10, 8)
 #        #g.fig.savefig('jointplot.png')
 #        f.savefig('jointplot.png')
 #        plt.close(f)
 #        sys.exit()
 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    
    #x = np.ma.log(BR_f_gas_Ha_ols__rg)
    x = BR_f_gas_Ha_ols__rg
    #x = Area_GAL__g[np.newaxis, :] * (H.McorSD__Trg[0] + BR_SigmaGas_Ha_ols__rg)
    #y = -1. * (H.O_O3N2_M13__rg - 8.69) / np.ma.log(BR_f_gas_Ha_ols__rg)
    #z = np.ma.log10(H.aSFRSD__Trg[0])
    #x = 1. - np.ma.exp(-1 * H.McorSD__Trg[0] / BR_SigmaGas_Ha_ols__rg)
    #x = np.ma.log10(H.McorSD__Trg[0] / BR_SigmaGas_Ha_ols__rg)
    #y = H.O_O3N2_M13__rg - 8.69
    y = H.O_O3N2_M13__rg
    z = H.Rtoplot()
    #x = 1. - np.ma.exp(-1 * H.Mcor_GAL__g / integrated_M_HI_miguel__g)
    #x = np.ma.log10(H.Mcor_GAL__g / integrated_M_HI_miguel__g)
    #y = H.integrated_O_O3N2_M13__g - 8.69
    #z = np.ma.log10(H.integrated_SFRSD__Tg[0]) 
    #sc_kwargs = dict(cmap = plt.cm.get_cmap('jet_r', 6), marker = 'o', s = 10, alpha = 0.9, edgecolor = 'none')
    sc_kwargs = dict(marker = 'o', s = 10, alpha = 0.9, edgecolor = 'none', cmap = 'Spectral')
    #z = H.reply_arr_by_radius(H.morfType_GAL__g), 
    #zticks = [9., 9.5, 10, 10.5, 11., 11.5], 
    #zticklabels = ['Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd'], 
     
    plt.clf()
    f = plt.gcf()
    ax = f.gca()
    r_kw = plot_zbins(return_kwargs = True, 
                      debug = True, 
                      f = f, 
                      ax = ax, 
                      x = x, 
                      y = y,
                      xlim = (0, 1),
                      ylim = (8.2, 8.8),
                      z = z, 
                      #z = H.reply_arr_by_radius(H.morfType_GAL__g), 
                      #zticks = [9., 9.5, 10, 10.5, 11., 11.5], 
                      #zticklabels = ['Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd'], 
                      #kwargs_scatter = dict(cmap = plt.cm.get_cmap('jet_r', 6), marker = 'o', s = 10, alpha = 0.9, edgecolor = 'none'),
                      kwargs_scatter = sc_kwargs, 
                      #add_mask = ~maskRadiusOk__rg,
                      #xlabel = r'1 - $e^{- \Sigma_{\star} / \Sigma_{gas}}',
                      #ylabel = r'$12 + \log(O/H) [Z_\odot]$', 
                       
                      #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                      # x_major_locator = .5, 
                      # x_minor_locator = .125, 
                      # y_major_locator = 1, 
                      # y_minor_locator = .25,
                      # xlim = (0, 2),
                      # ylim = (-3, 0), 
                      # xlabel = r'$\log \Sigma_{gas}\ [M_\odot\ pc^{-2}]$', 
                      # ylabel = r'$\log \Sigma_{SFR}\ [M_\odot\ kpc^{-2}\ yr^{-1}]$ ', 
                      #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                       
                      #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                      # xlabel = r'$\ln\ f_{gas}$ ', 
                      # ylabel = r'$\log[\frac{(O/H)}{(O/H)_\odot}]$',
                      # zlabel = r'$\log \Sigma_{SFR}(R, t_\star = 32Myrs)\ [M_\odot\ kpc^{-2}\ yr^{-1}]$', 
                      # xlim = (-7, 0), 
                      # ylim = (-0.5, 0), 
                      # x_major_locator = 2, 
                      # x_minor_locator = .4, 
                      # y_major_locator = .1, 
                      # y_minor_locator = .05,
                      #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

                      xlabel = r'$f_{gas}$ ',
                      ylabel = r'$12 + \log(O/H)$',
                      zlabel = r'R [HLR]', 
                      x_major_locator = .2, 
                      x_minor_locator = .05, 
                      y_major_locator = .2, 
                      y_minor_locator = .05,
                       
                      #ols = True, 
                      #kwargs_ols = dict(fs = 20, c = 'b'), 
                      #kwargs_ols_plot = dict(c = 'b', label = '')
                      )
    ax = r_kw['ax']
    ax.xaxis.labelpad = 20
    xm = r_kw['xm']
    ym = r_kw['ym']
    ax.set_title('N = %d elliptical bins' % len(xm.compressed()))
    #p = np.polyfit(xm.compressed(), ym.compressed(), 1)
    #ax.plot(ax.get_xlim(), p[0] * np.asarray(ax.get_xlim()) + p[1], '--k', lw = 2., label = '1d fit')
    #r = ym - (p[0] * xm + p[1])
    #plot_text_ax(ax, r'Z = %.2f $\ln f_{gas}$ - %.2f (rms: %.3f)' % (p[0], -1 * p[1], r.std()), fs = 20, pos_y = 0.05)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if args.output is None:
        if args.maskradius is not None:
            output = 'gas_maskRadius%.1f.pdf' % args.maskradius
        else:
            output = 'gas.pdf'
    else:
        output = args.output
          
    with PdfPages(output) as pdf:
        ##########################
        ######### PAGE 1 #########
        ##########################
        NRows = 3
        NCols = 2
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches[::-1]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax.axhline(y = RR_DGR, c = 'y')
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_ols__rg, mask = tmp_mask)
        rs_BR_DGR_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        m_aux = tmp_mask | BR_DGR_up_ols__rg.mask | BR_DGR_up_ols__rg.mask
        xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_up_ols__rg, mask = m_aux)
        rs_BR_DGR_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_down_ols__rg, mask = m_aux)
        rs_BR_DGR_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_DGR__rg, mask = tmp_mask) 
        rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_DGR_Ha__rg, mask = tmp_mask) 
        rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = 'b', **sc_kwargs)
            ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
            ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'g', **sc_kwargs)
        ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
        ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
        ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
        ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        ax.set_xlim(minR, H.RbinFin)
        #ax.set_ylim(0, .5)
        ax.set_title('Radial bins')
        ax.set_ylabel(r'$\delta_{DGR}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
          
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # if args.debug is not True:
        #     ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        #     ax.set_axis_on()
        #     tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_ols__g, mask = tmp_mask) 
        #     rs_BR_DGR_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     m_aux = tmp_mask | BR_DGR_up_ols__g.mask | BR_DGR_up_ols__g.mask 
        #     xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_up_ols__g, mask = m_aux)
        #     rs_BR_DGR_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_down_ols__g, mask = m_aux) 
        #     rs_BR_DGR_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_DGR__g, mask = tmp_mask) 
        #     rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_DGR_Ha__g, mask = tmp_mask) 
        #     rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     if args.scatter is True:
        #         ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = 'b', **sc_kwargs)
        #         ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
        #         ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'g', **sc_kwargs)
        #     ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
        #     ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
        #     ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
        #     ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        #     ax.axhline(y = RR_DGR, label = 'RR (const.)', c = 'y')
        #     plt.setp(ax.get_xticklabels(), visible = False)
        #     #ax.set_ylim(0, .5)
        #     ax.set_xlim(minR,  H.Rbin__r[-1])
        #     ax.tick_params(axis = 'y', which = 'major', pad = 10)
        #     ax.set_title('Zones')
        #     ax.yaxis.set_major_locator(MaxNLocator(6))
        #     ax.grid()
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_ols__rg), mask = tmp_mask) 
        rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_up_ols__rg), mask = tmp_mask) 
        rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_down_ols__rg), mask = yup.mask) 
        rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha_ols__rg), mask = tmp_mask)
        rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha_up_ols__rg), mask = tmp_mask) 
        rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(BR_SigmaGas_Ha_down_ols__rg), mask = yup.mask) 
        rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(SK_SigmaGas__rg), mask = tmp_mask)
        rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(SK_SigmaGas_Ha__rg), mask = tmp_mask)
        rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = np.ma.log10(RR_SigmaGas__rg), mask = tmp_mask) 
        rs_RR_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_SigmaGas_ols.x, rs_BR_SigmaGas_ols.y, c = 'b', **sc_kwargs)
            ax.scatter(rs_BR_SigmaGas_Ha_ols.x, rs_BR_SigmaGas_Ha_ols.y, c = 'r', **sc_kwargs)
            ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, c = 'g', **sc_kwargs)
            ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, c = 'black', **sc_kwargs)
            ax.scatter(rs_RR_SigmaGas.x, rs_RR_SigmaGas.y, c = 'y', **sc_kwargs)
        ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b', label = r'BR from $\tau_V^\star$')
        ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r', label = r'BR from $\tau_V^{neb}$')
        ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g', label = r'SK from synt.')
        ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black', label = r'SK from H$\alpha$')
        ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y', label = r'RR (DGR const.)')
        ax.fill_between(rs_BR_SigmaGas_up_ols.xS, rs_BR_SigmaGas_up_ols.yS, rs_BR_SigmaGas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        ax.fill_between(rs_BR_SigmaGas_Ha_up_ols.xS, rs_BR_SigmaGas_Ha_up_ols.yS, rs_BR_SigmaGas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(0.4, 2)
        #ax.legend(loc = 'upper right')
        ax.legend(bbox_to_anchor = (2.2, 2), fontsize = 12, frameon = False, ncol = 2)#, loc = 'upper right')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
         
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # if args.debug is not True:
        #     ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        #     ax.set_axis_on()
        #     tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_ols__g), mask = tmp_mask) 
        #     rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_up_ols__g), mask = tmp_mask) 
        #     rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_down_ols__g), mask = yup.mask) 
        #     rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_ols__g), mask = tmp_mask)
        #     rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_up_ols__g), mask = tmp_mask) 
        #     rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_down_ols__g), mask = yup.mask) 
        #     rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(SK_SigmaGas__g), mask = tmp_mask)
        #     rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(SK_SigmaGas_Ha__g), mask = tmp_mask)
        #     rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(RR_SigmaGas__g), mask = tmp_mask) 
        #     rs_RR_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     if args.scatter is True:
        #         ax.scatter(rs_BR_SigmaGas_ols.x, rs_BR_SigmaGas_ols.y, c = 'b', **sc_kwargs)
        #         ax.scatter(rs_BR_SigmaGas_Ha_ols.x, rs_BR_SigmaGas_Ha_ols.y, c = 'r', **sc_kwargs)
        #         ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, c = 'g', **sc_kwargs)
        #         ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, c = 'black', **sc_kwargs)
        #         ax.scatter(rs_RR_SigmaGas.x, rs_RR_SigmaGas.y, c = 'y', **sc_kwargs)
        #     ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b', label = r'BR from $\tau_V^\star$')
        #     ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r', label = r'BR from $\tau_V^{neb}$')
        #     ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g', label = r'SK from synt.')
        #     ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black', label = r'SK from H$\alpha$')
        #     ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y', label = r'RR (DGR const.)')
        #     ax.fill_between(rs_BR_SigmaGas_up_ols.xS, rs_BR_SigmaGas_up_ols.yS, rs_BR_SigmaGas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        #     ax.fill_between(rs_BR_SigmaGas_Ha_up_ols.xS, rs_BR_SigmaGas_Ha_up_ols.yS, rs_BR_SigmaGas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        #     ax.set_xlim(minR,  H.Rbin__r[-1])
        #     ax.set_ylim(0.4, 2)
        #     ax.legend(bbox_to_anchor = (2.6, 2), fontsize = 10, frameon = False, ncol = 2)#, loc = 'upper right')
        #     ax.tick_params(axis = 'y', which = 'major', pad = 15)
        #     plt.setp(ax.get_xticklabels(), visible = False)
        #     ax.yaxis.set_major_locator(MaxNLocator(6))
        #     ax.grid()
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         
        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_ols__rg, mask = tmp_mask) 
        rs_BR_f_gas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_up_ols__rg, mask = tmp_mask) 
        rs_BR_f_gas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_down_ols__rg, mask = yup.mask) 
        rs_BR_f_gas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha_ols__rg, mask = tmp_mask)
        rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, yup = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha_up_ols__rg, mask = tmp_mask) 
        rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        xm, ydown = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_f_gas_Ha_down_ols__rg, mask = yup.mask) 
        rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_f_gas__rg, mask = tmp_mask)
        rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = SK_f_gas_Ha__rg, mask = tmp_mask)
        rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = RR_f_gas__rg, mask = tmp_mask) 
        rs_RR_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        if args.scatter is True:
            ax.scatter(rs_BR_f_gas_ols.x, rs_BR_f_gas_ols.y, c = 'b', **sc_kwargs)
            ax.scatter(rs_BR_f_gas_Ha_ols.x, rs_BR_f_gas_Ha_ols.y, c = 'r', **sc_kwargs)
            ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, c = 'g', **sc_kwargs)
            ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, c = 'black', **sc_kwargs)
            ax.scatter(rs_RR_f_gas.x, rs_RR_f_gas.y, c = 'y', **sc_kwargs)
        ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'b')
        ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'r')
        ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
        ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black')
        ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'y')
        ax.fill_between(rs_BR_f_gas_up_ols.xS, rs_BR_f_gas_up_ols.yS, rs_BR_f_gas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        ax.fill_between(rs_BR_f_gas_Ha_up_ols.xS, rs_BR_f_gas_Ha_up_ols.yS, rs_BR_f_gas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(0, .4)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylabel(r'f${}_{gas}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_yticklabels(), visible = False)
                     
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # if args.debug is not True:
        #     ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        #     ax.set_axis_on()
        #     tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_ols__g, mask = ~maskRadiusOk__g) 
        #     rs_BR_f_gas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_up_ols__g, mask = tmp_mask) 
        #     rs_BR_f_gas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_down_ols__g, mask = yup.mask) 
        #     rs_BR_f_gas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_ols__g, mask = tmp_mask)
        #     rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_up_ols__g, mask = tmp_mask) 
        #     rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_down_ols__g, mask = yup.mask) 
        #     rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_f_gas__g, mask = tmp_mask)
        #     rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_f_gas_Ha__g, mask = tmp_mask)
        #     rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = RR_f_gas__g, mask = tmp_mask) 
        #     rs_RR_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        #     if args.scatter is True:
        #         ax.scatter(rs_BR_f_gas_ols.x, rs_BR_f_gas_ols.y, c = 'b', **sc_kwargs)
        #         ax.scatter(rs_BR_f_gas_Ha_ols.x, rs_BR_f_gas_Ha_ols.y, c = 'r', **sc_kwargs)
        #         ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, c = 'g', **sc_kwargs)
        #         ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, c = 'black', **sc_kwargs)
        #         ax.scatter(rs_RR_f_gas.x, rs_RR_f_gas.y, c = 'y', **sc_kwargs)
        #     ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'b')
        #     ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'r')
        #     ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
        #     ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black')
        #     ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'y')
        #     ax.fill_between(rs_BR_f_gas_up_ols.xS, rs_BR_f_gas_up_ols.yS, rs_BR_f_gas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        #     ax.fill_between(rs_BR_f_gas_Ha_up_ols.xS, rs_BR_f_gas_Ha_up_ols.yS, rs_BR_f_gas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        #     ax.set_xlim(minR,  H.Rbin__r[-1])
        #     ax.set_ylim(0, .4)
        #     ax.tick_params(axis = 'y', which = 'major', pad = 15)
        #     ax.set_xlabel(r'R [HLR]')
        #     ax.yaxis.set_major_locator(MaxNLocator(6))
        #     ax.grid()
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         
        ax = plt.subplot2grid(grid_shape, loc = (0, NCols - 1), rowspan = NRows)
        ax.set_axis_off()
        txt = r'NGals:%d  tSF:%.2f Myr' % (N_gals, (H.tSF__T[iT] / 1e6))
        kw_text = dict(pos_x = 0, pos_y = 0.65, fs = 12, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        kw_text = dict(pos_x = 0, pos_y = 0.60, fs = 12, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$\delta_{DGR} = (%.2f\times 10^{-3} - %.2f\times 10^{-2}) (\frac{O/H}{(O/H)_\odot})$' % (DGR_conv_lim_inf / 1e-3, DGR_conv_lim_sup / 1e-2)
        kw_text = dict(pos_x = 0, pos_y = 0.53, fs = 12, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$f_{gas}\ =\ [1 + (\frac{\sigma_d}{m_d})(\delta_{DGR})(\frac{\mu_\star}{\tau_V})]^{-1}$'
        kw_text = dict(pos_x = 0, pos_y = 0.47, fs = 15, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        f.subplots_adjust(hspace = 0.2, wspace = 0.3)
        pdf.savefig(f)
        plt.close(f)


        ##########################
        ######### PAGE 2 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        #page_size_inches = (NCols * 3, NRows * 2.5)
        page_size_inches = A4Size_inches[::-1]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        yaxis = [ 'alogSFRSDHakpcR', 'alogSFRSDkpcR', 'atfluxR', 'xYR', 'logO3N2M13R', 'alogZmassR', 'logMcorSDR' ]
        row, col = 0, 0
        for yk in yaxis:
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            _, ydict = H.get_plot_dict(iT, -1, yk)
            x = H.Rtoplot()
            y = ydict['v']
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            xm, ym = C.ma_mask_xyz(x = x, y = y, mask = tmp_mask)
            #rs_kwargs['sigma'] = 1.2
            rs_kwargs['overlap'] = 0.1
            rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            #if args.scatter is True:
            ax.scatter(rs.x, rs.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.plot(rs.xS, rs.yS, '.-', c = 'k')
            ax.set_xlim(minR, H.RbinFin)
            ax.set_ylim(ydict['lim'])
            ax.set_xlabel(r'R [HLR]')
            ax.set_ylabel(ydict['label'])
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            #plt.setp(ax.get_yticklabels(), visible = False)
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)

        
        ##########################
        ######### PAGE 3 #########
        ##########################
        NRows = 2
        NCols = 3
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        #page_size_inches = (NCols * 3, NRows * 2.5)
        page_size_inches = A4Size_inches[::-1]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        xaxis = [ 'atfluxR', 'xYR', 'morfTypeR', 'baR', 'alogZmassR', 'logMcorSDR' ]
        #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR' ]
        row, col = 0, 0
        for xk in xaxis:
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            _, xdict = H.get_plot_dict(iT, -1, xk)
            x = xdict['v']
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            xm, ym_BR_DGR_ols = C.ma_mask_xyz(x = x, y = BR_DGR_ols__rg, mask = tmp_mask)
            rs_BR_DGR_ols = runstats(xm.compressed(), ym_BR_DGR_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_DGR_up_ols = C.ma_mask_xyz(x = x, y = BR_DGR_up_ols__rg, mask = tmp_mask)
            rs_BR_DGR_up_ols = runstats(xm.compressed(), ym_BR_DGR_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_DGR_down_ols = C.ma_mask_xyz(x = x, y = BR_DGR_down_ols__rg, mask = tmp_mask)
            rs_BR_DGR_down_ols = runstats(xm.compressed(), ym_BR_DGR_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_DGR = C.ma_mask_xyz(x = x, y = SK_DGR__rg, mask = tmp_mask)
            rs_SK_DGR = runstats(xm.compressed(), ym_SK_DGR.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_DGR_Ha = C.ma_mask_xyz(x = x, y = SK_DGR_Ha__rg, mask = tmp_mask)
            rs_SK_DGR_Ha = runstats(xm.compressed(), ym_SK_DGR_Ha.compressed(), nBox = 20, **rs_kwargs)
            if args.scatter is True:
                ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
            ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(xdict['label'])
            ax.set_ylim(0, 15e-3)
            if col == 0:
                ax.set_ylabel(r'$\delta_{DGR}$')
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)

        ##########################
        ######### PAGE 4 #########
        ##########################
        NRows = 2
        NCols = 3
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        #page_size_inches = (NCols * 3, NRows * 2.5)
        page_size_inches = A4Size_inches[::-1]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        xaxis = [ 'atfluxR', 'xYR', 'morfTypeR', 'baR', 'alogZmassR', 'logMcorSDR', ]
        #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR', ]
        row, col = 0, 0
        for xk in xaxis:
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            _, xdict = H.get_plot_dict(iT, -1, xk)
            x = xdict['v']
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            xm, ym_BR_SigmaGas_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym_BR_SigmaGas_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_up_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_up_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), ym_BR_SigmaGas_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_down_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_down_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ym_BR_SigmaGas_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_Ha_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_Ha_up_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_up_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_Ha_down_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_down_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas__rg), mask = tmp_mask)
            rs_SK_SigmaGas = runstats(xm.compressed(), ym_SK_SigmaGas.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_SigmaGas_Ha = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas_Ha__rg), mask = tmp_mask)
            rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym_SK_SigmaGas_Ha.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_RR_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(RR_SigmaGas__rg), mask = tmp_mask)
            rs_RR_SigmaGas = runstats(xm.compressed(), ym_RR_SigmaGas.compressed(), nBox = 20, **rs_kwargs)
            if args.scatter is True:
                ax.scatter(rs_BR_SigmaGas_ols.x, rs_BR_SigmaGas_ols.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_BR_SigmaGas_Ha_ols.x, rs_BR_SigmaGas_Ha_ols.y, marker = 'o', c = 'r', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_RR_SigmaGas.x, rs_RR_SigmaGas.y, marker = 'o', c = 'y', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b')
            ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r')
            ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
            ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black')
            ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y')
            ax.fill_between(rs_BR_SigmaGas_up_ols.xS, rs_BR_SigmaGas_up_ols.yS, rs_BR_SigmaGas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.fill_between(rs_BR_SigmaGas_Ha_up_ols.xS, rs_BR_SigmaGas_Ha_up_ols.yS, rs_BR_SigmaGas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)            
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(xdict['label'])
            ax.set_ylim(0.4, 2)
            #ax.set_ylim(0,15e-3)
            if col == 0:
                ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)

        ##########################
        ######### PAGE 5 #########
        ##########################
        NRows = 2
        NCols = 3
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        #page_size_inches = (NCols * 3, NRows * 2.5)
        page_size_inches = A4Size_inches[::-1]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        xaxis = [ 'atfluxR', 'xYR', 'morfTypeR', 'baR', 'alogZmassR', 'logMcorSDR' ]
        #xaxis = [ 'atfluxR', 'xYR', 'baR', 'alogZmassR', 'logMcorSDR' ]
        row, col = 0, 0
        for xk in xaxis:
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            _, xdict = H.get_plot_dict(iT, -1, xk)
            x = xdict['v']
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            xm, ym_BR_f_gas_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_ols = runstats(xm.compressed(), ym_BR_f_gas_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_up_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_up_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_up_ols = runstats(xm.compressed(), ym_BR_f_gas_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_down_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_down_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_down_ols = runstats(xm.compressed(), ym_BR_f_gas_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_Ha_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_Ha_up_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_up_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_Ha_down_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_down_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_f_gas = C.ma_mask_xyz(x = x, y = SK_f_gas__rg, mask = tmp_mask)
            rs_SK_f_gas = runstats(xm.compressed(), ym_SK_f_gas.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_f_gas_Ha = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__rg, mask = tmp_mask)
            rs_SK_f_gas_Ha = runstats(xm.compressed(), ym_SK_f_gas_Ha.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_RR_f_gas = C.ma_mask_xyz(x = x, y = RR_f_gas__rg, mask = tmp_mask)
            rs_RR_f_gas = runstats(xm.compressed(), ym_RR_f_gas.compressed(), nBox = 20, **rs_kwargs)
            if args.scatter is True:
                ax.scatter(rs_BR_f_gas_ols.x, rs_BR_f_gas_ols.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_BR_f_gas_Ha_ols.x, rs_BR_f_gas_Ha_ols.y, marker = 'o', c = 'r', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
                ax.scatter(rs_RR_f_gas.x, rs_RR_f_gas.y, marker = 'o', c = 'y', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'b')
            ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'r')
            ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
            ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black')
            ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'y')
            ax.fill_between(rs_BR_f_gas_up_ols.xS, rs_BR_f_gas_up_ols.yS, rs_BR_f_gas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.fill_between(rs_BR_f_gas_Ha_up_ols.xS, rs_BR_f_gas_Ha_up_ols.yS, rs_BR_f_gas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)            
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(xdict['label'])
            ax.set_ylim(0, 0.4)
            #ax.set_ylim(0,15e-3)
            if col == 0:
                ax.set_ylabel(r'f${}_{gas}$')
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
