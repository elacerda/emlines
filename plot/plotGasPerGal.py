#!/usr/bin/python
#
# Lacerda@Granada - 1/Jul/2015
#
from CALIFAUtils.plots import DrawHLRCircleInSDSSImage
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
from CALIFAUtils.plots import DrawHLRCircle
from CALIFAUtils.plots import plot_text_ax
from matplotlib.ticker import MaxNLocator
#from CALIFAUtils.plots import plot_zbins
from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

RNuc = 0.5

mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'scatter' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'maskradius' : None,
        'califaID' : 'K0073',
        'vrunstr' : 'v20_q050.d15a',
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--scatter', '-s',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--hdf5', '-H',
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
    parser.add_argument('--califaID', '-g',
                        metavar = 'FILE',
                        type = str,
                        default = default['califaID'])
    parser.add_argument('--vrunstr',
                        metavar = 'RUNCODE',
                        type = str,
                        default = default['vrunstr'])

    return parser.parse_args()

def plot_ma_map(ax, map, p, draw_HLR = False, K = None, **kwargs_imshow):
    vmin, vmax = np.percentile(map.compressed(), p)
    kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = mpl.cm.spectral, vmin = vmin, vmax = vmax)
    kwargs.update(kwargs_imshow)
    im = ax.imshow(map, **kwargs)
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    if (draw_HLR == True) and K is not None:
        DrawHLRCircle(ax, K)
    return im

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    paths = C.CALIFAPaths()
    
    H = C.H5SFRData(args.hdf5)
    iT = args.itSF
    iU = -1

    minR = 0
    
    if (len(np.where(H.califaIDs == args.califaID)[0]) == 0):
        exit('<<< plot: %s: no data.' % args.califaID)
        
    K = C.read_one_cube(args.califaID, EL = True, v_run = args.vrunstr)
    galaxyImgFile = paths.get_image_file(args.califaID)
    
    # Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
        
    iGal = H.califaIDs.tolist().index(args.califaID)
    iGal_all = H.califaIDs_all.tolist().index(args.califaID)

    # Integrated or single-value measurements or calculations
    Mcor_tot = H.Mcor_GAL__g[iGal_all]
    McorSD = H.McorSD_GAL__g[iGal_all]
    mtype = H.morfType_GAL__g[iGal_all]
    N_zone = H.N_zones__g[iGal]
    at_flux = H.at_flux_GAL__g[iGal_all]
    ba = H.ba_GAL__g[iGal_all]
    ba_pyc = H.ba_PyCASSO_GAL__g[iGal_all]
    tau_V = H.integrated_tau_V__g[iGal_all]
    tau_V_neb = H.integrated_tau_V_neb__g[iGal_all]
    SFR__T = H.integrated_SFR__Tg[:, iGal_all]
    SFRSD__T = H.integrated_SFR__Tg[:, iGal_all] 
    SFR_Ha = H.integrated_SFR_Ha__g[iGal_all]
    SFRSD_Ha = H.integrated_SFR_Ha__g[iGal_all] 
        
    zone_dist_HLR__z = H.get_prop_gal(H.zone_dist_HLR__g, args.califaID)
    zone_area_pc2__z = H.get_prop_gal(H.zone_area_pc2__g, args.califaID)

    ## Synt.
    #### zones
    Mcor__z = H.get_prop_gal(H.Mcor__g, args.califaID)
    Mcor__Tz = H.get_prop_gal(H.Mcor__Tg, args.califaID)
    McorSD__z = H.get_prop_gal(H.McorSD__g, args.califaID)
    McorSD__Tz = H.get_prop_gal(H.McorSD__Tg, args.califaID)
    SFR__Tz = H.get_prop_gal(H.SFR__Tg, args.califaID)
    SFRSD__Tz = H.get_prop_gal(H.SFRSD__Tg, args.califaID)
    tau_V__Tz = H.get_prop_gal(H.tau_V__Tg, args.califaID)
    alogZ_mass__Uz = H.get_prop_gal(H.alogZ_mass__Ug, args.califaID)
    at_flux__Tz = H.get_prop_gal(H.at_flux__Tg, args.califaID)
    x_Y__Tz = H.get_prop_gal(H.x_Y__Tg, args.califaID)
    #### Radial averages
    McorSD__r = H.get_prop_gal(H.McorSD__Tg, args.califaID)
    McorSD__Tr = H.get_prop_gal(H.McorSD__Trg, args.califaID)
    tau_V__Tr = H.get_prop_gal(H.tau_V__Trg, args.califaID)
    aSFRSD__Tr = H.get_prop_gal(H.aSFRSD__Trg, args.califaID)
    alogZ_mass__Ur = H.get_prop_gal(H.alogZ_mass__Urg, args.califaID)
    at_flux__Tr = H.get_prop_gal(H.at_flux__Trg, args.califaID)
    x_Y__Tr = H.get_prop_gal(H.x_Y__Trg, args.califaID)

    ## Nebul.
    #### zones
    SFR_Ha__z = H.get_prop_gal(H.SFR_Ha__g, args.califaID)
    SFRSD_Ha__z = H.get_prop_gal(H.SFRSD_Ha__g, args.califaID)
    tau_V_neb__z = H.get_prop_gal(H.tau_V_neb__g, args.califaID)
    O_O3N2_M13__z = H.get_prop_gal(H.O_O3N2_M13__g, args.califaID)
    #### Radial averages
    aSFRSD_Ha__r = H.get_prop_gal(H.aSFRSD_Ha__rg, args.califaID)
    tau_V_neb__r = H.get_prop_gal(H.tau_V_neb__rg, args.califaID)
    O_O3N2_M13__r = H.get_prop_gal(H.O_O3N2_M13__rg, args.califaID)
    
    if args.maskradius is None:
        maskRadiusOk__z = np.ones(K.N_zone, dtype = np.bool)
        maskRadiusOk__r = np.ones(H.NRbins, dtype = np.bool)
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
        maskRadiusOk__z = np.greater_equal(K.zoneDistance_HLR, minR)
        maskRadiusOk__r = np.greater_equal(H.RbinCenter__r, minR)
        maskRadiusOk__g = (H.zone_dist_HLR__g >= args.maskradius) & (H.zone_dist_HLR__g <= H.Rbin__r[-1]) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= args.maskradius)).T
    
    dustdim = 0.2 # md / rhod

    ######################
    # Galaxy values
    ######################
    # SK Law 
    ######################
    SK_zero = 1.6e-4
    SK_slope = 1.4
    aux = 1e6 ** (1. / SK_slope)
    SK_SigmaGas__z = aux * (SFRSD__Tz[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__z = aux * (SFRSD_Ha__z / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas__r = aux * (aSFRSD__Tr[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__r = aux * (aSFRSD_Ha__r / SK_zero) ** (1. / SK_slope)
    SK_DGR__z = dustdim * tau_V__Tz[iT] / SK_SigmaGas__z
    SK_DGR_Ha__z = dustdim * tau_V_neb__z / SK_SigmaGas_Ha__z
    SK_DGR__r = dustdim * tau_V__Tr[iT] / SK_SigmaGas__r
    SK_DGR_Ha__r = dustdim * tau_V_neb__r / SK_SigmaGas_Ha__r
    SK_GSR__z = SK_SigmaGas__z / McorSD__Tz[iT]
    SK_GSR_Ha__z = SK_SigmaGas_Ha__z / McorSD__Tz[iT]
    SK_GSR__r = SK_SigmaGas__r / McorSD__Tr[iT]
    SK_GSR_Ha__r = SK_SigmaGas_Ha__r / McorSD__Tr[iT]
    SK_f_gas__z = 1. / (1. + 1. / SK_GSR__z)
    SK_f_gas_Ha__z = 1. / (1. + 1. / SK_GSR_Ha__z)
    SK_f_gas__r = 1. / (1. + 1. / SK_GSR__r)
    SK_f_gas_Ha__r = 1. / (1. + 1. / SK_GSR_Ha__r)
    ######################
    
    ######################
    # Remy-Ruyer
    ######################
    RR_DGR = 10. ** (-2.21) # 0.006166
    RR_SigmaGas__z = dustdim * tau_V__Tz[iT] / RR_DGR
    RR_SigmaGas_Ha__z = dustdim * tau_V_neb__z / RR_DGR
    RR_SigmaGas__r = dustdim * tau_V__Tr[iT] / RR_DGR
    RR_SigmaGas_Ha__r = dustdim * tau_V_neb__r / RR_DGR
    RR_GSR__z = RR_SigmaGas__z / McorSD__Tz[iT]
    RR_GSR_Ha__z = RR_SigmaGas_Ha__z / McorSD__Tz[iT]
    RR_GSR__r = RR_SigmaGas__r / McorSD__Tr[iT]
    RR_GSR_Ha__r = RR_SigmaGas_Ha__r / McorSD__Tr[iT]
    RR_f_gas__z = 1. / (1. + 1. / RR_GSR__z)
    RR_f_gas_Ha__z = 1. / (1. + 1. / RR_GSR_Ha__z)
    RR_f_gas__r = 1. / (1. + 1. / RR_GSR__r)
    RR_f_gas_Ha__r = 1. / (1. + 1. / RR_GSR_Ha__r)
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
    BR_OHBrinch_ols__z = np.ma.masked_all((O_O3N2_M13__z.shape))
    BR_OHBrinch_ols__z[~O_O3N2_M13__z.mask] = np.polyval(p_ols, O_O3N2_M13__z.compressed()) 
    BR_OHBrinch_ols__r = np.ma.masked_all((O_O3N2_M13__r.shape))
    BR_OHBrinch_ols__r[~O_O3N2_M13__r.mask] = np.polyval(p_ols, O_O3N2_M13__r.compressed()) 
    BR_DGR_up_ols__z = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__z - 12) * OHSunBrinch_inv)
    BR_DGR_up_ols__r = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__r - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__z = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__z - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__r = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__r - 12) * OHSunBrinch_inv)
    BR_DGR_ols__z = DGR_cte * (10 ** (BR_OHBrinch_ols__z - 12) * OHSunBrinch_inv)
    BR_DGR_ols__r = DGR_cte * (10 ** (BR_OHBrinch_ols__r - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_ols__z = dustdim * tau_V__Tz[iT] / BR_DGR_up_ols__z
    BR_SigmaGas_up_ols__r = dustdim * tau_V__Tr[iT] / BR_DGR_up_ols__r
    BR_SigmaGas_Ha_up_ols__z = dustdim * tau_V_neb__z / BR_DGR_up_ols__z
    BR_SigmaGas_Ha_up_ols__r = dustdim * tau_V_neb__r / BR_DGR_up_ols__r
    BR_SigmaGas_down_ols__z = dustdim * tau_V__Tz[iT] / BR_DGR_down_ols__z
    BR_SigmaGas_down_ols__r = dustdim * tau_V__Tr[iT] / BR_DGR_down_ols__r
    BR_SigmaGas_Ha_down_ols__z = dustdim * tau_V_neb__z / BR_DGR_down_ols__z
    BR_SigmaGas_Ha_down_ols__r = dustdim * tau_V_neb__r / BR_DGR_down_ols__r
    BR_SigmaGas_ols__z = dustdim * tau_V__Tz[iT] / BR_DGR_ols__z
    BR_SigmaGas_ols__r = dustdim * tau_V__Tr[iT] / BR_DGR_ols__r
    BR_SigmaGas_Ha_ols__z = dustdim * tau_V_neb__z / BR_DGR_ols__z
    BR_SigmaGas_Ha_ols__r = dustdim * tau_V_neb__r / BR_DGR_ols__r
    BR_GSR_up_ols__z = BR_SigmaGas_up_ols__z / McorSD__Tz[iT]
    BR_GSR_up_ols__r = BR_SigmaGas_up_ols__r / McorSD__Tr[iT]
    BR_GSR_Ha_up_ols__z = BR_SigmaGas_Ha_up_ols__z / McorSD__Tz[iT]
    BR_GSR_Ha_up_ols__r = BR_SigmaGas_Ha_up_ols__r / McorSD__Tr[iT]
    BR_GSR_down_ols__z = BR_SigmaGas_down_ols__z / McorSD__Tz[iT]
    BR_GSR_down_ols__r = BR_SigmaGas_down_ols__r / McorSD__Tr[iT]
    BR_GSR_Ha_down_ols__z = BR_SigmaGas_Ha_down_ols__z / McorSD__Tz[iT]
    BR_GSR_Ha_down_ols__r = BR_SigmaGas_Ha_down_ols__r / McorSD__Tr[iT]
    BR_GSR_ols__z = BR_SigmaGas_ols__z / McorSD__Tz[iT]
    BR_GSR_ols__r = BR_SigmaGas_ols__r / McorSD__Tr[iT]
    BR_GSR_Ha_ols__z = BR_SigmaGas_Ha_ols__z / McorSD__Tz[iT]
    BR_GSR_Ha_ols__r = BR_SigmaGas_Ha_ols__r / McorSD__Tr[iT]
    BR_f_gas_up_ols__z = 1. / (1. + 1. / BR_GSR_up_ols__z)
    BR_f_gas_up_ols__r = 1. / (1. + 1. / BR_GSR_up_ols__r)
    BR_f_gas_Ha_up_ols__z = 1. / (1. + 1. / BR_GSR_Ha_up_ols__z)
    BR_f_gas_Ha_up_ols__r = 1. / (1. + 1. / BR_GSR_Ha_up_ols__r)
    BR_f_gas_down_ols__z = 1. / (1. + 1. / BR_GSR_down_ols__z)
    BR_f_gas_down_ols__r = 1. / (1. + 1. / BR_GSR_down_ols__r)
    BR_f_gas_Ha_down_ols__z = 1. / (1. + 1. / BR_GSR_Ha_down_ols__z)
    BR_f_gas_Ha_down_ols__r = 1. / (1. + 1. / BR_GSR_Ha_down_ols__r)
    BR_f_gas_ols__z = 1. / (1. + 1. / BR_GSR_ols__z)
    BR_f_gas_ols__r = 1. / (1. + 1. / BR_GSR_ols__r)
    BR_f_gas_Ha_ols__z = 1. / (1. + 1. / BR_GSR_Ha_ols__z)
    BR_f_gas_Ha_ols__r = 1. / (1. + 1. / BR_GSR_Ha_ols__r)
    ######################
    # from cubic polynomial fit
    p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    BR_OHBrinch_cubic__z = np.ma.masked_all((O_O3N2_M13__z.shape))
    BR_OHBrinch_cubic__z[~O_O3N2_M13__z.mask] = np.polyval(p_cubic, O_O3N2_M13__z.compressed()) 
    BR_OHBrinch_cubic__r = np.ma.masked_all((O_O3N2_M13__r.shape))
    BR_OHBrinch_cubic__r[~O_O3N2_M13__r.mask] = np.polyval(p_cubic, O_O3N2_M13__r.compressed()) 
    BR_DGR_up_cubic__z = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__z - 12) * OHSunBrinch_inv)
    BR_DGR_up_cubic__r = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__r - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__z = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__z - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__r = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__r - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__z = DGR_cte * (10 ** (BR_OHBrinch_cubic__z - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__r = DGR_cte * (10 ** (BR_OHBrinch_cubic__r - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_cubic__z = dustdim * tau_V__Tz[iT] / BR_DGR_up_cubic__z
    BR_SigmaGas_up_cubic__r = dustdim * tau_V__Tr[iT] / BR_DGR_up_cubic__r
    BR_SigmaGas_Ha_up_cubic__z = dustdim * tau_V_neb__z / BR_DGR_up_cubic__z
    BR_SigmaGas_Ha_up_cubic__r = dustdim * tau_V_neb__r / BR_DGR_up_cubic__r
    BR_SigmaGas_down_cubic__z = dustdim * tau_V__Tz[iT] / BR_DGR_down_cubic__z
    BR_SigmaGas_down_cubic__r = dustdim * tau_V__Tr[iT] / BR_DGR_down_cubic__r
    BR_SigmaGas_Ha_down_cubic__z = dustdim * tau_V_neb__z / BR_DGR_down_cubic__z
    BR_SigmaGas_Ha_down_cubic__r = dustdim * tau_V_neb__r / BR_DGR_down_cubic__r
    BR_SigmaGas_cubic__z = dustdim * tau_V__Tz[iT] / BR_DGR_cubic__z
    BR_SigmaGas_cubic__r = dustdim * tau_V__Tr[iT] / BR_DGR_cubic__r
    BR_SigmaGas_Ha_cubic__z = dustdim * tau_V_neb__z / BR_DGR_cubic__z
    BR_SigmaGas_Ha_cubic__r = dustdim * tau_V_neb__r / BR_DGR_cubic__r
    BR_GSR_up_cubic__z = BR_SigmaGas_up_cubic__z / McorSD__Tz[iT]
    BR_GSR_up_cubic__r = BR_SigmaGas_up_cubic__r / McorSD__Tr[iT]
    BR_GSR_Ha_up_cubic__z = BR_SigmaGas_Ha_up_cubic__z / McorSD__Tz[iT]
    BR_GSR_Ha_up_cubic__r = BR_SigmaGas_Ha_up_cubic__r / McorSD__Tr[iT]
    BR_GSR_down_cubic__z = BR_SigmaGas_down_cubic__z / McorSD__Tz[iT]
    BR_GSR_down_cubic__r = BR_SigmaGas_down_cubic__r / McorSD__Tr[iT]
    BR_GSR_Ha_down_cubic__z = BR_SigmaGas_Ha_down_cubic__z / McorSD__Tz[iT]
    BR_GSR_Ha_down_cubic__r = BR_SigmaGas_Ha_down_cubic__r / McorSD__Tr[iT]
    BR_GSR_cubic__z = BR_SigmaGas_cubic__z / McorSD__Tz[iT]
    BR_GSR_cubic__r = BR_SigmaGas_cubic__r / McorSD__Tr[iT]
    BR_GSR_Ha_cubic__z = BR_SigmaGas_Ha_cubic__z / McorSD__Tz[iT]
    BR_GSR_Ha_cubic__r = BR_SigmaGas_Ha_cubic__r / McorSD__Tr[iT]
    BR_f_gas_up_cubic__z = 1. / (1. + 1. / BR_GSR_up_cubic__z)
    BR_f_gas_up_cubic__r = 1. / (1. + 1. / BR_GSR_up_cubic__r)
    BR_f_gas_Ha_up_cubic__z = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__z)
    BR_f_gas_Ha_up_cubic__r = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__r)
    BR_f_gas_down_cubic__z = 1. / (1. + 1. / BR_GSR_down_cubic__z)
    BR_f_gas_down_cubic__r = 1. / (1. + 1. / BR_GSR_down_cubic__r)
    BR_f_gas_Ha_down_cubic__z = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__z)
    BR_f_gas_Ha_down_cubic__r = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__r)
    BR_f_gas_cubic__z = 1. / (1. + 1. / BR_GSR_cubic__z)
    BR_f_gas_cubic__r = 1. / (1. + 1. / BR_GSR_cubic__r)
    BR_f_gas_Ha_cubic__z = 1. / (1. + 1. / BR_GSR_Ha_cubic__z)
    BR_f_gas_Ha_cubic__r = 1. / (1. + 1. / BR_GSR_Ha_cubic__r)
    ######################

    ######################
    # Global values
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
    SK_integrated_SigmaGas_Ha = aux * (H.integrated_SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope) 
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
    ######################
    # from cubic polynomial fit
    p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    BR_OHBrinch_cubic__g = np.ma.masked_all((H.O_O3N2_M13__g.shape))
    BR_OHBrinch_cubic__g[~H.O_O3N2_M13__g.mask] = np.polyval(p_cubic, H.O_O3N2_M13__g.compressed()) 
    BR_OHBrinch_cubic__rg = np.ma.masked_all((H.O_O3N2_M13__rg.shape))
    BR_OHBrinch_cubic__rg[~H.O_O3N2_M13__rg.mask] = np.polyval(p_cubic, H.O_O3N2_M13__rg.compressed()) 
    BR_DGR_up_cubic__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_cubic__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__g = DGR_cte * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__rg = DGR_cte * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_up_cubic__g
    BR_SigmaGas_up_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_up_cubic__rg
    BR_SigmaGas_Ha_up_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_up_cubic__g
    BR_SigmaGas_Ha_up_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_up_cubic__rg
    BR_SigmaGas_down_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_down_cubic__g
    BR_SigmaGas_down_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_down_cubic__rg
    BR_SigmaGas_Ha_down_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_down_cubic__g
    BR_SigmaGas_Ha_down_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_down_cubic__rg
    BR_SigmaGas_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_cubic__g
    BR_SigmaGas_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_cubic__rg
    BR_SigmaGas_Ha_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_cubic__g
    BR_SigmaGas_Ha_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_cubic__rg
    BR_GSR_up_cubic__g = BR_SigmaGas_up_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_up_cubic__rg = BR_SigmaGas_up_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_up_cubic__g = BR_SigmaGas_Ha_up_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_up_cubic__rg = BR_SigmaGas_Ha_up_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_down_cubic__g = BR_SigmaGas_down_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_down_cubic__rg = BR_SigmaGas_down_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_down_cubic__g = BR_SigmaGas_Ha_down_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_down_cubic__rg = BR_SigmaGas_Ha_down_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_cubic__g = BR_SigmaGas_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_cubic__rg = BR_SigmaGas_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_cubic__g = BR_SigmaGas_Ha_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_cubic__rg = BR_SigmaGas_Ha_cubic__rg / H.McorSD__Trg[iT]
    BR_f_gas_up_cubic__g = 1. / (1. + 1. / BR_GSR_up_cubic__g)
    BR_f_gas_up_cubic__rg = 1. / (1. + 1. / BR_GSR_up_cubic__rg)
    BR_f_gas_Ha_up_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__g)
    BR_f_gas_Ha_up_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__rg)
    BR_f_gas_down_cubic__g = 1. / (1. + 1. / BR_GSR_down_cubic__g)
    BR_f_gas_down_cubic__rg = 1. / (1. + 1. / BR_GSR_down_cubic__rg)
    BR_f_gas_Ha_down_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__g)
    BR_f_gas_Ha_down_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__rg)
    BR_f_gas_cubic__g = 1. / (1. + 1. / BR_GSR_cubic__g)
    BR_f_gas_cubic__rg = 1. / (1. + 1. / BR_GSR_cubic__rg)
    BR_f_gas_Ha_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_cubic__g)
    BR_f_gas_Ha_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_cubic__rg)
    ######################
    
    if args.output is None:
        if args.maskradius is not None:
            output = '%s_gas_maskRadius%.1f.pdf' % (K.califaID, args.maskradius)
        else:
            output = '%s_gas.pdf' % K.califaID
    else:
        output = args.output
          
    default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)
    default_im_kwargs = dict(interpolation = 'nearest', origin = 'lower', aspect = 'auto', cmap = 'soectral')
    
    phys__z = [
        tau_V__Tz[iT], 
        tau_V_neb__z,
        SFRSD__Tz[iT] * 1e6, 
        SFRSD_Ha__z * 1e6,
        McorSD__Tz[iT], 
        O_O3N2_M13__z, 
        alogZ_mass__Uz[-1],
    ]
    
    phys__yx = [
        K.zoneToYX(prop, extensive = False, surface_density = False) 
            for prop in phys__z
    ]
    
    phys__r = [
        np.ma.log10(tau_V__Tr[iT]), 
        np.ma.log10(tau_V_neb__r),
        np.ma.log10(aSFRSD__Tr[iT] * 1e6), 
        np.ma.log10(aSFRSD_Ha__r * 1e6),
        np.ma.log10(McorSD__Tr[iT]), 
        O_O3N2_M13__r, 
        alogZ_mass__Ur[-1],
    ]
    
    labels__r = [
        r'$\log\ \tau_V^\star (R)$',
        r'$\log\ \tau_V^{neb} (R)$',
        r'$\log\ \Sigma_{SFR}^\star (t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$',
        r'$\log\ \Sigma_{SFR}^{neb} (R)\ [M_\odot yr^{-1} kpc^{-2}]$',
        r'$\log\ \mu_\star (R)$ [$M_\odot \ pc^{-2}$]',
        r'12 + $\log\ O/H\ (R)$',
        r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9),
    ]

    labels__z = [
        r'$\log\ \tau_V^\star$',
        r'$\log\ \tau_V^{neb}$',
        r'$\log\ \Sigma_{SFR}^\star (t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$',
        r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$',
        r'$\log\ \mu_\star$ [$M_\odot \ pc^{-2}$]',
        r'12 + $\log\ O/H$',
        r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[-1] / 1e9),
    ]
    
    if args.debug: 
        sys.exit('bai')

    mtype_str = ['E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'S0', 'S0a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm', 'Sm', 'Ir']
    mtype_num = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 14])

    with PdfPages(output) as pdf:
        #########################
        ######## PAGE 1 #########
        ##########################
        NRows = 4
        NCols = 4
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        page_size_inches = (NCols * 3, NRows * 2)
        f.set_size_inches(page_size_inches)
        for ax in f.axes:
            ax.set_axis_off()
        grid_shape = (NRows, NCols)
        ax = plt.subplot2grid(grid_shape, loc = (0, 0), rowspan = 2, colspan = 2)
        ax.set_axis_on()
        galimg = plt.imread(galaxyImgFile)[::-1, :, :]
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.imshow(galimg, origin = 'lower')
        
        txt = r'%s' % K.califaID
        kw_text = dict(pos_x = 0.98, pos_y = 0.98, fs = 20, va = 'top', ha = 'right', c = 'white')
        plot_text_ax(ax, txt, **kw_text)
           
        txt = r'%s' % mtype_str[np.where(mtype_num == mtype)[0][0]]
        kw_text = dict(pos_x = 0.02, pos_y = 0.02, fs = 20, va = 'bottom', ha = 'left', c = 'white')
        plot_text_ax(ax, txt, **kw_text)

        txt = r'b/a %.2f [pyc: %.2f]' % (ba, ba_pyc)
        kw_text = dict(pos_x = 0.99, pos_y = 0.01, fs = 15, va = 'bottom', ha = 'right', c = 'white')
        plot_text_ax(ax, txt, **kw_text)
        
        DrawHLRCircleInSDSSImage(ax, K.HLR_pix, K.pa, K.ba)
        loc = [ [0, 2], [1, 2], [2, 0], [2, 2], [3, 0], [3, 2] ]
        N = len(loc)
        for i in xrange(N):
            pos = loc[i][:]
            ax = plt.subplot2grid(grid_shape, loc = pos)
            ax.set_axis_on()
            if i < N - 1:
                y = np.ma.log10(phys__yx[i])
            else:
                y = phys__yx[i]
            im = plot_ma_map(ax = ax, map = y, p = [1, 99], K = K, **dict(cmap = mpl.cm.hot))
            cb = f.colorbar(ax = ax, mappable = im)
            ax.set_title(labels__r[i], stretch = 'condensed', y = 1.01, fontsize = 10)
            pos[-1] += 1
            ax = plt.subplot2grid(grid_shape, loc = pos)
            ax.set_axis_on()
            ax.plot(H.RbinCenter__r, phys__r[i], 'o-')
            ax.set_xlim(minR, H.Rbin__r[-1])
        f.subplots_adjust(hspace = 0.2, wspace = 0.3)
        pdf.savefig(f)
        plt.close(f)

        ##########################
        ######### PAGE 2 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        #f, axArr = plt.subplots(NRows, NCols)
        page_size_inches = (NCols * 3, NRows * 1.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        x = H.RbinCenter__r
        xm, ym_BR_DGR_ols__r = C.ma_mask_xyz(x, y = BR_DGR_ols__r, mask = ~maskRadiusOk__r)
        ax.plot(xm, ym_BR_DGR_ols__r, '.-', c = 'b')
        m_aux = ~maskRadiusOk__r | BR_DGR_up_ols__r.mask | BR_DGR_up_ols__r.mask
        xm, ym_BR_DGR_up_ols__r = C.ma_mask_xyz(x, y = BR_DGR_up_ols__r, mask = m_aux)
        xm, ym_BR_DGR_down_ols__r = C.ma_mask_xyz(x, y = BR_DGR_down_ols__r, mask = m_aux)
        ax.fill_between(xm, ym_BR_DGR_up_ols__r, ym_BR_DGR_down_ols__r, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        xm, ym_SK_DGR__r = C.ma_mask_xyz(x, y = SK_DGR__r, mask = ~maskRadiusOk__r) 
        ax.plot(xm, ym_SK_DGR__r, '.-', c = 'g')
        xm, ym_SK_DGR_Ha_r = C.ma_mask_xyz(x, y = SK_DGR_Ha__r, mask = ~maskRadiusOk__r) 
        ax.plot(xm, ym_SK_DGR_Ha_r, '.-', c = 'black')
        ax.axhline(y = RR_DGR, c = 'y')
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(0, .025)
        ax.set_title('Radial bins')
        ax.set_ylabel(r'$\delta_{DGR}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
          
        if args.debug is not True:
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.3
            ax = plt.subplot2grid(grid_shape, loc = (0, 1))
            ax.set_axis_on()
            x = K.zoneDistance_HLR
            xm, ym_BR_DGR_ols__z = C.ma_mask_xyz(x, y = BR_DGR_ols__z, mask = ~maskRadiusOk__z)
            rs_BR_DGR_ols = runstats(xm.compressed(), ym_BR_DGR_ols__z.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_DGR_up_ols__z, mask = ~maskRadiusOk__z) 
            rs_BR_DGR_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_DGR_down_ols__z, mask = yup.mask) 
            rs_BR_DGR_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = SK_DGR__z, mask = ~maskRadiusOk__z)
            rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = SK_DGR_Ha__z, mask = ~maskRadiusOk__z)
            rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            if args.scatter is True:
                ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = 'b', **sc_kwargs)
                ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
                ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'black', **sc_kwargs)
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b', label = r'BR from $\tau_V^\star$')
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g', label = r'SK from synt.')
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black', label = r'SK from H$\alpha$')
            ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
            ax.axhline(y = RR_DGR, label = 'RR (const.)', c = 'y')
            plt.setp(ax.get_xticklabels(), visible = False)
            ax.set_ylim(0, .025)
            ax.set_xlim(minR, H.Rbin__r[-1])
            ax.tick_params(axis = 'y', which = 'major', pad = 10)
            ax.set_title('Zones')
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
         
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_ols__r), mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'b')
        xm, yup = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_up_ols__r), mask = ~maskRadiusOk__r) 
        xm, ydown = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_down_ols__r), mask = yup.mask) 
        ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_Ha_ols__r), mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'r')
        xm, yup = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_Ha_up_ols__r), mask = ~maskRadiusOk__r) 
        xm, ydown = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(BR_SigmaGas_Ha_down_ols__r), mask = yup.mask) 
        ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(SK_SigmaGas__r), mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'g')
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(SK_SigmaGas_Ha__r), mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'k')
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = np.ma.log10(RR_SigmaGas__r), mask = ~maskRadiusOk__r) 
        ax.plot(xm, ym, '.-', c = 'y')
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(0.4, 2)
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
         
        if args.debug is not True:
            ax = plt.subplot2grid(grid_shape, loc = (1, 1))
            ax.set_axis_on()
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_ols__z), mask = ~maskRadiusOk__z) 
            rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_up_ols__z), mask = ~maskRadiusOk__z) 
            rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_down_ols__z), mask = yup.mask) 
            rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_Ha_ols__z), mask = ~maskRadiusOk__z)
            rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_Ha_up_ols__z), mask = ~maskRadiusOk__z) 
            rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(BR_SigmaGas_Ha_down_ols__z), mask = yup.mask) 
            rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(SK_SigmaGas__z), mask = ~maskRadiusOk__z)
            rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(SK_SigmaGas_Ha__z), mask = ~maskRadiusOk__z)
            rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = np.ma.log10(RR_SigmaGas__z), mask = ~maskRadiusOk__z) 
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
            ax.set_xlim(minR, H.Rbin__r[-1])
            ax.set_ylim(0.4, 2)
            ax.legend(bbox_to_anchor = (2.6, 2), fontsize = 10, frameon = False, ncol = 2)#, loc = 'upper right')
            ax.tick_params(axis = 'y', which = 'major', pad = 15)
            plt.setp(ax.get_xticklabels(), visible = False)
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
        
        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_ols__r, mask = ~maskRadiusOk__r) 
        ax.plot(xm, ym, '.-', c = 'b')
        xm, yup = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_up_ols__r, mask = ~maskRadiusOk__r) 
        xm, ydown = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_down_ols__r, mask = yup.mask) 
        ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_Ha_ols__r, mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'r')
        xm, yup = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_Ha_up_ols__r, mask = ~maskRadiusOk__r) 
        xm, ydown = C.ma_mask_xyz(x = H.RbinCenter__r, y = BR_f_gas_Ha_down_ols__r, mask = yup.mask) 
        ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = SK_f_gas__r, mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'g')
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = SK_f_gas_Ha__r, mask = ~maskRadiusOk__r)
        ax.plot(xm, ym, '.-', c = 'k')
        xm, ym = C.ma_mask_xyz(x = H.RbinCenter__r, y = RR_f_gas__r, mask = ~maskRadiusOk__r) 
        ax.plot(xm, ym, '.-', c = 'y')
        ax.set_xlim(minR, H.RbinFin)
        ax.set_ylim(0, .3)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylabel(r'f${}_{gas}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_yticklabels(), visible = False)
                     
        if args.debug is not True:
            ax = plt.subplot2grid(grid_shape, loc = (2, 1))
            ax.set_axis_on()
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_ols__z, mask = ~maskRadiusOk__z) 
            rs_BR_f_gas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_up_ols__z, mask = ~maskRadiusOk__z) 
            rs_BR_f_gas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_down_ols__z, mask = yup.mask) 
            rs_BR_f_gas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_Ha_ols__z, mask = ~maskRadiusOk__z)
            rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_Ha_up_ols__z, mask = ~maskRadiusOk__z) 
            rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = BR_f_gas_Ha_down_ols__z, mask = yup.mask) 
            rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = SK_f_gas__z, mask = ~maskRadiusOk__z)
            rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = SK_f_gas_Ha__z, mask = ~maskRadiusOk__z)
            rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = K.zoneDistance_HLR, y = RR_f_gas__z, mask = ~maskRadiusOk__z) 
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
            ax.set_xlim(minR, H.Rbin__r[-1])
            ax.set_ylim(0, .3)
            ax.tick_params(axis = 'y', which = 'major', pad = 15)
            ax.set_xlabel(r'R [HLR]')
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
     
        ax = plt.subplot2grid(grid_shape, loc = (0, NCols - 1), rowspan = NRows)
        ax.set_axis_off()
        txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%' % (H.N_gals, (H.tSF__T[iT] / 1e6), H.xOkMin * 100.)
        kw_text = dict(pos_x = -0.2, pos_y = 0.65, fs = 11, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        kw_text = dict(pos_x = -0.2, pos_y = 0.55, fs = 11, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$\delta_{DGR} = (%.2f\times 10^{-3} - %.2f\times 10^{-2}) (\frac{O/H}{(O/H)_\odot})$' % (DGR_conv_lim_inf / 1e-3, DGR_conv_lim_sup / 1e-2)
        kw_text = dict(pos_x = -0.2, pos_y = 0.40, fs = 12, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        txt = r'$f_{gas}\ =\ [1 + (\frac{\sigma_d}{m_d})(\delta_{DGR})(\frac{\mu_\star}{\tau_V})]^{-1}$'
        kw_text = dict(pos_x = -0.2, pos_y = 0.25, fs = 15, va = 'bottom', ha = 'left', c = 'k')
        plot_text_ax(ax, txt, **kw_text)
        f.subplots_adjust(hspace = 0.2, wspace = 0.3)
        pdf.savefig(f)
        plt.close(f)

        ##########################
        ##########################
        ##########################
        
        xaxis__r = phys__r[:]
        xaxis__r.append(at_flux__Tr[iT])
        labels = labels__r[:]
        labels.append(r'$\langle \log\ t \rangle_L\ (R)$ [yr]')
        
        xaxis__z = phys__z[:]
        xaxis__z.append(at_flux__Tz[iT])
        lbls__z = labels__z[:]
        lbls__z.append(r'$\langle \log\ t \rangle_L$ [yr]')
        
        ##########################
        ######### PAGE 3 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\delta_{DGR}\ \times$ properties by radius')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        for i, x in enumerate(xaxis__r):
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 10
            sc_kwargs['alpha'] = 0.9
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_DGR_ols__r, mask = ~maskRadiusOk__r)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_DGR_up_ols__r, mask = ~maskRadiusOk__r)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_DGR_down_ols__r, mask = ~maskRadiusOk__r)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            #ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR__r, mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR_Ha__r, mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0, 0.025)
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
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\delta_{DGR}\ \times$ properties by zones')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.4
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_DGR_ols__z, mask = ~maskRadiusOk__z)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_DGR_up_ols__z, mask = ~maskRadiusOk__z)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_DGR_down_ols__z, mask = ~maskRadiusOk__z)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            #ax.fill_between(xm, yup, ydown, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR__z, mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR_Ha__z, mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0, 0.025)
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
        ######### PAGE 5 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\delta_{DGR}\ \times$ properties by zones with runstats')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.3
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_DGR_ols__z, mask = ~maskRadiusOk__z)
            rs_BR_DGR_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
            ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = 'b', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR__z, mask = ~maskRadiusOk__z)
            rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
            ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_DGR_Ha__z, mask = ~maskRadiusOk__z)
            rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'k')
            ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'k', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0, 0.025)
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
        ######### PAGE 6 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\log\ \Sigma_{gas}\ \times$ properties by radius')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        for i, x in enumerate(xaxis__r):
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 10
            sc_kwargs['alpha'] = 0.9
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_ols__r), mask = ~maskRadiusOk__r)
            xm, yup = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_up_ols__r), mask = ~maskRadiusOk__r)
            xm, ydown = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_down_ols__r), mask = ~maskRadiusOk__r)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_ols__r), mask = ~maskRadiusOk__r)
            xm, yup = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_up_ols__r), mask = ~maskRadiusOk__r)
            xm, ydown = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_down_ols__r), mask = ~maskRadiusOk__r)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'r', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas__r), mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas_Ha__r), mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.log10(RR_SigmaGas__r), mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'y', **sc_kwargs)

            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
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
        ######### PAGE 7 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\log\ \Sigma_{gas}\ \times$ properties by zones')
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.4
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_ols__z), mask = ~maskRadiusOk__z)
            xm, yup = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_up_ols__z), mask = ~maskRadiusOk__z)
            xm, ydown = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_down_ols__z), mask = ~maskRadiusOk__z)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_Ha_ols__z), mask = ~maskRadiusOk__z)
            xm, yup = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_Ha_up_ols__z), mask = ~maskRadiusOk__z)
            xm, ydown = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_Ha_down_ols__z), mask = ~maskRadiusOk__z)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'r', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(SK_SigmaGas__z), mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(SK_SigmaGas_Ha__z), mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(RR_SigmaGas__z), mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'y', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0.4, 2)
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
        ######### PAGE 8 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$\log\ \Sigma_{gas}\ \times$ properties by zones with runstats')
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.3
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_ols__z), mask = ~maskRadiusOk__z)
            rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b')
            ax.scatter(rs_BR_SigmaGas_ols.x, rs_BR_SigmaGas_ols.y, c = 'b', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(BR_SigmaGas_Ha_ols__z), mask = ~maskRadiusOk__z)
            rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r')
            ax.scatter(rs_BR_SigmaGas_Ha_ols.x, rs_BR_SigmaGas_Ha_ols.y, c = 'r', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(SK_SigmaGas__z), mask = ~maskRadiusOk__z)
            rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
            ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(SK_SigmaGas_Ha__z), mask = ~maskRadiusOk__z)
            rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'k')
            ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, c = 'k', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = np.ma.log10(RR_SigmaGas__z), mask = ~maskRadiusOk__z)
            rs_RR_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y')
            ax.scatter(rs_RR_SigmaGas.x, rs_RR_SigmaGas.y, c = 'y', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0.4, 2)
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

        #########################
        ######## PAGE 9 #########
        #########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$f_{gas}\ \times$ properties by radius')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        for i, x in enumerate(xaxis__r):
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 10
            sc_kwargs['alpha'] = 0.9
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_ols__r, mask = ~maskRadiusOk__r)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_f_gas_up_ols__r, mask = ~maskRadiusOk__r)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_f_gas_down_ols__r, mask = ~maskRadiusOk__r)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_ols__r, mask = ~maskRadiusOk__r)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_up_ols__r, mask = ~maskRadiusOk__r)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_down_ols__r, mask = ~maskRadiusOk__r)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'r', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas__r, mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__r, mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = RR_f_gas__r, mask = ~maskRadiusOk__r)
            ax.scatter(xm, ym, c = 'y', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
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

        ##########################
        ######### PAGE 10 ########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$f_{gas}\ \times$ properties by zones')
        #f, axArr = plt.subplots(NRows, NCols)
        #f.set_size_inches((NCols * 3, NRows * 1.5))
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.4
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_ols__z, mask = ~maskRadiusOk__z)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_f_gas_up_ols__z, mask = ~maskRadiusOk__z)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_f_gas_down_ols__z, mask = ~maskRadiusOk__z)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'b', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_ols__z, mask = ~maskRadiusOk__z)
            xm, yup = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_up_ols__z, mask = ~maskRadiusOk__z)
            xm, ydown = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_down_ols__z, mask = ~maskRadiusOk__z)
            ax.errorbar(xm, ym, yerr = [ ym - ydown, yup - ym ], c = 'r', linestyle='None', fmt='.')
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas__z, mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__z, mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'black', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = RR_f_gas__z, mask = ~maskRadiusOk__z)
            ax.scatter(xm, ym, c = 'y', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0, 0.4)
            if col == 0:
                ax.set_ylabel(r'$f_{gas}$')
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)       
 
        ##########################
        ######### PAGE 11 ########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        f.suptitle(r'$f_{gas}\ \times$ properties by zones with runstats')
        page_size_inches = (NCols * 3, NRows * 2.5)
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        row, col = 0, 0
        N = len(xaxis__z)
        for i, x in enumerate(xaxis__z):
            if i < N - 3:
                x = np.ma.log10(x)
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
            sc_kwargs['s'] = 5
            sc_kwargs['alpha'] = 0.3
            ax = plt.subplot2grid(grid_shape, loc = (row, col))
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_ols__z, mask = ~maskRadiusOk__z)
            rs_BR_f_gas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'b')
            ax.scatter(rs_BR_f_gas_ols.x, rs_BR_f_gas_ols.y, c = 'b', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_ols__z, mask = ~maskRadiusOk__z)
            rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'r')
            ax.scatter(rs_BR_f_gas_Ha_ols.x, rs_BR_f_gas_Ha_ols.y, c = 'r', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas__z, mask = ~maskRadiusOk__z)
            rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
            ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, c = 'g', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__z, mask = ~maskRadiusOk__z)
            rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'k')
            ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, c = 'k', **sc_kwargs)
            xm, ym = C.ma_mask_xyz(x = x, y = RR_f_gas__z, mask = ~maskRadiusOk__z)
            rs_RR_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'y')
            ax.scatter(rs_RR_f_gas.x, rs_RR_f_gas.y, c = 'y', **sc_kwargs)
            ax.xaxis.set_major_locator(MaxNLocator(6))
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
            ax.set_xlabel(labels[i])
            ax.set_ylim(0, 0.4)
            if col == 0:
                ax.set_ylabel(r'$f_{gas}$')
            if col == (NCols - 1):
                row += 1
                col = 0
            else:
                col += 1
        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)    

