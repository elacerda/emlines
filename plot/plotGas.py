#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
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
        'slice_gals' : None, 
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
    SK_f_gas__g = 1. / (1. + 1./SK_GSR__g)
    SK_f_gas_Ha__g = 1. / (1. + 1./SK_GSR_Ha__g)
    SK_f_gas__rg = 1. / (1. + 1./SK_GSR__rg)
    SK_f_gas_Ha__rg = 1. / (1. + 1./SK_GSR_Ha__rg)
    SK_integrated_f_gas = 1. / (1. + 1./SK_integrated_GSR)
    SK_integrated_f_gas_Ha = 1. / (1. + 1./SK_integrated_GSR_Ha)
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
    RR_f_gas__g = 1. / (1. + 1./RR_GSR__g)
    RR_f_gas_Ha__g = 1. / (1. + 1./RR_GSR_Ha__g)
    RR_f_gas__rg = 1. / (1. + 1./RR_GSR__rg)
    RR_f_gas_Ha__rg = 1. / (1. + 1./RR_GSR_Ha__rg)
    RR_integrated_f_gas = 1. / (1. + 1./RR_integrated_GSR)
    RR_integrated_f_gas_Ha = 1. / (1. + 1./RR_integrated_GSR_Ha)
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
    BR_f_gas_up_ols__g = 1. / (1. + 1./BR_GSR_up_ols__g)
    BR_f_gas_up_ols__rg = 1. / (1. + 1./BR_GSR_up_ols__rg)
    BR_f_gas_Ha_up_ols__g = 1. / (1. + 1./BR_GSR_Ha_up_ols__g)
    BR_f_gas_Ha_up_ols__rg = 1. / (1. + 1./BR_GSR_Ha_up_ols__rg)
    BR_f_gas_down_ols__g = 1. / (1. + 1./BR_GSR_down_ols__g)
    BR_f_gas_down_ols__rg = 1. / (1. + 1./BR_GSR_down_ols__rg)
    BR_f_gas_Ha_down_ols__g = 1. / (1. + 1./BR_GSR_Ha_down_ols__g)
    BR_f_gas_Ha_down_ols__rg = 1. / (1. + 1./BR_GSR_Ha_down_ols__rg)
    BR_f_gas_ols__g = 1. / (1. + 1./BR_GSR_ols__g)
    BR_f_gas_ols__rg = 1. / (1. + 1./BR_GSR_ols__rg)
    BR_f_gas_Ha_ols__g = 1. / (1. + 1./BR_GSR_Ha_ols__g)
    BR_f_gas_Ha_ols__rg = 1. / (1. + 1./BR_GSR_Ha_ols__rg)
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
    BR_f_gas_up_cubic__g = 1. / (1. + 1./BR_GSR_up_cubic__g)
    BR_f_gas_up_cubic__rg = 1. / (1. + 1./BR_GSR_up_cubic__rg)
    BR_f_gas_Ha_up_cubic__g = 1. / (1. + 1./BR_GSR_Ha_up_cubic__g)
    BR_f_gas_Ha_up_cubic__rg = 1. / (1. + 1./BR_GSR_Ha_up_cubic__rg)
    BR_f_gas_down_cubic__g = 1. / (1. + 1./BR_GSR_down_cubic__g)
    BR_f_gas_down_cubic__rg = 1. / (1. + 1./BR_GSR_down_cubic__rg)
    BR_f_gas_Ha_down_cubic__g = 1. / (1. + 1./BR_GSR_Ha_down_cubic__g)
    BR_f_gas_Ha_down_cubic__rg = 1. / (1. + 1./BR_GSR_Ha_down_cubic__rg)
    BR_f_gas_cubic__g = 1. / (1. + 1./BR_GSR_cubic__g)
    BR_f_gas_cubic__rg = 1. / (1. + 1./BR_GSR_cubic__rg)
    BR_f_gas_Ha_cubic__g = 1. / (1. + 1./BR_GSR_Ha_cubic__g)
    BR_f_gas_Ha_cubic__rg = 1. / (1. + 1./BR_GSR_Ha_cubic__rg)
    ######################
    
    ######################
    # Running Stats
    
    if args.output is None:
        if args.maskradius is not None:
            output = 'gas_maskRadius%.1f.pdf' % args.maskradius
        else:
            output = 'gas.pdf'
    else:
        output = args.output
          
    default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)

    with PdfPages(output) as pdf:
        ##########################
        ######### PAGE 1 #########
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
        tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
        rs_kwargs = default_rs_kwargs.copy()
        sc_kwargs = default_sc_kwargs.copy()
        ax.axhline(y = RR_DGR, c = 'y')
        xm, ym = C.ma_mask_xyz(x = H.Rtoplot(), y = BR_DGR_ols__rg, mask = tmp_mask)
        rs_BR_DGR_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
        m_aux = tmp_mask| BR_DGR_up_ols__rg.mask | BR_DGR_up_ols__rg.mask
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
        ax.set_xlim(minR,  H.RbinFin)
        #ax.set_ylim(0, .5)
        ax.set_title('Radial bins')
        ax.set_ylabel(r'$\delta_{DGR}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
          
        if args.debug is not True:
            ax = plt.subplot2grid(grid_shape, loc = (0, 1))
            ax.set_axis_on()
            tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_ols__g, mask = tmp_mask) 
            rs_BR_DGR_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            m_aux = tmp_mask | BR_DGR_up_ols__g.mask | BR_DGR_up_ols__g.mask 
            xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_up_ols__g, mask = m_aux)
            rs_BR_DGR_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_DGR_down_ols__g, mask = m_aux) 
            rs_BR_DGR_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_DGR__g, mask = tmp_mask) 
            rs_SK_DGR = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_DGR_Ha__g, mask = tmp_mask) 
            rs_SK_DGR_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            if args.scatter is True:
                ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = 'b', **sc_kwargs)
                ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = 'g', **sc_kwargs)
                ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = 'g', **sc_kwargs)
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
            ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
            ax.axhline(y = RR_DGR, label = 'RR (const.)', c = 'y')
            plt.setp(ax.get_xticklabels(), visible = False)
            #ax.set_ylim(0, .5)
            ax.set_xlim(minR,  H.Rbin__r[-1])
            ax.tick_params(axis = 'y', which = 'major', pad = 10)
            ax.set_title('Zones')
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
         
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
        ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b')
        ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r')
        ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
        ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black')
        ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y')
        ax.fill_between(rs_BR_SigmaGas_up_ols.xS, rs_BR_SigmaGas_up_ols.yS, rs_BR_SigmaGas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)
        ax.fill_between(rs_BR_SigmaGas_Ha_up_ols.xS, rs_BR_SigmaGas_Ha_up_ols.yS, rs_BR_SigmaGas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
        ax.set_xlim(minR,  H.RbinFin)
        ax.set_ylim(0.4, 2)
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
         
        if args.debug is not True:
            ax = plt.subplot2grid(grid_shape, loc = (1, 1))
            ax.set_axis_on()
            tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_ols__g), mask = tmp_mask) 
            rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_up_ols__g), mask = tmp_mask) 
            rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_down_ols__g), mask = yup.mask) 
            rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_ols__g), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_up_ols__g), mask = tmp_mask) 
            rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(BR_SigmaGas_Ha_down_ols__g), mask = yup.mask) 
            rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(SK_SigmaGas__g), mask = tmp_mask)
            rs_SK_SigmaGas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(SK_SigmaGas_Ha__g), mask = tmp_mask)
            rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = np.ma.log10(RR_SigmaGas__g), mask = tmp_mask) 
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
            ax.set_xlim(minR,  H.Rbin__r[-1])
            ax.set_ylim(0.4, 2)
            ax.legend(bbox_to_anchor = (2.6, 2), fontsize = 10, frameon = False, ncol = 2)#, loc = 'upper right')
            ax.tick_params(axis = 'y', which = 'major', pad = 15)
            plt.setp(ax.get_xticklabels(), visible = False)
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
         
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
        ax.set_xlim(minR,  H.RbinFin)
        ax.set_ylim(0, .4)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylabel(r'f${}_{gas}$')
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.grid()
        plt.setp(ax.get_yticklabels(), visible = False)
                     
        if args.debug is not True:
            ax = plt.subplot2grid(grid_shape, loc = (2, 1))
            ax.set_axis_on()
            tmp_mask = ~(maskRadiusOk__g & gals_slice__g)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_ols__g, mask = ~maskRadiusOk__g) 
            rs_BR_f_gas_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_up_ols__g, mask = tmp_mask) 
            rs_BR_f_gas_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_down_ols__g, mask = yup.mask) 
            rs_BR_f_gas_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_ols__g, mask = tmp_mask)
            rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, yup = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_up_ols__g, mask = tmp_mask) 
            rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), yup.compressed(), nBox = 20, **rs_kwargs)
            xm, ydown = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = BR_f_gas_Ha_down_ols__g, mask = yup.mask) 
            rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ydown.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_f_gas__g, mask = tmp_mask)
            rs_SK_f_gas = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = SK_f_gas_Ha__g, mask = tmp_mask)
            rs_SK_f_gas_Ha = runstats(xm.compressed(), ym.compressed(), nBox = 20, **rs_kwargs)
            xm, ym = C.ma_mask_xyz(x = H.zone_dist_HLR__g, y = RR_f_gas__g, mask = tmp_mask) 
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
            ax.set_xlim(minR,  H.Rbin__r[-1])
            ax.set_ylim(0, .4)
            ax.tick_params(axis = 'y', which = 'major', pad = 15)
            ax.set_xlabel(r'R [HLR]')
            ax.yaxis.set_major_locator(MaxNLocator(6))
            ax.grid()
         
        ax = plt.subplot2grid(grid_shape, loc = (0, NCols - 1), rowspan = NRows)
        ax.set_axis_off()
        txt = r'NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%' % (N_gals, (H.tSF__T[iT] / 1e6), H.xOkMin * 100.)
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
        ######### PAGE 2 #########
        ##########################
        NRows = 3
        NCols = 3
        f = plt.figure()
        page_size_inches = (NCols * 3, NRows * 2.5)
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
            ax.set_xlim(minR,  H.RbinFin)
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
        page_size_inches = (NCols * 3, NRows * 2.5)
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
            ax.set_ylim(0,15e-3)
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
        page_size_inches = (NCols * 3, NRows * 2.5)
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
        page_size_inches = (NCols * 3, NRows * 2.5)
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
