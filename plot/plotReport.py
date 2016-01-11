#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014 
#                 - 10/Sept/2015
#
from matplotlib.ticker import MultipleLocator
import sys
import CALIFAUtils as C
from scipy import stats as st
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
#from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.plots import plot_linreg_params
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.scripts import OLS_bisector
#from CALIFAUtils.plots import plot_zbins
#from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
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
        'output' : '',
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
                        default = default['scatter'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
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
    
    minR = 0
    fnamesuffix = '.pdf'
    
    if args.maskradius is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
        maxR = H.Rbin__r[-1]
        maxR = 3
        maskRadiusOk__g = (H.zone_dist_HLR__g >= minR) & (H.zone_dist_HLR__g <= maxR) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * ((H.RbinCenter__r >= minR) & (H.RbinCenter__r <= maxR))).T
        fnamesuffix = '_maskradius.pdf'
        
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1]
        fnamesuffix = '_%s%s' % ('.'.join(gals_txt.split('.')[:-1]), fnamesuffix)
        #integrated
        l_gals, _ = C.sort_gals(args.slice_gals)
        l_gals = sorted(l_gals)
        gals_slice__integr = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
        for g in l_gals:
            i = H.califaIDs_all.tolist().index(g)
            gals_slice__integr[i] = True

    txt_suptitle = r'$\Longrightarrow$ %s  NGals:%d  $x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f ' % (gals_txt, N_gals, H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)    
    #tSF__T              = H.tSF__T[0:20]
    tSF__T = H.tSF__T
    RRange = H.RRange
    RColor = H.RColor
    SFR__Tg = H.get_data_h5('SFR__Tg')
    SFR_Ha__g = H.get_data_h5('SFR_Ha__g')
    aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
    aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
    SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
    SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
    dist_zone__g = H.get_data_h5('dist_zone__g')
     
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    f.suptitle(txt_suptitle, fontsize = 11)
    filename = 'SFR_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFR__Tg[iT])
            y = np.ma.log10(SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-6, 0]
            yran = [-6, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
            
    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            x = np.ma.log10(SFR__Tg[i])
            y = np.ma.log10(SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
            
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    xlabel = r'$\log\ t_\star$ [yr]'  
    x = np.log10(tSF__T)
    plot_linreg_params(a, x, xlabel,
                       'a', 'SFR_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age%s' % fnamesuffix)
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical') 
    filename = 'SFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)  
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            x = np.ma.log10(SFRSD__Tg[i] * 1e6)
            y = np.ma.log10(SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'SFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'SFRSD_Rs_age%s' % fnamesuffix, 1., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')   
    filename = 'aSFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(aSFRSD__Trg[iT] * 1e6)
            y = np.ma.log10(aSFRSD_Ha__rg * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            x = np.ma.log10(aSFRSD__Trg[i] * 1e6)
            y = np.ma.log10(aSFRSD_Ha__rg * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())
             
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'aSFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'aSFRSD_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
    ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}(R)}{\Sigma_{SFR}^{neb}(@1HLR)}$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    fnamesuftmp = '_norm%s' % fnamesuffix   
    filename = 'aSFRSD_linregress_report%s' % fnamesuftmp
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j]
            age = tSF__T[iT]
            aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]
            aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / H.aSFRSD_Ha_oneHLR__g
            xran = [-1.5, 1.5]
            yran = [-1.5, 1.5]
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # binsx = np.linspace(-4.5, 1., 51)
            # binsy = np.linspace(min(ym),max(ym), 51)
            # density_contour(xm, ym, binsx, binsy, ax=ax)
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            x = np.ma.log10(aSFRSD_norm__rg)
            y = np.ma.log10(aSFRSD_Ha_norm__rg)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.6)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            aSFRSD_norm__rg = H.aSFRSD__Trg[i] / H.aSFRSD_oneHLR__Tg[i]
            aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / H.aSFRSD_Ha_oneHLR__g
            x = np.ma.log10(aSFRSD_norm__rg)
            y = np.ma.log10(aSFRSD_Ha_norm__rg)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'aSFRSD_linregress_slope_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep2_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'aSFRSD_Rs_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuftmp, 1., 16) 
    ############# Integrated #############
    ############# Integrated #############
    ############# Integrated #############
    integrated_SFR__Tg = H.integrated_SFR__Tg
    integrated_SFR_Ha__g = H.integrated_SFR_Ha__g
    integrated_SFRSD__Tg = H.integrated_SFRSD__Tg
    integrated_SFRSD_Ha__g = H.integrated_SFRSD_Ha__g
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star^{int}(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}^{int}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFR_linregress_report.png'
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFR__Tg[iT])
            y = np.ma.log10(integrated_SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-5, 2]
            yran = [-5, 2]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            x = np.ma.log10(integrated_SFR__Tg[i])
            y = np.ma.log10(integrated_SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'integrated_SFR_linregress_slope_age.png', 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFR_linregress_intercep_age.png', 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFR_linregress_intercep2_age.png', 0., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rs', 'SFR_Rs_age.png', 1., 16) 
     
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age.png')
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age.png', 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'integrated_SFR_Rs_age.png', 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(int, t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(int)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFRSD_linregress_report.png'
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)  
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-5, 0]
            yran = [-5, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1

    if iT < (H.N_T - 1):
        for i in range(iT, H.N_T):
            x = np.ma.log10(integrated_SFRSD__Tg[i] * 1e6)
            y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            a[i], b[i], sigma_a, sigma_b = OLS_bisector(xm, ym)
            b2[i] = (ym - xm).mean()
            Rs[i], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[i], _ = st.pearsonr(xm.compressed(), ym.compressed())

    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'integrated_SFRSD_linregress_slope_age.png', 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFRSD_linregress_intercep_age.png', 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFRSD_linregress_intercep2_age.png', 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'integrated_SFRSD_Rs_age.png', 1., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age.png', 1., 16) 
  
