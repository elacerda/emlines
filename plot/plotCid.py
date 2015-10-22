#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import sys
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
from os.path import basename
from matplotlib import pyplot as plt
from CALIFAUtils.objects import runstats 
from matplotlib.pyplot import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plotOLSbisectorAxis, plot_text_ax

#mask_radius = False
RNuc = 0.7
maskradius = RNuc
#maskradius = None
if maskradius is None: RNuc = 0 
#zlim = [RNuc, 3]
A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

def iTGen(tSF__T, iT_values = [ 11, 17 ]):
    for iT in iT_values:
        yield iT, tSF__T[iT]

mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)

H = C.H5SFRData(h5file)
tSF__T = H.tSF__T
#iT_values = [ 2, 3 ]
iT_values = [ 2 ]
zlim = [RNuc, H.Rbin__r[-1]]

if maskradius is None:
    namepref = ''
    maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
    maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
else:
    namepref = '_R%.1f' % RNuc
    minR = maskradius
    maskRadiusOk__g = (H.zone_dist_HLR__g >= maskradius) & (H.zone_dist_HLR__g <= H.Rbin__r[-1]) 
    maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= maskradius)).T

default_kwargs_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 50)
default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.6, label = '')
default_kwargs_zbins = dict(
    debug = True,
    zlabel = r'R (HLR)',
    zmask = None,
    #xlim = [-1.5, 0.5],
    #ylim = [-3.5, 1],
    zlim = zlim,
    x_major_locator = 1.,
    x_minor_locator = 0.2,
    y_major_locator = 1.,
    y_minor_locator = 0.2,
    ols = True,
    running_stats = True,
    rs_gaussian_smooth = True,
    kwargs_rs = default_kwargs_rs,
    kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)', alpha = 0.9),
    kwargs_figure = dict(figsize = (10, 8), dpi = 100),
    kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.45, label = ''),
    kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 6, rms = True, text = True),
    kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS'),
    write_N = True,
)

kwargs_ols_plot = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
kwargs_ols = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 6, rms = True, text = True, kwargs_plot = kwargs_ols_plot)


for iT, tSF in iTGen(H.tSF__T, iT_values):
    ba_max = 0
    #mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.less(H.ba_GAL__g, ba_max))
    #mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, ba_max))
    mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.zeros_like(H.ba_GAL__g, dtype = np.bool))
    #mask_GAL__g = np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool)
    
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.O_O3N2_M13__g.mask)
    mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    #mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), ba_max))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    #mask__g = ~maskRadiusOk__g
    
    mask__rg = np.bitwise_or(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Trg[iT]).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.aSFRSD_Ha__rg * 1e6).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.tau_V_neb__rg).mask)
    mask__rg = np.bitwise_or(mask__rg, H.O_O3N2_M13__rg.mask)
    mask__rg = np.bitwise_or(mask__rg, np.less(H.EW_Ha__rg, 3.))
    #mask__rg = np.bitwise_or(mask__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    mask__rg = np.bitwise_or(mask__rg, ~maskRadiusOk__rg)
    #mask__rg = ~maskRadiusOk__rg
            
    with PdfPages('SK_%.2fMyr%s_%s.pdf' % ((tSF/1e6), namepref, basename(h5file).replace('SFR_', '').replace('.h5', ''))) as pdf:
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)

        if maskradius is None:
            suptitle = r'NGals:%d  Nzones:%d  NRbins:%d  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (~mask__g).sum(), (~mask__rg).sum(), (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        else:
            suptitle = r'NGals:%d  Nzones:%d  NRbins:%d  R > %.1fHLR  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (~mask__g).sum(), (~mask__rg).sum(), RNuc, (tSF/1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    
        f.suptitle(suptitle, fontsize = 9)
        #f, axArr = plt.subplots(1, 3)
        #f.set_dpi(100)
        #f.set_size_inches(15, 8)
         
        tau_V_norm_GAL__g = (H.tau_V__Trg[iT][9:11, :]).mean(axis = 0)
        tau_V_norm__rg = H.tau_V__Trg[iT] / tau_V_norm_GAL__g
        SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][9:11, :]).mean(axis = 0)
        aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g

        x1 = np.ma.masked_array(np.ma.log10(H.tau_V__Tg[iT]), mask = mask__g)
        y1 = np.ma.masked_array(np.ma.log10(H.SFRSD__Tg[iT] * 1e6), mask = mask__g)
        x2 = np.ma.masked_array(np.ma.log10(H.tau_V__Trg[iT]), mask = mask__rg)
        y2 = np.ma.masked_array(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6), mask = mask__rg) 
        x3 = np.ma.masked_array(np.ma.log10(tau_V_norm__rg), mask = mask__rg)
        y3 = np.ma.masked_array(np.ma.log10(aSFRSD_norm__rg), mask = mask__rg) 
        x4 = np.ma.masked_array(np.ma.log10(tau_V_norm_GAL__g), mask = mask_GAL__g) 
        y4 = np.ma.masked_array(np.ma.log10(SFRSD_norm_GAL__g * 1e6), mask = mask_GAL__g)
        x5 = np.ma.masked_array(np.ma.log10(H.integrated_tau_V__g), mask = mask_GAL__g) 
        y5 = np.ma.masked_array(np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6), mask = mask_GAL__g)
        
        x1label = r'$\log\ \tau_V^{\star}$'
        y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{\star}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
        x4label = r'$\log\ \tau_V^\star(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^\star(@1HLR)$'
        x5label = r'$\log\ \tau_V^\star (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^\star (int)$'
    
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label, z = H.zone_dist_HLR__g, **kwargs_zbins)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(),  z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)    

        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label, **kwargs_zbins)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        ax.grid()
        
        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
        
        ############## Neb ##############
        ############## Neb ##############
        ############## Neb ##############
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
    
        f.suptitle(suptitle, fontsize = 9)

        tau_V_neb_norm_GAL__g = (H.tau_V_neb__rg[9:11, :]).mean(axis = 0)
        tau_V_neb_norm__rg = H.tau_V_neb__rg / tau_V_neb_norm_GAL__g
        SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[9:11, :]).mean(axis = 0)
        aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g
     
        x1 = np.ma.masked_array(np.ma.log10(H.tau_V_neb__g), mask = mask__g)
        y1 = np.ma.masked_array(np.ma.log10(H.SFRSD_Ha__g * 1e6), mask = mask__g) 
        x2 = np.ma.masked_array(np.ma.log10(H.tau_V_neb__rg), mask = mask__rg)
        y2 = np.ma.masked_array(np.ma.log10(H.aSFRSD_Ha__rg * 1e6), mask = mask__rg) 
        x3 = np.ma.masked_array(np.ma.log10(tau_V_neb_norm__rg), mask = mask__rg)
        y3 = np.ma.masked_array(np.ma.log10(aSFRSD_Ha_norm__rg), mask = mask__rg) 
        x4 = np.ma.masked_array(np.ma.log10(tau_V_neb_norm_GAL__g), mask = mask_GAL__g) 
        y4 = np.ma.masked_array(np.ma.log10(SFRSD_Ha_norm_GAL__g * 1e6), mask = mask_GAL__g)
        x5 = np.ma.masked_array(np.ma.log10(H.integrated_tau_V_neb__g), mask = mask_GAL__g)
        y5 = np.ma.masked_array(np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6), mask = mask_GAL__g)
      
        x1label = r'$\log\ \tau_V^{neb}$'
        y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{neb}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'
        x4label = r'$\log\ \tau_V^{neb}(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^{H\alpha}(@1HLR)$'
        x5label = r'$\log\ \tau_V^{neb} (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^{H\alpha} (int)$'
      
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, z = H.zone_dist_HLR__g, xlabel = x1label, ylabel = y1label, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), xlabel = x2label, ylabel = y2label, z = H.Rtoplot(x2.shape).flatten(), **kwargs_zbins)    
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), xlabel = x3label, ylabel = y3label, z = H.Rtoplot(x3.shape).flatten(), **kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()

        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
      
        ############## Mixed ##############
        ############## Mixed ##############
        ############## Mixed ############## 
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)

        f.suptitle(suptitle, fontsize = 9)
           
        tau_V_neb_norm_GAL__g = (H.tau_V_neb__rg[9:11, :]).mean(axis = 0)
        tau_V_neb_norm__rg = H.tau_V_neb__rg / tau_V_neb_norm_GAL__g
        SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][9:11, :]).mean(axis = 0)
        aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
       
        x1 = np.ma.masked_array(np.ma.log10(H.tau_V_neb__g), mask = mask__g)
        y1 = np.ma.masked_array(np.ma.log10(H.SFRSD__Tg[iT] * 1e6), mask = mask__g) 
        x2 = np.ma.masked_array(np.ma.log10(H.tau_V_neb__rg), mask = mask__rg)
        y2 = np.ma.masked_array(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6), mask = mask__rg) 
        x3 = np.ma.masked_array(np.ma.log10(tau_V_neb_norm__rg), mask = mask__rg)
        y3 = np.ma.masked_array(np.ma.log10(aSFRSD_norm__rg), mask = mask__rg) 
        x4 = np.ma.masked_array(np.ma.log10(tau_V_neb_norm_GAL__g), mask = mask_GAL__g) 
        y4 = np.ma.masked_array(np.ma.log10(SFRSD_norm_GAL__g * 1e6), mask = mask_GAL__g)
        x5 = np.ma.masked_array(np.ma.log10(H.integrated_tau_V_neb__g), mask = mask_GAL__g) 
        y5 = np.ma.masked_array(np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6), mask = mask_GAL__g)
 
        x1label = r'$\log\ \tau_V^{neb}$'
        y1label = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{neb}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^{neb}(R)}{\tau_V^{neb}(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
        x4label = r'$\log\ \tau_V^{neb}(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^\star(@1HLR)$'
        x5label = r'$\log\ \tau_V^{neb} (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^\star (int)$'

        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label, z = H.zone_dist_HLR__g, **kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()
      
        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)
      
        ############## Mixed ##############
        ############## Mixed ##############
        ############## Mixed ############## 
        NRows = 3
        NCols = 4
        f = plt.figure()
        kwargs_zbins = default_kwargs_zbins.copy()
        kwargs_zbins.update(f = f)
        #f, axArr = plt.subplots(NRows, NCols)
        #page_size_inches = (NCols * 3, NRows * 1.5)
        page_size_inches = A4Size_inches
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        
        f.suptitle(suptitle, fontsize = 9)

        tau_V_norm_GAL__g = (H.tau_V__Trg[iT][9:11, :]).mean(axis = 0)
        tau_V_norm__rg = H.tau_V__Trg[iT] / tau_V_norm_GAL__g
        SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[9:11, :]).mean(axis = 0)
        aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g
       
        x1 = np.ma.masked_array(np.ma.log10(H.tau_V__Tg[iT]), mask = mask__g)
        y1 = np.ma.masked_array(np.ma.log10(H.SFRSD_Ha__g * 1e6), mask = mask__g) 
        x2 = np.ma.masked_array(np.ma.log10(H.tau_V__Trg[iT]), mask = mask__rg)
        y2 = np.ma.masked_array(np.ma.log10(H.aSFRSD_Ha__rg * 1e6), mask = mask__rg) 
        x3 = np.ma.masked_array(np.ma.log10(tau_V_norm__rg), mask = mask__rg)
        y3 = np.ma.masked_array(np.ma.log10(aSFRSD_Ha_norm__rg), mask = mask__rg) 
        x4 = np.ma.masked_array(np.ma.log10(tau_V_norm_GAL__g), mask = mask_GAL__g) 
        y4 = np.ma.masked_array(np.ma.log10(SFRSD_Ha_norm_GAL__g * 1e6), mask = mask_GAL__g)
        x5 = np.ma.masked_array(np.ma.log10(H.integrated_tau_V__g), mask = mask_GAL__g) 
        y5 = np.ma.masked_array(np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6), mask = mask_GAL__g)
      
        x1label = r'$\log\ \tau_V^{\star}$'
        y1label = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$'
        x2label = r'$\log\ \tau_V^{\star}(R)$'
        y2label = r'$\log\ \Sigma_{SFR}^{H\alpha}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        x3label = r'$\log\ \frac{\tau_V^\star(R)}{\tau_V^\star(@1HLR)}$'
        y3label = r'$\log\ \frac{\Sigma_{SFR}^{H\alpha}(R)}{\Sigma_{SFR}^{H\alpha}(@1HLR)}$'
        x4label = r'$\log\ \tau_V^\star(@1HLR)$'
        y4label = r'$\log\ \Sigma_{SFR}^{H\alpha}(@1HLR)$'
        x5label = r'$\log\ \tau_V^\star (int)$'
        y5label = r'$\log\ \Sigma_{SFR}^{H\alpha} (int)$'
       
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 0.5])
        kwargs_zbins.update(ylim = [-3.5, 1])
        kwargs_zbins.update(cb = False)
        C.plot_zbins(x = x1, y = y1, xlabel = x1label, ylabel = y1label,z = H.zone_dist_HLR__g, **kwargs_zbins)    
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        C.plot_zbins(x = x2.flatten(), y = y2.flatten(), z = H.Rtoplot(x2.shape).flatten(), xlabel = x2label, ylabel = y2label, **kwargs_zbins)
         
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        ax.set_axis_on()
        kwargs_zbins.update(ax = ax)
        kwargs_zbins.update(xlim = [-1.5, 1])
        kwargs_zbins.update(ylim = [-2.5, 2])
        kwargs_zbins.update(cb = True)
        C.plot_zbins(x = x3.flatten(), y = y3.flatten(), z = H.Rtoplot(x3.shape).flatten(), xlabel = x3label, ylabel = y3label,**kwargs_zbins)
 
        ax = plt.subplot2grid(grid_shape, loc = (0, 3))
        ax.set_axis_on()
        ax.scatter(x4, y4, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x4.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax2 = ax.twinx()
        ax.set_xlabel(x4label)
        ax2.set_ylabel(y4label)
        plt.setp(ax2.get_yticklabels(), visible = False)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x4, y4, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (1, 3))
        ax.set_axis_on()
        ax.scatter(x5, y5, **default_sc_kwargs)
        plot_text_ax(ax, 'N:%d' % x5.count() , 0.98, 0.98, 8, 'top', 'right', 'k')
        ax.set_xlim(-1.5, 0.5)
        ax.set_ylim(-3.5, 1)
        ax.set_xlabel(x5label)
        ax.set_ylabel(y5label)
        a, b, sa, sb = plotOLSbisectorAxis(ax, x5, y5, **kwargs_ols)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.grid()

        ax = plt.subplot2grid(grid_shape, loc = (2, 3))
        ax.set_axis_on()
        ax.plot(H.Rtoplot()[:, 0], x2.mean(axis = 1), label = x2label)
        ax.plot(H.Rtoplot()[:, 0], x3.mean(axis = 1), label = x3label)
        ax.plot(H.Rtoplot()[:, 0], y2.mean(axis = 1), label = y2label.replace('\ [M_\odot yr^{-1} kpc^{-2}]', ''))
        ax.plot(H.Rtoplot()[:, 0], y3.mean(axis = 1), label = y3label)
        ax.set_xlabel(r'R [HLR]')
        ax.set_ylim(-3.5, 0.5)
        ax.set_ylim(-3.5, 0.5)
        x_major_locator = 1
        x_minor_locator = 0.2
        y_major_locator = 1.
        y_minor_locator = 0.2
        ax.xaxis.set_major_locator(MultipleLocator(x_major_locator))
        ax.xaxis.set_minor_locator(MultipleLocator(x_minor_locator))
        ax.yaxis.set_major_locator(MultipleLocator(y_major_locator))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_locator))
        ax.legend(bbox_to_anchor = (1.0, 1.15), fontsize = 6, frameon = False, ncol = 4)
        #ax.legend(fontsize = 8, frameon = False)
        ax.grid()

        Nbins = 20
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        ax.set_axis_on()
        ax.hist(x1.compressed(), bins = Nbins)
        ax.set_xlabel(x1label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        ax.set_axis_on()
        ax.hist(x2.compressed(), bins = Nbins)
        ax.set_xlabel(x2label)

        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        ax.set_axis_on()
        ax.hist(x3.compressed(), bins = Nbins)
        ax.set_xlabel(x3label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 0))
        ax.set_axis_on()
        ax.hist(y1.compressed(), bins = Nbins)
        ax.set_xlabel(y1label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 1))
        ax.set_axis_on()
        ax.hist(y2.compressed(), bins = Nbins)
        ax.set_xlabel(y2label)

        ax = plt.subplot2grid(grid_shape, loc = (2, 2))
        ax.set_axis_on()
        ax.hist(y3.compressed(), bins = Nbins)
        ax.set_xlabel(y3label)

        f.subplots_adjust(hspace = 0.4, wspace = 0.4)
        pdf.savefig(f)
        plt.close(f)