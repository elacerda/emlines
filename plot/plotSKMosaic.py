#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import sys
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from CALIFAUtils.objects import H5SFRData
from CALIFAUtils.plots import OLS_bisector
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.plots import plot_gal_img_ax
from CALIFAUtils.plots import plotOLSbisectorAxis

#debug = False
debug = True

mpl.rcParams['font.size'] = 24
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['xtick.labelsize'] = 24
mpl.rcParams['ytick.labelsize'] = 24 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

sorted_args = {
    'Mcor' : 'Mcor_GAL__g',
    'McorSD' : 'McorSD_GAL__g',
    'morf' : 'morfType_GAL__g',
    'Mr' : 'Mr_GAL__g',
    'ba' : 'ba_GAL__g',
}

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
        iT = np.int(sys.argv[2])
        sort_by = sys.argv[3]
    except IndexError:
        print 'usage: %s HDF5FILE age_i %s' % (sys.argv[0], sorted_args.keys())
        exit(1)
        
    H = H5SFRData(h5file)
    
    tSF__T = H.tSF__T
    
    xOkMin = H.xOkMin
    tauVOkMin = H.tauVOkMin
    tauVNebOkMin = H.tauVNebOkMin
    tauVNebErrMax = H.tauVNebErrMax

    aSFRSD__rg = H.aSFRSD__Trg[iT]
    tau_V__rg = H.tau_V__Trg[iT]
    x = np.ma.log10(tau_V__rg.flatten())
    y = np.ma.log10(aSFRSD__rg.flatten() * 1e6)
    mask = x.mask | y.mask
    xm = np.ma.masked_array(x, mask = mask)
    ym = np.ma.masked_array(y, mask = mask)
    a_ols, b_ols, sigma_a_ols, sigma_b_ols = OLS_bisector(xm, ym)
    step = (x.max() - x.min()) / len(x)
    ALL_X = np.linspace(x.min(), x.max() + step, len(x))
    ALL_Y = a_ols * ALL_X + b_ols
    Yrms = (y - (a_ols * xm + b_ols)).std()
    Yrms_str = r' : $y_{rms}$:%.2f' % Yrms
    if b_ols > 0:
        txt_ols = r'$y_{OLS}$ = %.2f$x$ + %.2f%s' %  (a_ols, b_ols, Yrms_str)
    else:
        txt_ols = r'$y_{OLS}$ = %.2f$x$ - %.2f%s' %  (a_ols, b_ols * -1., Yrms_str)
    
    gals, arg_sorted = H.sort_gal_by_prop(sorted_args[sort_by])
    
    if debug:
        gals = gals[0:30]
    
    ###################################################################################
    xname = 'atauV'
    yname = 'aSFRSD'
    fname_suffix = 'sorted_by_%s' % sort_by
    xlabel = r'$\log\ \tau_V^{\star}(R)$'
    ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
    
    newImage = True
    NGal = len(gals)
    NRows = 4
    NCols = 8
    iGal = 0
    i = 0
    j = 0
    k = 0
    last_row = 0
    
    while iGal < NGal:
        if newImage:
            f, axArr = plt.subplots(NRows, NCols)
            f.set_size_inches((NCols * 5.34, NRows * 5.))
            plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
            plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
            f.suptitle(r'sorted by:%s  age:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f  %s' % (sort_by, tSF__T[iT]/1e6, xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax, txt_ols), fontsize = 24)
            for ax in f.axes:
                ax.set_axis_off()
            f.text(0.5, 0.02, xlabel, ha = 'center', va = 'center', fontsize = 30)
            f.text(0.03, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical', fontsize = 30)
            newImage = False
            
        gal = gals[iGal]
        
        if gal in H.califaIDs:
            last_row = i
    
            ax = axArr[i, j]
            j += 1
            ax_img = axArr[i, j]
            ax.set_axis_on()
            ax_img.set_axis_on()
    
            aSFRSD__r = getattr(H, '%s_aSFRSD__Trg' % gal)[iT]
            atau_V__r = getattr(H, '%s_tau_V__Trg' % gal)[iT]
            
            xlim = [np.log10(tauVOkMin), 0.5]
            ylim = [-3.5, 1]
            
            ##########################
            x = np.ma.log10(atau_V__r)
            y = np.ma.log10(aSFRSD__r * 1e6)
            mask = x.mask | y.mask
            xm = np.ma.masked_array(x, mask = mask)
            ym = np.ma.masked_array(y, mask = mask)
            sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'jet_r', marker = 'o', s = 50, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
            plotOLSbisectorAxis(ax, xm, ym, 0.98, 0.05, 23, color = 'b', rms = True)
            ax.plot(ALL_X, ALL_Y, c = 'grey', ls = '--', lw = 2)
            ##########################
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.125))
            ax.yaxis.set_major_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.125))
            ax.grid(which = 'major')
            
            ax = plot_gal_img_ax(ax_img, '/Users/lacerda/CALIFA/images/%s.jpg' % gal, gal, 0.02, 0.98, 30)
            txt = '%s' % str(arg_sorted[i])
            plot_text_ax(ax, txt, 0.98, 0.02, 30, 'bottom', 'right', color = 'w')
            
            if i == NRows - 1 and j == 1:
                plt.setp(ax.get_xticklabels(), visible = True, rotation = 90)
                plt.setp(ax.get_yticklabels(), visible = True)
                
            if j == NCols - 1:
                if i == NRows - 1:
                    i = 0
                    newImage = True
                    plt.subplots_adjust(wspace = 0, hspace = 0, left = 0.08, bottom = 0.05, right = 0.95, top = 0.95)
                    f.savefig('%s_%s_%s_%02d.png' % (xname, yname, fname_suffix, k))
                    plt.close(f)
                    k += 1
                else:
                    i += 1
                j = 0
            else:
                if newImage == False and iGal == NGal - 1:
                    ax1 = axArr[last_row, 0]
                    plt.setp(ax1.get_xticklabels(), visible = True, rotation = 90)
                    plt.setp(ax1.get_yticklabels(), visible = True)
                    
                    plt.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
                    f.savefig('%s_%s_%s_%02d.png' % (xname, yname, fname_suffix, k))
                    plt.close(f)
                else:
                    j += 1
            
        iGal += 1
