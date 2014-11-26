#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5, density_contour, \
                     list_gal_sorted_by_data, calcRunningStats
from matplotlib.ticker import MultipleLocator

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')

# zones
ALL_SFR__Tg = get_attrib_h5(h5, 'ALL_SFR__Tg')
ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
ALL_SFRSD__Tg = get_attrib_h5(h5, 'ALL_SFRSD__Tg')
ALL_SFRSD_Ha__g = get_attrib_h5(h5, 'ALL_SFRSD_Ha__g')
ALL_SFRSD_Ha_kpc__g = get_attrib_h5(h5, 'ALL_SFRSD_Ha_kpc__g')
ALL_dist_zone__g = get_attrib_h5(h5, 'ALL_dist_zone__g')
ALL_tau_V__Tg = get_attrib_h5(h5, 'ALL_tau_V__Tg')
ALL_tau_V_neb__g = get_attrib_h5(h5, 'ALL_tau_V_neb__g')
ALL_L_int_Ha__g = get_attrib_h5(h5, 'ALL_L_int_Ha__g')
ALL_F_obs_Ha__g = get_attrib_h5(h5, 'ALL_F_obs_Ha__g')
ALL_Mcor__g = get_attrib_h5(h5, 'ALL_Mcor__g')
ALL_McorSD__g = get_attrib_h5(h5, 'ALL_McorSD__g')

# galaxy wide quantities replicated by zones
ALL_Mcor_GAL_zones__g = get_attrib_h5(h5, 'ALL_Mcor_GAL_zones__g')
ALL_McorSD_GAL_zones__g = get_attrib_h5(h5, 'ALL_McorSD_GAL_zones__g')
ALL_morfType_GAL_zones__g = get_attrib_h5(h5, 'ALL_morfType_GAL_zones__g')
ALL_at_flux_GAL_zones__g = get_attrib_h5(h5, 'ALL_at_flux_GAL_zones__g')
ALL_califaID_GAL_zones__g = get_attrib_h5(h5, 'ALL_califaID_GAL_zones__g')
ALL_Mr_GAL_zones__g = get_attrib_h5(h5, 'ALL_Mr_GAL_zones__g')
ALL_ur_GAL_zones__g = get_attrib_h5(h5, 'ALL_ur_GAL_zones__g')

listGal = list_gal_sorted_by_data(ALL_califaID_GAL_zones__g, ALL_Mcor_GAL_zones__g, 0)
NGal = len(listGal)

fname = 'SFRNeb_tauV_tauVNeb'
newImage = True
NRows = 7
NCols = 10
iGal = 0
i = 0
j = 0
k = 0
min_pixel_to_plot = 5

while iGal < NGal:
    if newImage:
        f, axArr = plt.subplots(NRows, NCols)
        f.set_dpi(300)
        f.set_size_inches(11.69,8.27) 
        plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
        for ax in f.axes:
            ax.set_axis_off()
        newImage = False
        
    gal = listGal[iGal]
    where_slice = np.where(ALL_califaID_GAL_zones__g == gal)[0]
    N_zone = len(where_slice)
    
    SFR_Ha__z = ALL_SFR_Ha__g[where_slice] 
    
    tau_V_neb__z = ALL_tau_V_neb__g[where_slice]

    tau_V__Tz = np.ma.masked_all((len(tSF__T), N_zone))
    for iT, tSF in enumerate(tSF__T):
        tau_V__g = ALL_tau_V__Tg[iT]
        tau_V__Tz[iT, :] = tau_V__g[where_slice]
        
    iT = 4
    x = np.ma.log10(SFR_Ha__z)
    y1 = tau_V__Tz[iT, :]
    y2 = tau_V_neb__z
    mask = ~(x.mask | y1.mask | y2.mask)
    xm = x[mask]
    y1m = y1[mask]
    y2m = y2[mask]    

    N_not_masked = mask.sum()
     
    if N_not_masked >= min_pixel_to_plot: 
        print '%s %d %d' % (gal, N_zone, N_not_masked)
            
        xlabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'
        ylabel = r'$\tau_V$'

        ax1 = axArr[i, j]
        j += 1
        ax2 = axArr[i, j]
        
        ax1.set_axis_on()
        ax2.set_axis_on()
        
        ax1.scatter(xm, y1m, c = 'b', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^\star$')
        ax1.scatter(xm, y2m, c = 'r', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^{neb}$')
        ax1.set_ylim(0., 2.)
        ax1.set_xlim(-4., -0.5)
        #ax1.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
        #ax1.legend(fontsize = 14, frameon = False, loc = 'upper right')
        
        imgfile = '/Users/lacerda/CALIFA/images/%s.jpg' % gal
        
        galimg = plt.imread(imgfile)
        ax2.imshow(galimg)
        txt = '%s' %  gal
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax2.text(0.05, 0.92, txt, fontsize = 8, color = 'white',
                 transform = ax2.transAxes,
                 verticalalignment = 'top', horizontalalignment = 'left',
                 bbox = textbox)
    
        if i == NRows - 1 and j == 1:
            plt.setp(ax1.get_xticklabels(), visible = True, rotation = 90)
            plt.setp(ax1.get_yticklabels(), visible = True)
            
        if i == NRows - 1 and j == 5:
            ax2.set_xlabel(xlabel)
            
        if i == 3 and j == 1:
            ax1.set_ylabel(ylabel)
        
        if j == NCols - 1:
            if i == NRows -1:
                i = 0
                newImage = True
                plt.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
                f.savefig('%s_%d.png' %(fname, k))
                plt.close(f)
                k += 1
            else:
                i += 1
            j = 0
        else:
            j += 1
    
    if newImage == False and iGal == NGal - 1:
        plt.subplots_adjust(wspace=0, hspace=0, left=0.05, bottom=0.05, right=0.95, top=0.95)
        f.savefig('%s_%d.png' %(fname, k))
        plt.close(f)
        
    iGal += 1
