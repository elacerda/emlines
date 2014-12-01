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
from scipy import stats as st


def plot_gal_img_ax(ax, imgfile, gal):
    galimg = plt.imread(imgfile)
    ax.imshow(galimg)
    txt = '%s' %  gal
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.05, 0.92, txt, fontsize = 8, color = 'white',
            transform = ax.transAxes,
            verticalalignment = 'top', horizontalalignment = 'left',
            bbox = textbox)

def plot_reglin_ax(ax, SFR, tau_V):
    step = (SFR.max() - SFR.min()) / len(SFR)
    X = np.linspace(SFR.min(), SFR.max() + step, len(SFR))

    y_slope, y_intercept, y_r_value, y_p_value, y_std_err = st.linregress(SFR, tau_V)
    ax.plot(X, y_slope * X + y_intercept, 'b-', lw = 2)
    txt = 'y = %.2fx+%.2f' %  (y_slope, y_intercept)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.05, 0.82, txt, fontsize = 8, color = 'blue',
            transform = ax.transAxes,
            verticalalignment = 'top', horizontalalignment = 'left',
            bbox = textbox)

    return y_slope, y_intercept, y_r_value, y_p_value, y_std_err


def plot_tau_SFR_ax(ax, SFR, tau_V, tau_V_neb, SFRlim, taulim, ticks): 
    ax.scatter(SFR, tau_V, c = 'b', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^\star$')
    ax.scatter(SFR, tau_V_neb, c = 'r', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^{neb}$')
        
    ax.set_ylim(taulim)
    if ticks:
        ax.set_xlim(SFRlim)
    #ax.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
    #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')


def new_img_mosaic(NRows, NCols, age, ticks, sorted_by):
        f, axArr = plt.subplots(NRows, NCols)
        f.set_dpi(300)
        f.set_size_inches(11.69,8.27) 
        plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
        plot_suptitle = '%.2f Myr' % (age/1e6)
        if ticks:
            plot_suptitle = '%s xlim fixed' % plot_suptitle
        plot_suptitle = '%s - sorted by %s' % (plot_suptitle, sorted_by)             
        f.suptitle(plot_suptitle)
        for ax in f.axes:
            ax.set_axis_off()
        return f, axArr

    
def save_img_mosaic(f, fname, ticks):
    plt.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
    if ticks:
        fname = '%s_xfix' % (fname)
    f.savefig('%s.png' % fname)

sorted = {
    'Mcor' : False,
    'McorSD' : False,
    'MorphType' : False,
    'Mr' : False,
    'u-r' : False,
}

try:
    h5file = sys.argv[1]
    iT = np.int(sys.argv[2])
    sorted_by = sys.argv[3]
except IndexError:
    print 'usage: %s HDF5FILE [0-39] %s' % (sys.argv[0], sorted.keys())
    exit(1)

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

h5 = h5py.File(h5file, 'r')

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

sorted['Mcor'] = ALL_Mcor_GAL_zones__g
sorted['McorSD'] = ALL_McorSD_GAL_zones__g
sorted['MorphType'] = ALL_morfType_GAL_zones__g
sorted['Mr'] = ALL_Mr_GAL_zones__g
sorted['u-r'] = ALL_ur_GAL_zones__g

listGal = list_gal_sorted_by_data(ALL_califaID_GAL_zones__g, sorted[sorted_by], 0)
NGal = len(listGal)

fname = 'SFRNeb_tauV_tauVNeb_sor%s' % sorted_by
newImage = True
NRows = 7
NCols = 10
iGal = 0
i = 0
j = 0
k = 0
min_pixel_to_plot = 5
last_row = 0
ticks = True
#ticks = False
tau_V_slope = np.ma.masked_all((NGal))
tau_V_x0 = np.ma.masked_all((NGal))
tau_V_stderr = np.ma.masked_all((NGal))
tau_V_neb_slope = np.ma.masked_all((NGal))
tau_V_neb_x0 = np.ma.masked_all((NGal))
tau_V_neb_stderr = np.ma.masked_all((NGal))

while iGal < NGal:
    if newImage:
        f, axArr = new_img_mosaic(NRows, NCols, tSF__T[iT], ticks, sorted_by)
        newImage = False
        
    gal = listGal[iGal]
    where_slice = np.where(ALL_califaID_GAL_zones__g == gal)[0]
    N_zone = len(where_slice)
    
    x = np.ma.log10(ALL_SFR_Ha__g[where_slice])
    y1 = ALL_tau_V__Tg[iT][where_slice]
    y2 = ALL_tau_V_neb__g[where_slice]
    
    mask = ~(x.mask | y1.mask | y2.mask)
    xm = x[mask]
    y1m = y1[mask]
    y2m = y2[mask]    

    N_not_masked = mask.sum()
     
    if N_not_masked >= min_pixel_to_plot:
        last_row = i
        ax1 = axArr[i, j]
        j += 1
        ax2 = axArr[i, j]
        ax1.set_axis_on()
        ax2.set_axis_on()

        print '%s %d %d' % (gal, N_zone, N_not_masked)
        
        xlabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'
        ylabel = r'$\tau_V$'

        plot_tau_SFR_ax(ax1, xm, y1m, y2m, [-4., -0.5], [ 0.,  2. ], ticks)
        
        aux = plot_reglin_ax(ax1, xm, y1m)
        tau_V_slope[iGal] = aux[0]
        tau_V_x0[iGal] = aux[1]
        tau_V_stderr[iGal] = aux[4]
        
        aux = plot_reglin_ax(ax1, xm, y2m)            
        tau_V_neb_slope[iGal] = aux[0]
        tau_V_neb_x0[iGal] = aux[1]
        tau_V_neb_stderr[iGal] = aux[4]

        plot_gal_img_ax(ax2, '/Users/lacerda/CALIFA/images/%s.jpg' % gal, gal)
        
        if ticks and i == NRows - 1 and j == 1:
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
                save_img_mosaic(f, '%s_%.2fMyr_%d' % (fname, tSF__T[iT]/1e6, k), ticks)
                plt.close(f)
                k += 1
            else:
                i += 1
            j = 0
        else:
            j += 1
    
    if newImage == False and iGal == NGal - 1:
        ax1 = axArr[last_row, 0]
        if ticks:
            plt.setp(ax1.get_xticklabels(), visible = True, rotation = 90)
            plt.setp(ax1.get_yticklabels(), visible = True)
        
        if last_row < NRows - 1:
            ax = axArr[last_row, 5]
            ax.set_xlabel(xlabel)
        
            if last_row < 3:
                ax1.set_ylabel(ylabel)
        
        save_img_mosaic(f, '%s_%.2fMyr_%d' % (fname, tSF__T[iT]/1e6, k), ticks)
        plt.close(f)
        
    iGal += 1
