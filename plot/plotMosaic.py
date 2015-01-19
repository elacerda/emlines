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
                     list_gal_sorted_by_data, calcRunningStats, \
                     data_uniq, plot_text_ax, plot_linreg_params
from matplotlib.ticker import MultipleLocator
from scipy import stats as st


def plot_gal_img_ax(ax, imgfile, gal):
    galimg = plt.imread(imgfile)
    ax.imshow(galimg)
    txt = '%s' %  gal
    plot_text_ax(ax, txt, 0.05, 0.92, 8, 'top', 'left')


def plot_reglin_ax(ax, x, y, txt_x_pos, txt_y_pos, color):
    y_slope, y_intercept, y_r_value, y_p_value, y_std_err = st.linregress(x, y)

    step = (x.max() - x.min()) / len(x)
    X = np.linspace(x.min(), x.max() + step, len(x))
    Y = y_slope * X + y_intercept
    ax.plot(X, Y, c = color, ls = '-', lw = 2)
    txt = 'y = %.2fx+%.2f' %  (y_slope, y_intercept)
    plot_text_ax(ax, txt, txt_x_pos, txt_y_pos, 8, 'top', 'left', color)

    return y_slope, y_intercept, y_r_value, y_p_value, y_std_err


def new_img_mosaic(NRows, NCols, age, ticks, sorted_by):
        f, axArr = plt.subplots(NRows, NCols)
        f.set_dpi(300)
        f.set_size_inches(11.69,8.27) 
        plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
        plot_suptitle = '%.2f Myr' % (age/1e6)
        if ticks:
            plot_suptitle = '%s xlim fixed' % plot_suptitle
        plot_suptitle = '%s - sort_param: %s' % (plot_suptitle, sorted_by)             
        f.suptitle(plot_suptitle)
        for ax in f.axes:
            ax.set_axis_off()
        return f, axArr

    
def save_img_mosaic(f, fname, ticks):
    plt.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
    if ticks:
        fname = '%s_xfix' % (fname)
    f.savefig('%s.png' % fname)


sort_param = {
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
    print 'usage: %s HDF5FILE [0-39] %s' % (sys.argv[0], sort_param.keys())
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
SFR__Tg = get_attrib_h5(h5, 'SFR__Tg')
SFR_Ha__g = get_attrib_h5(h5, 'SFR_Ha__g')
SFRSD__Tg = get_attrib_h5(h5, 'SFRSD__Tg')
SFRSD_Ha__g = get_attrib_h5(h5, 'SFRSD_Ha__g')
SFRSD_Ha_kpc__g = get_attrib_h5(h5, 'SFRSD_Ha_kpc__g')
dist_zone__g = get_attrib_h5(h5, 'dist_zone__g')
tau_V__Tg = get_attrib_h5(h5, 'tau_V__Tg')
tau_V_neb__g = get_attrib_h5(h5, 'tau_V_neb__g')
L_int_Ha__g = get_attrib_h5(h5, 'L_int_Ha__g')
F_obs_Ha__g = get_attrib_h5(h5, 'F_obs_Ha__g')
Mcor__g = get_attrib_h5(h5, 'Mcor__g')
McorSD__g = get_attrib_h5(h5, 'McorSD__g')

# galaxy wide quantities replicated by zones
Mcor_GAL_zones__g = get_attrib_h5(h5, 'Mcor_GAL_zones__g')
McorSD_GAL_zones__g = get_attrib_h5(h5, 'McorSD_GAL_zones__g')
morfType_GAL_zones__g = get_attrib_h5(h5, 'morfType_GAL_zones__g')
at_flux_GAL_zones__g = get_attrib_h5(h5, 'at_flux_GAL_zones__g')
califaID_GAL_zones__g = get_attrib_h5(h5, 'califaID_GAL_zones__g')
Mr_GAL_zones__g = get_attrib_h5(h5, 'Mr_GAL_zones__g')
ur_GAL_zones__g = get_attrib_h5(h5, 'ur_GAL_zones__g')

sort_param['Mcor'] = np.log10(Mcor_GAL_zones__g)
sort_param['McorSD'] = np.log10(McorSD_GAL_zones__g)
sort_param['MorphType'] = morfType_GAL_zones__g
sort_param['Mr'] = Mr_GAL_zones__g
sort_param['u-r'] = ur_GAL_zones__g

label = { '%s' % k : False for k in sort_param.keys()} 

label['Mcor'] = r'$\log\ M_\star^{gal}\ [M_\odot]$'
label['McorSD'] = r'$\log\ \mu_\star^{gal}\ [M_\odot\ pc^{-2}]$'
label['MorphType'] = r'Morphological type'
label['Mr'] = r'$M_r$'
label['u-r'] = r'u - r'

NGal, listGal_ns, sorted_data__g = data_uniq(califaID_GAL_zones__g, sort_param[sorted_by])
listGal = list_gal_sorted_by_data(listGal_ns, sorted_data__g, -1)
#ticks = True
ticks = False
min_pixel_to_plot = 5
fname_suffix = 'sor%s_%.2fMyr' % (sorted_by, tSF__T[iT]/1e6)

###################################################################################
xname = 'dtau'
yname = 'SFR_Ha'

newImage = True
NRows = 7
NCols = 10
iGal = 0
i = 0
j = 0
k = 0
last_row = 0

dtau_slope = np.ma.masked_all((NGal))
dtau_x0 = np.ma.masked_all((NGal))
dtau_stderr = np.ma.masked_all((NGal))
dtau_r = np.ma.masked_all((NGal))

while iGal < NGal:
    if newImage:
        f, axArr = new_img_mosaic(NRows, NCols, tSF__T[iT], ticks, sorted_by)
        newImage = False
        
    gal = listGal[iGal]
    where_slice = np.where(califaID_GAL_zones__g == gal)[0]
    N_zone = len(where_slice)
    
    x1 = tau_V__Tg[iT][where_slice]
    x2 = tau_V_neb__g[where_slice]
    y = np.ma.log10(SFR_Ha__g[where_slice])
    
    mask = ~(x1.mask | x2.mask | y.mask)
    x1m = x1[mask]
    x2m = x2[mask]
    ym = y[mask]
    xm = x2m - x1m    

    N_not_masked = mask.sum()
     
    if N_not_masked >= min_pixel_to_plot:
        last_row = i
        ax1 = axArr[i, j]
        j += 1
        ax2 = axArr[i, j]
        ax1.set_axis_on()
        ax2.set_axis_on()

        #print '%s %d %d' % (gal, N_zone, N_not_masked)
        
        xlabel = r'$\delta\tau\ [\tau_V^{neb}\ -\ \tau_V^\star]$'
        ylabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'

        ax1.scatter(xm, ym, c = 'k', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6)
        xlim = [ -1,  2. ]    
        ylim = [-4., -0.5]
        
        ax1.set_xlim(xlim)
        
        if ticks:
            ax1.set_ylim(ylim)
        #ax.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
        #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')
        
        # notice the change in y and x axis so:
        # log SFR = A log_dtau + B 
        aux = plot_reglin_ax(ax1, xm, ym, 0.05, 0.92, 'k')
        dtau_slope[iGal] = aux[0]
        dtau_x0[iGal] = aux[1]
        dtau_r[iGal] = aux[2]
        dtau_stderr[iGal] = aux[4]

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
                save_img_mosaic(f, '%s_%s_%s_%d' % (xname, yname, fname_suffix, k), ticks)
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
        
        save_img_mosaic(f, '%s_%s_%s_%d' % (xname, yname, fname_suffix, k), ticks)
        plt.close(f)
        
    iGal += 1

ylab = r'$\delta\tau$'
plot_linreg_params(dtau_slope, sorted_data__g, label[sorted_by], 
                   r'%s slope' % ylab, '%s_slope_%s.png' % (xname, fname_suffix),
                   best_param = 1., fontsize = 8)
plot_linreg_params(dtau_x0, sorted_data__g, label[sorted_by], 
                   r'%s x0' % ylab, '%s_x0_%s.png' % (xname, fname_suffix), 
                   best_param = 0., fontsize = 8)
plot_linreg_params(dtau_stderr, sorted_data__g, label[sorted_by], 
                   r'%s stderr' % ylab, '%s_stderr_%s.png' % (xname, fname_suffix))
plot_linreg_params(dtau_r**2., sorted_data__g, label[sorted_by], 
                   r'%s $r^2$' % ylab, '%s_sqrcor_%s.png' % (xname, fname_suffix),
                   best_param = 1., fontsize = 8)

##########################################################################################
x1name = 'tau_V'
x2name = 'tau_V_neb'
yname = 'SFR_Ha'
newImage = True
NRows = 7
NCols = 10
iGal = 0
i = 0
j = 0
k = 0
min_pixel_to_plot = 5
last_row = 0

tau_V_slope = np.ma.masked_all((NGal))
tau_V_x0 = np.ma.masked_all((NGal))
tau_V_stderr = np.ma.masked_all((NGal))
tau_V_r = np.ma.masked_all((NGal))
tau_V_neb_slope = np.ma.masked_all((NGal))
tau_V_neb_x0 = np.ma.masked_all((NGal))
tau_V_neb_stderr = np.ma.masked_all((NGal))
tau_V_neb_r = np.ma.masked_all((NGal))
 
while iGal < NGal:
    if newImage:
        f, axArr = new_img_mosaic(NRows, NCols, tSF__T[iT], ticks, sorted_by)
        newImage = False
         
    gal = listGal[iGal]
    where_slice = np.where(califaID_GAL_zones__g == gal)[0]
    N_zone = len(where_slice)
     
    x1 = tau_V__Tg[iT][where_slice]
    x2 = tau_V_neb__g[where_slice]
    y = np.ma.log10(SFR_Ha__g[where_slice])
     
    mask = ~(y.mask | x1.mask | x2.mask)
    x1m = x1[mask]
    x2m = x2[mask]    
    ym = y[mask]
 
    N_not_masked = mask.sum()
      
    if N_not_masked >= min_pixel_to_plot:
        last_row = i
        ax1 = axArr[i, j]
        j += 1
        ax2 = axArr[i, j]
        ax1.set_axis_on()
        ax2.set_axis_on()

        #print '%s %d %d' % (gal, N_zone, N_not_masked)
        
        xlabel = r'$\tau_V$'
        ylabel = r'$\log\ SFR_{neb}\ $[M${}_\odot$ yr${}^{-1}]$'
 
        ax1.scatter(x1m, ym, c = 'r', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^\star$')
        ax1.scatter(x2m, ym, c = 'b', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6, label = r'$\tau_V^{neb}$')
        xlim = [ 0.,  2. ]    
        ylim = [-4., -0.5]

        ax1.set_xlim(xlim)
        
        if ticks:
            ax1.set_ylim(ylim)
        
        #ax.set_title(r'%.2f Myr' % (tSF__T[iT] / 1.e6))
        #ax.legend(fontsize = 14, frameon = False, loc = 'upper right')
         
        aux = plot_reglin_ax(ax1, x1m, ym, 0.05, 0.82, 'r')
        tau_V_slope[iGal] = aux[0]
        tau_V_x0[iGal] = aux[1]
        tau_V_r[iGal] = aux[2]
        tau_V_stderr[iGal] = aux[4]
         
        aux = plot_reglin_ax(ax1, x2m, ym, 0.05, 0.92, 'b')            
        tau_V_neb_slope[iGal] = aux[0]
        tau_V_neb_x0[iGal] = aux[1]
        tau_V_neb_r[iGal] = aux[2]
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
                save_img_mosaic(f, '%s_%s_%s_%s_%d' % (x1name, x2name, yname, fname_suffix, k), ticks)
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
         
        save_img_mosaic(f, '%s_%s_%s_%s_%d' % (x1name, x2name, yname, fname_suffix, k), ticks)
        plt.close(f)
         
    iGal += 1

y1lab = r'$\tau_V^\star$'
plot_linreg_params(tau_V_slope, sorted_data__g, label[sorted_by], 
                   r'%s slope' % y1lab, '%s_slope_%s.png' % (x1name, fname_suffix),
                   best_param = 1., fontsize = 8)
plot_linreg_params(tau_V_x0, sorted_data__g, label[sorted_by], 
                   r'%s x0' % y1lab, '%s_x0_%s.png' % (x1name, fname_suffix), 
                   best_param = 0., fontsize = 8)
plot_linreg_params(tau_V_r**2., sorted_data__g, label[sorted_by], 
                   r'%s $r^2$' % y1lab, '%s_sqrcor_%s.png' % (x1name, fname_suffix),
                   best_param = 1., fontsize = 8)
plot_linreg_params(tau_V_stderr, sorted_data__g, label[sorted_by], 
                   r'%s stderr' % y1lab, '%s_stderr_%s.png' % (x1name, fname_suffix))

y2lab = r'$\tau_V^{neb}$' 
plot_linreg_params(tau_V_neb_slope, sorted_data__g, label[sorted_by], 
                   r'%s slope' % y2lab, '%s_slope_%s.png' % (x2name, fname_suffix),
                   best_param = 1., fontsize = 8)
plot_linreg_params(tau_V_neb_x0, sorted_data__g, label[sorted_by], 
                   r'%s x0' % y2lab, '%s_x0_%s.png' % (x2name, fname_suffix), 
                   best_param = 0., fontsize = 8)
plot_linreg_params(tau_V_neb_r**2., sorted_data__g, label[sorted_by], 
                   r'%s $r^2$' % y2lab, '%s_sqrcor_%s.png' % (x2name, fname_suffix),
                   best_param = 1., fontsize = 8)
plot_linreg_params(tau_V_neb_stderr, sorted_data__g, label[sorted_by], 
                   r'%s stderr' % y2lab, '%s_stderr_%s.png' % (x2name, fname_suffix))
