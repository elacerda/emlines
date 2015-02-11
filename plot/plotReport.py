#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys
from plot_aux import H5SFRData, density_contour, plot_text_ax, \
                     plot_linreg_params, OLS_bisector, \
                     plotOLSbisectorAxis
from scipy import stats as st

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)
    
H = H5SFRData(h5file)    

tSF__T              = H.tSF__T[0:20]
RbinCenter__r       = H.RbinCenter__r
RRange              = H.RRange
RColor              = H.RColor
SFR__Tg             = H.get_data_h5('SFR__Tg')
correl_SFR__T       = H.get_data_h5('correl_SFR__T')
correl_SFRSD__T     = H.get_data_h5('correl_SFRSD__T')
correl_aSFRSD__rT   = H.get_data_h5('correl_aSFRSD__rT')
SFR_Ha__g           = H.get_data_h5('SFR_Ha__g')
aSFRSD__Trg         = H.get_data_h5('aSFRSD__Trg')
aSFRSD_Ha__rg       = H.get_data_h5('aSFRSD_Ha__rg')
SFRSD__Tg           = H.get_data_h5('SFRSD__Tg')
SFRSD_Ha__g         = H.get_data_h5('SFRSD_Ha__g')
dist_zone__g        = H.get_data_h5('dist_zone__g')

NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(300)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
   
xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
   
iT = 0
  
a = np.ones_like(tSF__T)
b = np.ones_like(tSF__T)
b2 = np.ones_like(tSF__T)
   
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFR__Tg[iT])
        y = np.ma.log10(SFR_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        age = tSF__T[iT]
        print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-6, 0]
        yran = [-6, 0]
        scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
        
        a[iT], b[iT] = plotOLSbisectorAxis(ax, xm, ym, 0.92, 0.05, 8)

        b2[iT] = (ym - xm).mean()
        Y2 = xm + b2[iT]
        Yrms = (ym - Y2).std()
        ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
        txt = r'y = x + (%.2f) $Y_{rms}$:%.2f' %  (b2[iT], Yrms)
        plot_text_ax(ax, txt, 0.92, 0.15, 8, 'bottom', 'right', color = 'b')
        
        txt = '%.2f Myr' % (age / 1e6)
        plot_text_ax(ax, txt, 0.05, 0.92, 8, 'top', 'left')
  
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
               
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
               
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
           
        iT += 1
   
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('SFR_linregress_report.png')
plt.close(f)
  
x = np.log10(tSF__T)
xlabel = r'$\log\ t_\star$ [yr]'

plot_linreg_params(a, x, xlabel, 
                   'a', 'SFR_linregress_slope_age.png', 1., 16) 
plot_linreg_params(b, x, xlabel, 
                   r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep_age.png', 0., 16)
plot_linreg_params(b2, x, xlabel, 
                   r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep2_age.png', 0., 16)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# plot_linreg_params(sigma, x, xlabel, 
#                    r'$\sigma$', 'SFR_linregress_sigma_age.png')
# plot_linreg_params(r**2., x, xlabel, 
#                    r'$r^2$', 'SFR_linregress_sqrcorr_age.png', 1., 16) 
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(300)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
      
xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    
iT = 0
    
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFR__Tg[iT])
        y = np.ma.log10(SFR_Ha__g)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        age = tSF__T[iT]
        print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-6, 0]
        yran = [-6, 0]
        scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
        binsx = np.linspace(-6., 0., 31)
        binsy = np.linspace(min(ym), max(ym), 31)
        density_contour(xm, ym, binsx, binsy, ax = ax)
        txt = '%.2f Myr' % (age / 1e6)
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.05, 0.92, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top', horizontalalignment = 'left',
                bbox = textbox)
        txt = '$R_p:%.2f$ $R_S:%.2f$' %  (st.pearsonr(xm, ym)[0], correl_SFR__T[iT])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.92, 0.05, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'bottom', horizontalalignment = 'right',
                bbox = textbox)
        #ax.grid()
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
                
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
                
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
            
        iT += 1
    
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('SFR_report.png')
plt.close(f)
   
NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(300)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
   
xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
   
iT = 0
   
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFR__Tg[iT])
        y = np.ma.log10(SFR_Ha__g)
        z = dist_zone__g
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        zm = z[~mask]
        age = tSF__T[iT]
        print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-6, 0]
        yran = [-6, 0]
        scat = ax.scatter(xm, ym, c = zm, cmap = 'hot_r', vmax = 2., marker = 'o', s = 0.3, edgecolor = 'none') #, alpha = 0.4)
        binsx = np.linspace(-6., 0., 31)
        binsy = np.linspace(min(ym), max(ym), 31)
        density_contour(xm, ym, binsx, binsy, ax = ax)
        txt = '%.2f Myr' % (age / 1e6)
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.05, 0.92, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top', horizontalalignment = 'left',
                bbox = textbox)
        txt = '$R_p:%.2f$ $R_S:%.2f$' %  (st.pearsonr(xm, ym)[0], correl_SFR__T[iT])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.92, 0.05, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'bottom', horizontalalignment = 'right',
                bbox = textbox)
        #ax.grid()
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
               
        if i == NRows - 1 and j == 4:
            ax.set_xlabel(xlabel)
               
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
           
        iT += 1
   
f.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.05, right=0.85, top = 0.95)
#f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
cbar_ax = f.add_axes([0.88, 0.1, 0.05, 0.85])
cb = f.colorbar(scat, ticks=[0, 0.5, 1., 1.5, 2.], cax=cbar_ax)
cb.ax.set_yticklabels(['0 HLR', '0.5 HLR', '1 HLR', '1.5 HLR', '2 HLR'])
   
f.savefig('SFR_dist_report.png')
plt.close(f)
   
NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
  
xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
  
iT = 0
a = np.ones_like(tSF__T)
b = np.ones_like(tSF__T)
b2 = np.ones_like(tSF__T)
  
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
        y = np.ma.log10(SFRSD_Ha__g  * 1e6)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        age = tSF__T[iT]
        print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-3.5, 1]
        yran = [-3.5, 1]
        scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)

        a[iT], b[iT] = plotOLSbisectorAxis(ax, xm, ym, 0.92, 0.05, 8)

        b2[iT] = (ym - xm).mean()
        Y2 = xm + b2[iT]
        Yrms = (ym - Y2).std()
        ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
        txt = r'y = x + (%.2f) $Y_{rms}$:%.2f' %  (b2[iT], Yrms)
        plot_text_ax(ax, txt, 0.92, 0.15, 8, 'bottom', 'right', color = 'b')
        
        txt = '%.2f Myr' % (age / 1e6)
        plot_text_ax(ax, txt, 0.05, 0.92, 8, 'top', 'left')
 
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
              
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
              
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
          
        iT += 1
  
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('SFRSD_linregress_report.png')
plt.close(f)
 
x = np.log10(tSF__T)
xlabel = r'$\log\ t_\star$ [yr]'

plot_linreg_params(a, x, xlabel, 
                   'a', 'SFRSD_linregress_slope_age.png', 1., 16) 
plot_linreg_params(b, x, xlabel, 
                   r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep_age.png', 0., 16)
plot_linreg_params(b2, x, xlabel, 
                   r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep2_age.png', 0., 16)

NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
   
xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
   
iT = 0
   
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
        y = np.ma.log10(SFRSD_Ha__g * 1e6)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        age = tSF__T[iT]
        print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-3.5, 1]
        yran = [-3.5, 1]
        scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
        binsx = np.linspace(-6., 0., 31)
        binsy = np.linspace(min(ym), max(ym), 31)
        density_contour(xm, ym, binsx, binsy, ax = ax)
        txt = '%.2f Myr' % (age / 1e6)
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.05, 0.92, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top', horizontalalignment = 'left',
                bbox = textbox)
        txt = '$R_p:%.2f$ $R_S:%.2f$' %  (st.pearsonr(xm, ym)[0], correl_SFRSD__T[iT])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.92, 0.05, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'bottom', horizontalalignment = 'right',
                bbox = textbox)
        #ax.grid()
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
               
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
               
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
           
        iT += 1
   
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('SFRSD_report.png')
plt.close(f)
   
NRows = 4
NCols = 5
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(11.69,8.27) 
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
   
xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
   
iT = 0
   
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j] 
        x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
        y = np.ma.log10(SFRSD_Ha__g * 1e6)
        z = dist_zone__g
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        zm = z[~mask]
        age = tSF__T[iT]
        print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
        xran = [-3.5, 1]
        yran = [-3.5, 1]
        scat = ax.scatter(xm, ym, c = zm, cmap = 'hot_r', vmax = 2., marker = 'o', s = 0.3, edgecolor = 'none') #, alpha = 0.4)
        binsx = np.linspace(-6., 0., 31)
        binsy = np.linspace(min(ym), max(ym), 31)
        density_contour(xm, ym, binsx, binsy, ax = ax)
        txt = '%.2f Myr' % (age / 1e6)
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.05, 0.92, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top', horizontalalignment = 'left',
                bbox = textbox)
        txt = '$R_p:%.2f$ $R_S:%.2f$' %  (st.pearsonr(xm, ym)[0], correl_SFRSD__T[iT])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.92, 0.05, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'bottom', horizontalalignment = 'right',
                bbox = textbox)
        #ax.grid()
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
               
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
               
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
           
        iT += 1
   
f.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.05, right=0.85, top = 0.95)
#f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
cbar_ax = f.add_axes([0.88, 0.1, 0.05, 0.85])
cb = f.colorbar(scat, ticks=[0, 0.5, 1., 1.5, 2.], cax=cbar_ax)
cb.ax.set_yticklabels(['0 HLR', '0.5 HLR', '1 HLR', '1.5 HLR', '2 HLR'])
   
f.savefig('SFRSD_dist_report.png')
plt.close(f)
   
NRows = 4
NCols = 5
    
pos_y_ini = 0.38
pos_step = 0.09
Rfontsize = 12
      
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(11.69,8.27)
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
       
xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
       
NAxes = len(f.axes) 
iT = 0
jump = 0
           
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j]
    
        age = tSF__T[iT]
        n_mask = n_tot = 0
         
        for iR, RUp in enumerate(RRange):
            if iR == 0:
                RMask = RbinCenter__r <= RUp
                txt = 'R <= %.1f HLR' % RUp
            else:
                RDown = RRange[iR - 1]
                RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                txt = '%.1f < R <= %.1f HLR' % (RDown, RUp)
                     
            c = RColor[iR] 
            x = np.ma.log10(aSFRSD__Trg[iT, RMask, :].flatten() * 1e6)
            y = np.ma.log10(aSFRSD_Ha__rg[RMask, :].flatten() * 1e6)
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            n_mask += mask.sum()
            n_tot += len(x)
   
            if i == 0 and j == 0:
                if iR > 0:
                    jump = 0.25 * (iR - 1.)
                pos_x = iR + jump
                #pos_y = pos_y_ini - (iR * pos_step)
                pos_y = 1.8
                textbox = dict(alpha = 0.)
                ax.text(pos_x, 1.1, txt,
                        fontsize = Rfontsize, color = c,
                        transform = ax.transAxes,
                        va = 'top', ha = 'left',
                        bbox = textbox)
               
            scat = ax.scatter(xm, ym, c = c, marker = 'o', s = 1., edgecolor = 'none', alpha = 1.)
           
        print 'SigmaSFR x SigmaSFR_Ha Age: %.2f Myr: masked %d points of %d (Total: %d)' % (age / 1e6, n_mask, n_tot, n_tot - n_mask)
        #ax.legend(loc = 'lower left', fontsize = 12, frameon = False)
        age = tSF__T[iT]
        xran = [-3.5, 1.]
        yran = [-3.5, 1.]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # binsx = np.linspace(-4.5, 1., 51)
        # binsy = np.linspace(min(ym),max(ym), 51)
        # density_contour(xm, ym, binsx, binsy, ax=ax)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        txt = '%.2f Myr' % (age / 1e6)
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.05, 0.92, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top', horizontalalignment = 'left',
                bbox = textbox)
        txt = '$R_p:%.2f$ $R_S:%.2f$' %  (st.pearsonr(xm, ym)[0], correl_aSFRSD__rT[0, iT])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        ax.text(0.92, 0.05, txt, fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'bottom', horizontalalignment = 'right',
                bbox = textbox)
        #ax.grid()
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
               
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
               
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
   
        iT += 1
   
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('aSFRSD_report.png')
plt.close(f)
 
NRows = 4
NCols = 5
   
pos_y_ini = 0.38
pos_step = 0.09
Rfontsize = 12
     
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(11.69,8.27)
plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
      
xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
      
NAxes = len(f.axes) 
iT = 0
jump = 0
 
a = np.ones_like(tSF__T)
b = np.ones_like(tSF__T)
b2 = np.ones_like(tSF__T)
          
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j]
   
        age = tSF__T[iT]
        n_mask = n_tot = 0
        
        for iR, RUp in enumerate(RRange):
            if iR == 0:
                RMask = RbinCenter__r <= RUp
                txt = 'R <= %.1f HLR' % RUp
            else:
                RDown = RRange[iR - 1]
                RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                txt = '%.1f < R <= %.1f HLR' % (RDown, RUp)
                    
            c = RColor[iR] 
            x = np.ma.log10(aSFRSD__Trg[iT, RMask, :].flatten() * 1e6)
            y = np.ma.log10(aSFRSD_Ha__rg[RMask, :].flatten() * 1e6)
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            n_mask += mask.sum()
            n_tot += len(x)
  
            if i == 0 and j == 0:
                if iR > 0:
                    jump = 0.25 * (iR - 1.)
                pos_x = iR + jump
                #pos_y = pos_y_ini - (iR * pos_step)
                pos_y = 1.8
                textbox = dict(alpha = 0.)
                ax.text(pos_x, 1.1, txt,
                        fontsize = Rfontsize, color = c,
                        transform = ax.transAxes,
                        va = 'top', ha = 'left',
                        bbox = textbox)
              
            scat = ax.scatter(xm, ym, c = c, marker = 'o', s = 1., edgecolor = 'none', alpha = 1.)
          
        print 'SigmaSFR x SigmaSFR_Ha Age: %.2f Myr: masked %d points of %d (Total: %d)' % (age / 1e6, n_mask, n_tot, n_tot - n_mask)
        #ax.legend(loc = 'lower left', fontsize = 12, frameon = False)
        age = tSF__T[iT]
        xran = [-3.5, 1.]
        yran = [-3.5, 1.]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # binsx = np.linspace(-4.5, 1., 51)
        # binsy = np.linspace(min(ym),max(ym), 51)
        # density_contour(xm, ym, binsx, binsy, ax=ax)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
 
        x = np.ma.log10(aSFRSD__Trg[iT].flatten() * 1e6)
        y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
         
        a[iT], b[iT] = plotOLSbisectorAxis(ax, xm, ym, 0.92, 0.05, 8)

        b2[iT] = (ym - xm).mean()
        Y2 = xm + b2[iT]
        Yrms = (ym - Y2).std()
        ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
        txt = r'y = x + (%.2f) $Y_{rms}$:%.2f' %  (b2[iT], Yrms)
        plot_text_ax(ax, txt, 0.92, 0.15, 8, 'bottom', 'right', color = 'b')
        
        txt = '%.2f Myr' % (age / 1e6)
        plot_text_ax(ax, txt, 0.05, 0.92, 8, 'top', 'left')
         
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
              
        if i == NRows - 1 and j == 3:
            ax.set_xlabel(xlabel)
              
        if i == 2 and j == 0:
            ax.set_ylabel(ylabel)
  
        iT += 1
  
f.subplots_adjust(wspace=0, hspace=0, left=0.1, bottom=0.1, right=0.9, top=0.95)
f.savefig('aSFRSD_linregress_report.png')
plt.close(f)
 
x = np.log10(tSF__T)
xlabel = r'$\log\ t_\star$ [yr]'
 
plot_linreg_params(a, x, xlabel, 
                   'a', 'aSFRSD_linregress_slope_age.png', 1., 16) 
plot_linreg_params(b, x, xlabel, 
                   r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep_age.png', 0., 16)
plot_linreg_params(b2, x, xlabel, 
                   r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep2_age.png', 0., 16)
