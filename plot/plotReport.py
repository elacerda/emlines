#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5, density_contour
from matplotlib.ticker import MultipleLocator


mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

tSF_to_plot = [0, 10, 14, 17, 20, 23, 26, 29, 32, 35, 39 ]

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')
RbinCenter__r = get_attrib_h5(h5, 'RbinCenter__r')
RRange = get_attrib_h5(h5, 'RRange')
RColor = get_attrib_h5(h5, 'RColor')
correl_SFR__T = get_attrib_h5(h5, 'correl_SFR')
correl_SFRSD__rT = get_attrib_h5(h5, 'correl_SFRSD__rT')
correl_SFRSD_kpc__rT = get_attrib_h5(h5, 'correl_SFRSD_kpc__rT')
ALL_SFR__Tg = get_attrib_h5(h5, 'ALL_SFR__Tg')
ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
ALL_aSFRSD_kpc__Trg = get_attrib_h5(h5, 'ALL_aSFRSD_kpc__Trg')
ALL_aSFRSD_Ha_kpc__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha_kpc__rg')

h5.close()

NRows = 4
NCols = 3
  
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(10, 12)
   
for ax in f.axes:
    ax.set_axis_off()

xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'

NAxes = len(f.axes) 
       
k = 0
  
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j]

        if k < len(tSF_to_plot):
            iT = tSF_to_plot[k]
          
        if i < (NRows - 1) or j < (NCols - 1):
            ax.set_axis_on()
            x = np.ma.log10(ALL_SFR__Tg[iT])
            y = np.ma.log10(ALL_SFR_Ha__g)
            mask = x.mask | y.mask
            xm = x[~mask]
            ym = y[~mask]
            age = tSF__T[iT]
            print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, mask.sum(), len(x), len(x) - mask.sum())
            xran = [-6, 0]
            yran = [-6, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.3)
            binsx = np.linspace(-6., 0., 101)
            binsy = np.linspace(min(ym), max(ym), 101)
            density_contour(xm, ym, binsx, binsy, ax = ax)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            txt = 'age: %.2d Myr - $R_S$: %.2f' % ((age / 1e6), correl_SFR__T[iT])
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.05, 0.92, txt, fontsize = 12,
                    transform = ax.transAxes,
                    verticalalignment = 'top', horizontalalignment = 'left',
                    bbox = textbox)
            ax.grid()
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            if j == 0 and i == 1:
                ax.set_ylabel(ylabel)
            if j == 1 and i == (NRows - 1):
                ax.set_xlabel(xlabel)
            k += 1
        elif i == (NRows - 1) and j == (NCols - 1):
            ax.set_axis_on()
            ax.plot(np.log10(tSF__T), correl_SFR__T, 'k-', label = r'$R_S$')
            ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
            ax.legend(fontsize = 12, frameon = False)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.set_ylim([0., 1.])
            ax.grid(which = 'minor')

f.savefig('SFR_report.png')
plt.close(f)

NCols = 3 
NRows = 4

pos_y_ini = 0.38
pos_step = 0.09
Rfontsize = 10
  
f, axArr = plt.subplots(NRows, NCols)
f.set_dpi(96)
f.set_size_inches(10, 12)
   
for ax in f.axes:
    ax.set_axis_off()

xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
   
NAxes = len(f.axes) 
k = 0
       
for i in range(0, NRows):
    for j in range(0, NCols):
        ax = axArr[i, j]

        if k < len(tSF_to_plot):
            iT = tSF_to_plot[k]
         
        if i < (NRows - 1) or j < (NCols - 1):
            ax.set_axis_on()
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
                x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, RMask, :].flatten())
                y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[RMask, :].flatten())
                mask = x.mask | y.mask
                xm = x[~mask]
                ym = y[~mask]
                n_mask += mask.sum()
                n_tot += len(x)
                 
                if i == 0 and j == 0:
                    pos_y = pos_y_ini - (iR * pos_step)
                    textbox = dict(alpha = 0.)
                    ax.text(0.05, pos_y, txt,
                            fontsize = Rfontsize, color = c,
                            transform = ax.transAxes,
                            va = 'top', ha = 'left',
                            bbox = textbox)
                scat = ax.scatter(xm, ym, c = c, marker = 'o', s = 1., edgecolor = 'none', alpha = 1.)
            print 'SigmaSFR x SigmaSFR_Ha Age: %.2f Myr: masked %d points of %d (Total: %d)' % (age / 1e6, n_mask, n_tot, n_tot - n_mask)
                     
            ax.legend(loc = 'lower left', fontsize = 12, frameon = False)
            age = tSF__T[iT]
            xran = [-3.5, 1.]
            yran = [-3.5, 1.]
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # binsx = np.linspace(-4.5, 1., 51)
            # binsy = np.linspace(min(ym),max(ym), 51)
            # density_contour(xm, ym, binsx, binsy, ax=ax)
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            txt = 'age: %.2d Myr - $R_S$: %.2f' % ((age / 1e6), correl_SFRSD__rT[0, iT])
            textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
            ax.text(0.05, 0.92, txt, fontsize = 12,
                    transform = ax.transAxes,
                    verticalalignment = 'top', horizontalalignment = 'left',
                    bbox = textbox)
            ax.grid()
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
         
            if j == 0 and i == 1:
                ax.set_ylabel(ylabel)
            if j == 1 and i == (NRows - 1):
                ax.set_xlabel(xlabel)
            k += 1
        elif i == (NRows - 1) and j == (NCols - 1):
            ax.set_axis_on()
            ax.plot(np.log10(tSF__T), correl_SFRSD__rT[0, :], 'k-', label = r'$R_S$')

            for iR, RUp in enumerate(RRange):
                if iR == 0:
                    RMask = RbinCenter__r <= RUp
                else:
                    RDown = RRange[iR - 1]
                    RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
                iiR = iR + 1
                 
                c = RColor[iR]
                ax.plot(np.log10(tSF__T), correl_SFRSD__rT[iiR, :], color = c, ls = '--', label = None)

            ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
            ax.legend(fontsize = 12, frameon = False)
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))
            ax.set_ylim([0., 1.])
            ax.grid(which = 'minor')

f.savefig('SFRSD_report.pdf')
plt.close(f)
  
f = plt.figure()
ax = f.gca()
ax.plot(np.log10(tSF__T), correl_SFRSD_kpc__rT[0, :], 'k-', label = r'$R_S(\Sigma_{\mathrm{SFR}})$')
ax.plot(np.log10(tSF__T), correl_SFR__T, 'k--', label = r'$R_S(\mathrm{SFR})$')
ax.set_xlabel(r'$\log\ t_\star\ [yr]$')
ax.set_ylabel(r'$R_S$')
ax.set_ylim([0, 1.])
ax.grid()
ax.legend(loc = 'best', fontsize = 12)
f.tight_layout()
f.savefig('Rs_SFR.pdf')
plt.close(f)
