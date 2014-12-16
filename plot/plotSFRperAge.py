#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Lacerda@Granada - 26/Nov/2014
#
import numpy as np
#import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import H5SFRData, plot_text_ax
#from matplotlib.ticker import MultipleLocator
from scipy import stats as st
from astropy.stats import sigma_clip

try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

K = H5SFRData(h5file)

min_pixel_to_plot = 5
tSF__T = K.get_data_h5('tSF__T')

xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'

###################################################################################
for iT, age in enumerate(tSF__T):
    x = np.ma.log10(K.get_data_h5('SFR__Tg')[iT])
    y = np.ma.log10(K.get_data_h5('SFR_Ha__g'))

    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]

    f = plt.figure()
    f.set_dpi(100)
    f.set_size_inches(11.69,8.27) 
    plot_suptitle = '%.2f Myr' % (age/1e6)
    f.suptitle(plot_suptitle)
    ax = f.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    xran = [-6, 0]
    yran = [-6, 0]
    scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 5, edgecolor = 'none', alpha = 0.5, label='')
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3", label = '')
      
    step = (xm.max() - xm.min()) / len(xm)
    X = np.linspace(xm.min(), xm.max() + step, len(xm))
    A, B, Rp, pval, std_err = st.linregress(xm, ym)
    ax.plot(X, A * X + B, c = 'k', ls = '--', lw = 2)
    txt = '%.2f Myr' % (age / 1e6)
    plot_text_ax(ax, txt, 0.05, 0.92, 14, 'top', 'left')
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20, 14, 'bottom', 'right')
    
    csig = 1
    c = 'r'
    xclip = sigma_clip(xm, csig)
    yclip = sigma_clip(ym, csig)
    maskclip = ~(xclip.mask | yclip.mask)
    xclipm = xclip[maskclip]
    yclipm = yclip[maskclip]
    A, B, Rp, pval, std_err = st.linregress(xclipm, yclipm)
    X = np.linspace(xm.min(), xm.max(), len(xm))
    ax.plot(X, A * X + B, c = c, ls = '--', lw = 2, label = r'$%d\sigma$' % csig)
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20 - (csig * 0.07), 14, 'bottom', 'right', c)

    csig = 2
    c = 'b'
    xclip = sigma_clip(xm, csig)
    yclip = sigma_clip(ym, csig)
    maskclip = ~(xclip.mask | yclip.mask)
    xclipm = xclip[maskclip]
    yclipm = yclip[maskclip]
    A, B, Rp, pval, std_err = st.linregress(xclipm, yclipm)
    X = np.linspace(xm.min(), xm.max(), len(xm))
    ax.plot(X, A * X + B, c = c, ls = '--', lw = 2, label = r'$%d\sigma$' % csig)
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20 - (csig * 0.07), 14, 'bottom', 'right', c)

    ax.set_xlim(xran)
    ax.set_ylim(yran)
            
    ax.legend()
    f.savefig('SFR_SFRHa_%.2fMyr.png' % (age / 1e6))


xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 

###################################################################################
for iT, age in enumerate(tSF__T):
    x = np.ma.log10(K.get_data_h5('aSFRSD_kpc__Trg')[iT].flatten())
    y = np.ma.log10(K.get_data_h5('aSFRSD_Ha_kpc__rg').flatten())

    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]

    f = plt.figure()
    f.set_dpi(100)
    f.set_size_inches(11.69,8.27) 
    plot_suptitle = '%.2f Myr' % (age/1e6)
    f.suptitle(plot_suptitle)
    ax = f.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    xran = [-3.5, 1.]
    yran = [-3.5, 1.]
    scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 5, edgecolor = 'none', alpha = 0.5, label='')
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3", label = '')
      
    step = (xm.max() - xm.min()) / len(xm)
    X = np.linspace(xm.min(), xm.max() + step, len(xm))
    A, B, Rp, pval, std_err = st.linregress(xm, ym)
    ax.plot(X, A * X + B, c = 'k', ls = '--', lw = 2)
    txt = '%.2f Myr' % (age / 1e6)
    plot_text_ax(ax, txt, 0.05, 0.92, 14, 'top', 'left')
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20, 14, 'bottom', 'right')
    
    csig = 1
    c = 'r'
    xclip = sigma_clip(xm, csig)
    yclip = sigma_clip(ym, csig)
    maskclip = ~(xclip.mask | yclip.mask)
    xclipm = xclip[maskclip]
    yclipm = yclip[maskclip]
    A, B, Rp, pval, std_err = st.linregress(xclipm, yclipm)
    X = np.linspace(xm.min(), xm.max(), len(xm))
    ax.plot(X, A * X + B, c = c, ls = '--', lw = 2, label = r'$%d\sigma$' % csig)
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20 - (csig * 0.07), 14, 'bottom', 'right', c)

    csig = 2
    c = 'b'
    xclip = sigma_clip(xm, csig)
    yclip = sigma_clip(ym, csig)
    maskclip = ~(xclip.mask | yclip.mask)
    xclipm = xclip[maskclip]
    yclipm = yclip[maskclip]
    A, B, Rp, pval, std_err = st.linregress(xclipm, yclipm)
    X = np.linspace(xm.min(), xm.max(), len(xm))
    ax.plot(X, A * X + B, c = c, ls = '--', lw = 2, label = r'$%d\sigma$' % csig)
    txt = r'$\log\ SFR_{neb} = %.2f\ \log\ \overline{SFR_\star}(t_\star)\ +\ %.2f\ (\sigma:%.4f)$' %  (A, B, std_err)
    plot_text_ax(ax, txt, 0.98, 0.20 - (csig * 0.07), 14, 'bottom', 'right', c)

    ax.set_xlim(xran)
    ax.set_ylim(yran)
    ax.legend()
    f.savefig('aSFRSD_aSFRSDHa_%.2fMyr.png' % (age / 1e6))
