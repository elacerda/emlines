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
from sklearn import linear_model

try:
    h5file = sys.argv[1]
    try:
        iT = np.int(sys.argv[2])
    except IndexError:
        iT = -1
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

def plotSFRage(x, y, xlabel, ylabel, xlim, ylim, fname):
    f = plt.figure()
    f.set_dpi(100)
    f.set_size_inches(11.69,8.27) 
    plot_suptitle = '%.2f Myr' % (age/1e6)
    f.suptitle(plot_suptitle)
    ax = f.gca()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    ax.scatter(xm, ym, c = 'black', marker = 'o', s = 5, edgecolor = 'none', alpha = 0.5, label='')
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3", label = '')
    
    c = 'b'
    step = (xm.max() - xm.min()) / len(xm)
    A1, B1, Rp, pval, std_err = st.linregress(xm, ym)
    X = np.linspace(xm.min(), xm.max() + step, len(xm))
    Y = A1 * X + B1
    Yrms = (Y - ym).std()
    ax.plot(X, Y, c = c, ls = '--', lw = 2, label = 'least squares')
    txt = '%.2f Myr' % (age / 1e6)
    plot_text_ax(ax, txt, 0.05, 0.92, 14, 'top', 'left')
    txt = r'$y = %.2f\ x\ +\ (%.2f)\ (y_{rms}:%.2f)$' %  (A1, B1, Yrms)
    plot_text_ax(ax, txt, 0.98, 0.21, 14, 'bottom', 'right', c)
    
    c = 'g'
    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
    model_ransac.fit(np.vstack(xm),np.vstack(ym))
    Y = model_ransac.predict(X[:, np.newaxis])
    Yrms = (Y - ym).std()
    ax.plot(X, Y, c = c, ls = '--', lw = 2, label = 'RANSAC')
    A = model_ransac.estimator_.coef_
    B = model_ransac.estimator_.intercept_
    inlier_mask = model_ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    txt = r'$y = %.2f\ x\ +\ (%.2f)\ (y_{rms}:%.2f)$' %  (A, B, Yrms)
    plot_text_ax(ax, txt, 0.98, 0.14, 14, 'bottom', 'right', c)
    ax.scatter(xm[inlier_mask], ym[inlier_mask], c = c, marker = 'x', s = 20, facecolor = 'k', edgecolor = c, alpha = 0.3, label='')
    
    c = 'r'
    A2, B2, Rp, pval, std_err = st.linregress(ym, xm)
    X = np.linspace(xm.min(), xm.max() + step, len(xm))
    A = ((A1 / A2 - 1. + ((1. + A1 ** 2.) * (1. + A2 ** -2.)) ** 0.5) / (A1 + (1. / A2)))
    B = ym.mean() - A * xm.mean()
    Y = A * X + B
    Yrms = (Y - ym).std()
    ax.plot(X, Y, c = c, ls = '--', lw = 2, label = 'OLS bisector')
    txt = r'$y = %.2f\ x\ +\ (%.2f)\ (y_{rms}:%.2f)$' %  (A, B, Yrms)
    plot_text_ax(ax, txt, 0.98, 0.07, 14, 'bottom', 'right')

    ax.legend()
    f.savefig(fname)
    plt.close(f)

K = H5SFRData(h5file)

min_pixel_to_plot = 5
###################################################################################
if iT < 0:
    for iT, age in enumerate(K.tSF__T):
        x = np.ma.log10(K.get_data_h5('SFR__Tg')[iT])
        y = np.ma.log10(K.get_data_h5('SFR_Ha__g'))
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        xran = [-5, 0]
        yran = [-5, 0]
        xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
        ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
        fname = 'SFR_SFRHa_%.2fMyr.png' % (age / 1e6)
        plotSFRage(xm, ym, xlabel, ylabel, xran, yran, fname)
    
        x = np.ma.log10(K.get_data_h5('aSFRSD_kpc__Trg')[iT].flatten())
        y = np.ma.log10(K.get_data_h5('aSFRSD_Ha_kpc__rg').flatten())
        mask = x.mask | y.mask
        xm = x[~mask]
        ym = y[~mask]
        xran = [-3.5, 1.]
        yran = [-3.5, 1.]
        xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
        ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
        fname = 'aSFRSD_aSFRSDHa_%.2fMyr.png' % (age / 1e6)
        plotSFRage(xm, ym, xlabel, ylabel, xran, yran, fname) 
else:
    age = K.tSF__T[iT]
    x = np.ma.log10(K.get_data_h5('SFR__Tg')[iT])
    y = np.ma.log10(K.get_data_h5('SFR_Ha__g'))
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    xran = [-5, 0]
    yran = [-5, 0]
    xlabel = r'$\log\ \overline{SFR_\star}(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    fname = 'SFR_SFRHa_%.2fMyr.png' % (age / 1e6)
    plotSFRage(xm, ym, xlabel, ylabel, xran, yran, fname)

    x = np.ma.log10(K.get_data_h5('aSFRSD_kpc__Trg')[iT].flatten())
    y = np.ma.log10(K.get_data_h5('aSFRSD_Ha_kpc__rg').flatten())
    mask = x.mask | y.mask
    xm = x[~mask]
    ym = y[~mask]
    xran = [-3.5, 1.]
    yran = [-3.5, 1.]
    xlabel = r'$\log\ \overline{\Sigma_{SFR}^\star}(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$'
    fname = 'aSFRSD_aSFRSDHa_%.2fMyr.png' % (age / 1e6)
    plotSFRage(xm, ym, xlabel, ylabel, xran, yran, fname) 