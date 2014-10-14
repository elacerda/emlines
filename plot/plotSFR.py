#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5
from scipy import stats as st


mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'


def plotSFR(x,y,xlabel,ylabel,xlim,ylim,age,fname):
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    yxlabel = r'$%s /\ %s $ ' % (ylabel.split('[')[0].strip('$ '), xlabel.split('[')[0].strip('$ '))
    txt = '%s mean:%.3f  median:%.3f  $\sigma(y/x)$:%.3f  Rs: %.2f' % (yxlabel, (y/x).mean(), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.03, 0.97, txt, fontsize = 16, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    if fname:
        f.savefig(fname)
    else:
        f.show()
    plt.close(f)


def plotTau(x,y,xlabel,ylabel,xlim,ylim,age,fname):
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    #ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    yxlabel = r'$%s /\ %s $ ' % (ylabel.split('[')[0].strip('$ '), xlabel.split('[')[0].strip('$ '))
    txt = '%s mean:%.3f  median:%.3f  $\sigma(y/x)$:%.3f  Rs: %.2f' % (yxlabel, (y/x).mean(), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(0.3, 0.97, txt, fontsize = 15, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    if fname:
        f.savefig(fname)
    else:
        f.show()
    plt.close(f)
    

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')

ALL_SFR__Tg = get_attrib_h5(h5, 'ALL_SFR__Tg')
ALL_SFR_Ha__g = get_attrib_h5(h5, 'ALL_SFR_Ha__g')
ALL_aSFRSD_kpc__Trg = get_attrib_h5(h5, 'ALL_aSFRSD_kpc__Trg')
ALL_aSFRSD_Ha_kpc__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha_kpc__rg')
ALL_tauV__Trg = get_attrib_h5(h5, 'ALL_tauV__Trg')
ALL_tau_V_neb__rg = get_attrib_h5(h5, 'ALL_tau_V_neb__rg')

h5.close()

for iT,tSF in enumerate(tSF__T):
    x = np.ma.log10(ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg[:, :].flatten())
    mask = ~(x.mask | y.mask)
    xm = x[mask]
    ym = y[mask]
    xlim = np.percentile(xm, [1, 100 * (xm.shape[0] - xm.mask.sum()) / xm.shape[0] - 1])
    ylim = np.percentile(ym, [1, 100 * (ym.shape[0] - ym.mask.sum()) / ym.shape[0] - 1])
    xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr.png' % str(tSF / 1.e6)
    plotSFR(xm,ym,xlabel,ylabel,xlim,ylim,tSF,fname)
    
    x = np.ma.log10(ALL_SFR__Tg[iT])
    y = np.ma.log10(ALL_SFR_Ha__g)
    xlabel = r'$\log\ SFR_\star\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$' 
    fname = 'logSFR_logSFR_neb_age_%sMyr.png' % str(tSF / 1.e6)
    xlim = np.percentile(x, [1, 100 * (x.shape[0] - x.mask.sum()) / x.shape[0] - 1])
    ylim = np.percentile(y, [1, 100 * (y.shape[0] - y.mask.sum()) / y.shape[0] - 1])
    plotSFR(x,y,xlabel,ylabel,xlim,ylim,tSF,fname)
        
    x = ALL_tauV__Trg[iT, :, :].flatten()
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    xlabel = r'$\tau_V^\star(R)$'
    ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    plotTau(x,y,xlabel,ylabel,None,None,tSF,fname) 
    
    x = ALL_tau_V_neb__rg.flatten()
    y = np.ma.log10(ALL_aSFRSD_Ha_kpc__rg.flatten() / ALL_aSFRSD_kpc__Trg[iT, :, :].flatten())
    xlabel = r'$\tau_V^{neb}(R)$'
    ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
    fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr.png' % str(tSF / 1.e6)
    plotTau(x,y,xlabel,ylabel,None,None,tSF,fname)
