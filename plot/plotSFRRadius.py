#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys
from plot_aux import get_attrib_h5, calcRunningStats, \
                     plotRunningStatsAxis
from scipy import stats as st


mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

h5 = h5py.File(sys.argv[1], 'r')

tSF__T = get_attrib_h5(h5, 'tSF__T')
RbinCenter__r = get_attrib_h5(h5, 'RbinCenter__r')
RRange = get_attrib_h5(h5, 'RRange')
RColor = get_attrib_h5(h5, 'RColor')
ALL_aSFRSD__Trg = get_attrib_h5(h5, 'ALL_aSFRSD__Trg')
ALL_aSFRSD_Ha__rg = get_attrib_h5(h5, 'ALL_aSFRSD_Ha__rg')
ALL_tauV__Trg = get_attrib_h5(h5, 'ALL_tauV__Trg')
ALL_tau_V_neb__rg = get_attrib_h5(h5, 'ALL_tau_V_neb__rg')
ALL_McorSD_GAL__rg = get_attrib_h5(h5, 'ALL_McorSD_GAL__rg')

h5.close()

for iT,tSF in enumerate(tSF__T):
   ###################### RADIUS COLOR ######################
   xlabel = r'$\log\ \Sigma_{SFR}^\star(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
   ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
   fname = 'alogSFRSD_alogSFRSD_neb_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
    
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
               
       x = np.ma.log10(ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten())
       mask = ~(x.mask | y.mask)
       xm = x[mask]
       ym = y[mask]
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])
    
   ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")

   x = np.ma.log10(ALL_aSFRSD__Trg[iT, :, :].flatten())
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten())
   mask = ~(x.mask | y.mask)
   xm = x[mask]
   ym = y[mask]
    
   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % ((ym/xm).mean(), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)

   xlabel = r'$\tau_V^\star(R)$'
   ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
   fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 

       x = ALL_tauV__Trg[iT, RMask, :].flatten()
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       mask = ~(x.mask | y.mask)
       xm = x[mask]
       ym = y[mask]
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])

   x = ALL_tauV__Trg[iT, :, :].flatten()
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
   mask = ~(x.mask | y.mask)
   xm = x[mask]
   ym = y[mask]

   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)
        
   xlabel = r'$\tau_V^{neb}(R)$'
   ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
   fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 
        
       x = ALL_tau_V_neb__rg[RMask, :].flatten()
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       mask = ~(x.mask | y.mask)
       xm = x[mask]
       ym = y[mask]
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])
        
   x = ALL_tau_V_neb__rg.flatten()
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
   mask = ~(x.mask | y.mask)
   xm = x[mask]
   ym = y[mask]
   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)
        
   xlabel = r'$\log\ \tau_V^\star(R)$'
   ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
   fname = 'tauV_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 

       x = np.ma.log10(ALL_tauV__Trg[iT, RMask, :].flatten())
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       mask = ~(x.mask | y.mask)
       xm = x[mask]
       ym = y[mask]
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])

   x = np.ma.log10(ALL_tauV__Trg[iT, :, :].flatten())
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
   mask = ~(x.mask | y.mask)
   xm = x[mask]
   ym = y[mask]
    
   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)
        
   xlabel = r'$\log\ \tau_V^{neb}(R)$'
   ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
   fname = 'tauVneb_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 

       x = np.ma.log10(ALL_tau_V_neb__rg[RMask, :].flatten())
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       mask = ~(x.mask | y.mask)
       xm = x[mask]
       ym = y[mask]
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])

   x = np.ma.log10(ALL_tau_V_neb__rg.flatten())
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
   mask = ~(x.mask | y.mask)
   xm = x[mask]
   ym = y[mask]
    
   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)
    
   xlabel = r'$\log\ \mathcal{M}_\star(R)$'
   ylabel = r'$\log\ (\Sigma_{SFR}^{neb}(R)/\Sigma_{SFR}^\star(R))$'
   fname = 'logMcorSD_SFRSDHa_SFRSD_age_%sMyr_RColor.png' % str(tSF / 1.e6)
   f = plt.figure(dpi = 96)
   f.set_size_inches(21.88,12.5)
   ax = f.gca()
   for iR, RUp in enumerate(RRange):
       if iR == 0:
           RMask = RbinCenter__r <= RUp
       else:
           RDown =  RRange[iR - 1]
           RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp) 

       x = np.ma.log10(ALL_McorSD_GAL__rg[RMask, :].flatten())
       y = np.ma.log10(ALL_aSFRSD_Ha__rg[RMask, :].flatten()/ALL_aSFRSD__Trg[iT, RMask, :].flatten())
       x.mask = y.mask
       xm = x
       ym = y
       scat = ax.scatter(xm, ym, c = RColor[iR], edgecolor = 'none', alpha = 0.8)
       plotRunningStatsAxis(ax, xm, ym, 'mean', RColor[iR])

   x = np.ma.log10(ALL_McorSD_GAL__rg.flatten())
   y = np.ma.log10(ALL_aSFRSD_Ha__rg.flatten()/ALL_aSFRSD__Trg[iT, :, :].flatten())
   x.mask = y.mask
   xm = x
   ym = y

   rhoSpearman, pvalSpearman = st.spearmanr(xm, ym)
   txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(ym/xm), np.ma.median((ym/xm)), np.ma.std(ym/xm), rhoSpearman)
   textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
   ax.text(0.3, 0.97, txt, fontsize = 28, transform = ax.transAxes, verticalalignment = 'top', bbox = textbox)
   ax.grid()
   ax.set_xlabel(xlabel)
   ax.set_ylabel(ylabel)
   ax.set_title(r'$%s$ Myr' % str(tSF / 1.e6))
   f.savefig(fname)
   plt.close(f)
