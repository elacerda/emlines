#!/usr/bin/python
#
# Lacerda@Granada - 14/Apr/2015
#
import sys
import numpy as np
import CALIFAUtils as C
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st

#debug = True
debug = False
#mask_radius = True
mask_radius = False

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

Rnuc = 0.5

def corr_xy_SFR(H, v, tSF__T, iT, mask_radius = False):
    C.debug_var(debug, corr = v)
    C.debug_var(debug, iT = iT)
    C.debug_var(debug, tSF = tSF__T[iT])
    if v == 'aSFRSD':
        dim = 'rg'
    else:
        dim = 'g'
    x = H.get_data_h5('%s__T%s' % (v, dim))[iT]
    y = H.get_data_h5('%s_Ha__%s' % (v, dim))
    xm, ym = C.ma_mask_xyz(x, y)
    if mask_radius is True:
        m = (H.zone_dist_HLR__g > 0.5) 
        xm[~m] = np.ma.masked
        ym[~m] = np.ma.masked 
    Rs, _ = st.spearmanr(xm.compressed(), ym.compressed())
    Rp, _ = st.pearsonr(xm.compressed(), ym.compressed())
    C.debug_var(debug, Rs = Rs)
    C.debug_var(debug, Rp = Rp)
    return Rs, Rp

if __name__ == '__main__':
    try:
        h5file = sys.argv[1]
    except IndexError:
        print 'usage: %s HDF5FILE' % (sys.argv[0])
        exit(1)
    
    H = C.H5SFRData(h5file)    
    
    tSF__T = H.tSF__T
    C.debug_var(debug, tSF__T = tSF__T)
    iT__t = np.arange(len(tSF__T))
    C.debug_var(debug, iTs = iT__t)
    RsSFR = np.empty(iT__t.shape, dtype = np.double)
    RpSFR = np.empty(iT__t.shape, dtype = np.double)
    RsSFRSD = np.empty(iT__t.shape, dtype = np.double)
    RpSFRSD = np.empty(iT__t.shape, dtype = np.double)
    RsaSFRSD = np.empty(iT__t.shape, dtype = np.double)
    RpaSFRSD = np.empty(iT__t.shape, dtype = np.double)
    RsaSFRSDn = np.empty(iT__t.shape, dtype = np.double)
    RpaSFRSDn = np.empty(iT__t.shape, dtype = np.double)
    RsiSFR = np.empty(iT__t.shape, dtype = np.double)
    RpiSFR = np.empty(iT__t.shape, dtype = np.double)
    RsiSFRSD = np.empty(iT__t.shape, dtype = np.double)
    RpiSFRSD = np.empty(iT__t.shape, dtype = np.double)

    for iT in iT__t:
        tSF = tSF__T[iT]
        corr = ['SFR', 'SFRSD', 'aSFRSD', 'integrated_SFR', 'integrated_SFRSD' ]
        R = [
            [RsSFR, RpSFR], [RsSFRSD, RpSFRSD], [RsaSFRSD, RpaSFRSD],
            [RsiSFR, RpiSFR], [RsiSFRSD, RpiSFRSD]
        ]
        for i, c in enumerate(corr):
            R[i][0][iT], R[i][1][iT] = corr_xy_SFR(H, c, tSF__T, iT, mask_radius)
        SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][10, :] + H.aSFRSD__Trg[iT][9, :] / 2.)
        aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
        aSFRSD_Ha_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
        C.debug_var(debug, iT = iT)
        C.debug_var(debug, tSF = tSF__T[iT])
        xm, ym = C.ma_mask_xyz(aSFRSD_norm__rg, aSFRSD_Ha_norm__rg)
        if mask_radius is True:
            m = (H.zone_dist_HLR__g > 0.5) 
            xm[~m] = np.ma.masked
            ym[~m] = np.ma.masked 
        Rs, _ = st.spearmanr(xm.compressed(), ym.compressed())
        Rp, _ = st.pearsonr(xm.compressed(), ym.compressed())
        C.debug_var(debug, Rs = Rs)
        C.debug_var(debug, Rp = Rp)
        RsaSFRSDn[iT] = Rs
        RpaSFRSDn[iT] = Rp
    c_label = ['SFR', r'$\Sigma_{SFR}$', r'$\Sigma_{SFR}$(R)', r'$SFR^{int}$', r'$\Sigma_{SFR}^{int}$' ]
    f = plt.figure(figsize = (10, 8), dpi = 100)
    ax = f.gca()
    for i, c in enumerate(corr):
        ax.plot(H.tSF__T, R[i][0], label = c_label[i])
    ax.plot(H.tSF__T, RsaSFRSDn, label = r'$\frac{\Sigma_{SFR}(R)}{\Sigma_{SFR}(1HLR)}$')
    plt.xscale('log')
    plt.legend(loc ='upper right', fontsize = 10)
    plt.xlabel(r'$t_{SF}\ [yr]$')
    plt.ylabel(r'$R_s$')
    plt.grid(which = 'minor') 
    plt.grid(which = 'major') 
    #plt.title(r'Stellar vs $H\alpha$')
    if mask_radius is True:
        suf = '_maskradius.png'
    else:
        suf = '.png'
    suptitle = r'NGals:%d  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    f.suptitle(suptitle)
    f.savefig('SFR_correlations%s' % suf)
    plt.close(f)
    ###################
    ###################
    ###################
