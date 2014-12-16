#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import h5py
from scipy import stats as st
import scipy.optimize as so
from matplotlib import pyplot as plt


def plot_linreg_params(param, x, xlabel, ylabel, fname, best_param = None, fontsize = 12):
    y = param
    xm = x
    ym = y
    #xm = x[~y.mask]
    #ym = y[~y.mask]
    f = plt.figure()
    ax = f.gca()
    ax.scatter(xm, ym, c = 'k', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.6)
    if best_param != None:
        ax.axhline(y = best_param, c = 'k', ls = '--')
        delta = np.abs(ym - best_param)
        where = np.where(delta == delta.min())
        txt = r'$x_{best}:$ %.2f ($\Delta$:%.2f)' % (x[where], delta[where])
        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
        x_pos = xm[where]
        y_pos = ym[where]
        xlim_inf, xlim_sup = ax.get_xlim()
        ylim_inf, ylim_sup = ax.get_ylim() 
        arrow_size_x = (xlim_sup - x_pos) / 6
        arrow_size_y = (ylim_sup - y_pos) / 3 
        ax.annotate(txt,
            xy=(x_pos, y_pos), xycoords='data',
            xytext=(x_pos + arrow_size_x, y_pos + arrow_size_y), 
            textcoords='data',
            verticalalignment = 'top', horizontalalignment = 'left',
            bbox = textbox,
            arrowprops=dict(arrowstyle="->", 
                            #linestyle="dashed",
                            color="0.5",
                            connectionstyle="angle3,angleA=90,angleB=0",
                            ),
            )
#        ax.text(xm[where], ym[where], txt, fontsize = fontsize,
#                verticalalignment = 'top', horizontalalignment = 'left',
#                bbox = textbox)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    f.savefig(fname)


def plot_text_ax(ax, txt, xpos, ypos, fontsize, va, ha, color = 'k'):
    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
    ax.text(xpos, ypos, txt, fontsize = fontsize, color = color,
            transform = ax.transAxes,
            verticalalignment = va, horizontalalignment = ha,
            bbox = textbox)


def plotRunningStatsAxis(ax, x, y, ylegend, plot_stats = 'mean', color = 'black', errorbar = True, nBox = 25):
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    
    if plot_stats == 'median':
        xx = xMedian
        yy = yMedian
    else:
        xx = xMean
        yy = yMean

    ax.plot(xx, yy, color, lw = 2, label = ylegend)
    
    if errorbar:
        ax.errorbar(xx, yy, yerr = yStd, xerr = xStd, c = color)


def get_attrib_h5(h5, attrib):
    if any([ attrib in s for s in h5['masked/mask'].keys() ]):
        node = '/masked/data/' + attrib
        ds = h5[node]
        if type(ds) == h5py.Dataset:
            data = h5.get('/masked/data/' + attrib).value
            mask = h5.get('/masked/mask/' + attrib).value
            arr = np.ma.masked_array(data, mask = mask)
        else:
            arr = []
            tSF__T = h5.get('/data/tSF__T').value
            
            for iT, tSF in enumerate(tSF__T):
                group = '%s/%d' % (attrib, iT)
                data = h5.get('/masked/data/' + group).value
                mask = h5.get('/masked/mask/' + group).value
                arr.append(np.ma.masked_array(data, mask = mask))
        return arr 
    else:
        return h5.get('/data/' + attrib).value


def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

 
def density_contour(xdata, ydata, binsx, binsy, ax = None, **contour_kwargs):
    """ Create a density contour plot.
 
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """    
    #nbins_x = len(binsx) - 1
    #nbins_y = len(binsy) - 1

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins = [binsx, binsy], normed = True)
    x_bin_sizes = (xedges[1:] - xedges[:-1])
    y_bin_sizes = (yedges[1:] - yedges[:-1])
 
    pdf = (H * (x_bin_sizes * y_bin_sizes))
 
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args = (pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args = (pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args = (pdf, 0.99))
    levels = [one_sigma, two_sigma, three_sigma]
 
    X, Y = 0.5 * (xedges[1:] + xedges[:-1]), 0.5 * (yedges[1:] + yedges[:-1])
    Z = pdf.T
 
    if ax == None:
        contour = plt.contour(X, Y, Z, levels = levels, origin = "lower", **contour_kwargs)
    else:
        contour = ax.contour(X, Y, Z, levels = levels, origin = "lower", **contour_kwargs)
 
    return contour


def calcRunningStats(x, y, dxBox = 0.3, xbinIni = 8.5, xbinFin = 12, xbinStep = 0.05):
    '''
    Statistics of x & y in equal size x-bins (dx-box).
    Note the mery small default xbinStep, so we have overlaping boxes.. so running stats..

    Cid@Lagoa -
    '''

    # Def x-bins
    xbin = np.arange(xbinIni, xbinFin + xbinStep, xbinStep)
    xbinCenter = (xbin[:-1] + xbin[1:]) / 2.0
    Nbins = len(xbinCenter)

    # Reset in-bin stats arrays
    xMedian , xMean , xStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    yMedian , yMean , yStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    nInBin = np.zeros(Nbins)

    # fill up in x & y stats for each x-bin
    for ixBin in np.arange(Nbins):
        isInBin = (np.abs(x - xbinCenter[ixBin]) <= dxBox / 2.)
        xx , yy = x[isInBin] , y[isInBin]
        xMedian[ixBin] , xMean[ixBin] , xStd[ixBin] = np.median(xx) , xx.mean() , xx.std()
        yMedian[ixBin] , yMean[ixBin] , yStd[ixBin] = np.median(yy) , yy.mean() , yy.std()
        nInBin[ixBin] = isInBin.sum()

    return xbinCenter, xMedian, xMean, xStd, yMedian, yMean, yStd, nInBin


def plotStatCorreAxis(ax, x, y, pos_x, pos_y, fontsize):
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = '<y/x>:%.3f - (y/x) median:%.3f - $\sigma(y/x)$:%.3f - Rs: %.2f' % (np.mean(y/x), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    plot_text_ax(ax, txt, pos_x, pos_y, fontsize, 'top', 'left')


def gaussSmooth_YofX(x, y, FWHM):
    '''
    Sloppy function to return the gaussian-smoothed version of an y(x) relation.
    Cid@Lagoa - 07/June/2014
    '''

    sig = FWHM / np.sqrt(8. * np.log(2))
    xS , yS = np.zeros_like(x), np.zeros_like(x)
    w__ij = np.zeros((len(x), len(x)))
    for i in np.arange(len(x)):
        # for j in np.arange(len(x)):
        #     w__ij[i,j] = np.exp( -0.5 * ((x[j] - x[i]) / sig)**2  )

        w__ij[i, :] = np.exp(-0.5 * ((x - x[i]) / sig) ** 2)
        w__ij[i, :] = w__ij[i, :] / w__ij[i, :].sum()

        xS[i] = (w__ij[i, :] * x).sum()
        yS[i] = (w__ij[i, :] * y).sum()

    return xS , yS


def plotSFR(x,y,xlabel,ylabel,xlim,ylim,age,fname):
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    scat = ax.scatter(x, y, c = 'black', edgecolor = 'none', alpha = 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), ls="--", c=".3")
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    yxlabel = r'$%s /\ %s $ ' % (ylabel.split('[')[0].strip('$ '), xlabel.split('[')[0].strip('$ '))
    txt = '%s mean:%.3f  median:%.3f  $\sigma(y/x)$:%.3f  Rs: %.2f' % (yxlabel, (y/x).mean(), np.ma.median((y/x)), np.ma.std(y/x), rhoSpearman)
    plot_text_ax(ax, txt, 0.03, 0.97, 16, 'top', 'left')
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
    plot_text_ax(ax, txt, 0.03, 0.97, 15, 'top', 'left')
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
    
    
def plotScatterColor(x, y, z, xlabel, ylabel, zlabel, xlim, ylim, zlim, age, fname):
    f = plt.figure()
    f.set_size_inches(10,9)
    ax = f.gca()
    if zlim != None:
        sc = ax.scatter(x, y, c = z, cmap = 'spectral_r', vmin = zlim[0], vmax = zlim[1], marker = 'o', s = 5., edgecolor = 'none')
    else:
        sc = ax.scatter(x, y, c = z, cmap = 'spectral_r', marker = 'o', s = 5., edgecolor = 'none')
    ax.plot((y - np.log10(1.6e-4))/1.4, y, 'k--')
    binsx = np.linspace(min(x), max(x), 21)
    binsy = np.linspace(min(y), max(y), 21)
    density_contour(x, y, binsx, binsy, ax = ax, color = 'k')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim != None:
        ax.set_xlim(xlim)
    if ylim != None:
        ax.set_ylim(ylim)
    plotStatCorreAxis(ax, x, y, 0.03, 0.97, 16)
    #plotRunningStatsAxis(ax, x, y, 'k')    
    nBox = 20
    dxBox       = (x.max() - x.min()) / (nBox - 1.)
    aux         = calcRunningStats(x, y, dxBox = dxBox, xbinIni = x.min(), xbinFin = x.max(), xbinStep = dxBox)
    xbinCenter  = aux[0]
    xMedian     = aux[1]
    xMean       = aux[2]
    xStd        = aux[3]
    yMedian     = aux[4]
    yMean       = aux[5]
    yStd        = aux[6]
    nInBin      = aux[7]
    ax.plot(xMedian, yMedian, 'k', lw = 2)
    cb = f.colorbar(sc)
    cb.set_label(zlabel)
    ax.set_title(r'$%s$ Myr' % str(age / 1.e6))
    f.savefig(fname)
    plt.close(f)
    

def calcYofXStats_EqNumberBins(x, y, nPerBin = 25):
    '''
    This gives the statistics of y(x) for x-bins of variable width, but all containing
    the same number of points.
    We 1st sort x, and the y accordingly. Then we compute the median, mean and std
    of x & y in contiguous x-bins in x defined to have nPerBin points each

    Cid@Lagoa - 05/June/2014
    '''

    ind_sx = np.argsort(x)
    xS , yS = x[ind_sx] , y[ind_sx]

    Nbins = len(x) - nPerBin + 1
    xMedian , xMean , xStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    yMedian , yMean , yStd = np.zeros(Nbins) , np.zeros(Nbins) , np.zeros(Nbins)
    nInBin = np.zeros(Nbins)

    for ixBin in np.arange(0, Nbins):
        xx , yy = xS[ixBin:ixBin + nPerBin] , yS[ixBin:ixBin + nPerBin]
        xMedian[ixBin] , xMean[ixBin] , xStd[ixBin] = np.median(xx) , xx.mean() , xx.std()
        yMedian[ixBin] , yMean[ixBin] , yStd[ixBin] = np.median(yy) , yy.mean() , yy.std()
        nInBin[ixBin] = len(xx)
    return xMedian, xMean, xStd, yMedian, yMean , yStd, nInBin


def data_uniq(list_gal, data):
    list_uniq_gal = np.unique(list_gal)
    NGal = len(list_uniq_gal)
    data__g = np.ones((NGal))
    
    for i, g in enumerate(list_uniq_gal):
        data__g[i] = np.unique(data[np.where(list_gal == g)])
        
    return NGal, list_uniq_gal, data__g


def list_gal_sorted_by_data(list_gal, data, type):
    '''
    type = 0 - sort asc
    type = 1 - sort desc
    type = -1 - receives list_gal and data as uniq
                aka. only sorts, list_gal by data
    
    '''
    if type >= 0:
        NGal, list_uniq_gal, data__g = data_uniq(list_gal, data)
    else:
        NGal = len(list_gal)
        list_uniq_gal = list_gal
        data__g = data
        
    iS = np.argsort(data__g)
    
    if type != 0:
        iS = iS[::-1]
    
    return list_uniq_gal[iS]
