#!/usr/bin/python
#
# Lacerda@Granada - 04/Mar/2015
#
from CALIFAUtils.scripts import sort_gals, debug_var, loop_cubes
from CALIFAUtils.plots import plot_zbins
import numpy as np
import sys


def mask_line(EL, line, dw, minSNR, maxsigma):
    mask = np.zeros((EL.N_zone), dtype = np.bool)
    np.bitwise_or(EL._setMaskLineFluxNeg(line), EL._setMaskLineDisplacement(line, dw), mask)
    np.bitwise_or(mask, EL._setMaskLineSNR(line, minSNR), mask)
    np.bitwise_or(mask, EL._setMaskLineSigma(line, maxsigma), mask)
    return mask


if __name__ == '__main__':
    debug = True
    kwargs = dict()
    try:
        filename = sys.argv[1]
        try:
            imax = np.int(sys.argv[2])
            kwargs.update(dict(imax = imax))
        except:
            pass
    except:
        filename = '/Users/lacerda/CALIFA/listOf300GalPrefixes.txt'
    kwargs.update(dict(EL = True))    
    debug_var(debug, kwargs = kwargs)
    gals = sort_gals(filename, return_data_sort = False)
    lines = {
        '5007' : dict(v = [], name = '[OIII]5007', ccm = 1.12022),
        '4959' : dict(v = [], name = '[OIII]4959', ccm = 1.13427),
        '6583' : dict(v = [], name = '[NII]6583', ccm = 0.81466),
        '6548' : dict(v = [], name = '[OIII]6548', ccm = 0.82006),
    }
    plot_ratios = {
        'OIII' : dict(keys = [ '5007', '4959' ], ilines = [], gals = []),
        'NII' : dict(keys = [ '6583', '6548' ], ilines = [], gals = []),
    }
    zoneToGal = []
    g = 0
    for iGal, K in loop_cubes(gals, **kwargs):
        minSNR = 3
        dw = 3.
        maxsigma = 3.5
        maskHa = mask_line(K.EL, '6563', 5., minSNR, maxsigma)
        maskHb = mask_line(K.EL, '4861', 3., minSNR, maxsigma)
        maskTau = np.bitwise_or(maskHa, maskHb)
        for k, v in plot_ratios.iteritems():
            kk = v['keys']
            Nkk = len(kk)
            try:
                tmp = []
                for line in kk:
                    tmp.append(K.EL.lines.index(line))
            except:
                print 'miss line %s' % line 
                debug_var(debug, gal = gals[iGal])
                pass
            v['ilines'].append(tmp)
            v['gals'].append(gals[iGal])
            mask = np.bitwise_or(mask_line(K.EL, kk[0], dw, minSNR, maxsigma), 
                                 mask_line(K.EL, kk[1], dw, minSNR, maxsigma)) 
            for i in xrange(Nkk):
                ext_law = lines[kk[i]]['ccm']
                fc = K.EL._getLineFlux(kk[i])
                fi = K.EL.intrinsic_flux__z(line = fc, ext_law = ext_law)
                tmp = np.ma.masked_where(np.bitwise_or(mask, maskTau), fi, copy=True)
                lines[v['keys'][i]]['v'].append(tmp)
                if debug and i == 1:
                    y = lines[kk[0]]['v'][g] / lines[kk[1]]['v'][g]
                    if (y > 7).any():
                        print gals[iGal], y.argmax(), y.max()
        # this g index is only for debug purpose. 
        g = g + 1
    for k, v in plot_ratios.iteritems():
        for line in v['keys']:
            lines[line]['v'] = np.ma.hstack(lines[line]['v'])
        x = np.ma.log10(lines[v['keys'][0]]['v'])
        y = lines[v['keys'][0]]['v'] / lines[v['keys'][1]]['v']
        xname = lines[v['keys'][0]]['name']
        yname = lines[v['keys'][1]]['name']
        plot_zbins(
            debug = debug,
            x = x,
            xlabel = r'$\log\ %s\ [erg\ cm^{-2}\ s^{-1}]$' % xname,
            y = y,
            ylabel = r'$\frac{%s}{%s}$' % (xname, yname),
            kwargs_scatter = dict(marker = 'o', s = 10, edgecolor = 'none', alpha = 0.5, label = ''),
            running_stats = True,
            rs_gaussian_smooth = True,
            rs_gs_fwhm = 0.4,
            kwargs_plot_rs = dict(c = 'k', lw = 2, label = 'Median (run. stats)'),
            rs_errorbar = False,
            suptitle = 'NGals = %d' % len(v['gals']), 
            filename = '%s.png' % k,
            kwargs_legend = dict(fontsize = 8),
        )