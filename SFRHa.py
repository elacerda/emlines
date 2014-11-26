#!/usr/bin/python
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# 
# Calculate the k for SFR(Halpha) = k . L(Halpha)
# 
#     Lacerda@Saco - 9/Jul/2014
#     
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
import numpy as np
from pystarlight.util.base import StarlightBase

#useTrapz = True
useTrapz = False
plot = True
#plot = False
outputImgSuffix = 'pdf'

def add_subplot_axes(ax, rect, axisbg = 'w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x, y, width, height], axisbg = axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2] ** 0.5
    y_labelsize *= rect[3] ** 0.5
    subax.xaxis.set_tick_params(labelsize = x_labelsize)
    subax.yaxis.set_tick_params(labelsize = y_labelsize)

    return subax

if plot:
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MultipleLocator  #, MaxNLocator
    import matplotlib.gridspec as gridspec
    
    mpl.rcParams['font.size']       = 20
    mpl.rcParams['axes.labelsize']  = 20
    mpl.rcParams['axes.titlesize']  = 22
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16 
    mpl.rcParams['font.family']     = 'serif'
    mpl.rcParams['font.serif']      = 'Times New Roman'
    
    plotConf__Z = [
        dict( c = 'b', lw = 0.5),
        dict( c = 'g', lw = 0.5),
        dict( c = 'r', lw = 0.5),
        dict( c = 'y', lw = 0.5),
        dict( c = 'k', lw = 2.),
        dict( c = 'c', lw = 0.5),
    ] 


bases = [ 'Padova1994.chab', 'Padova1994.salp', 'Padova2000.chab', 'Padova2000.salp' ]

baseFile    = '/Users/lacerda/LOCAL/data/Base.bc03.h5'


def create_dx(x):
    dx          = np.zeros_like(x)
    dx[1:]      = (x[1:] - x[:-1])/2.   # dl/2 from right neighbor
    dx[:-1]     += dx[1:]               # dl/2 from left neighbor
    dx[0]       = 2 * dx[0]
    dx[-1]      = 2 * dx[-1]
    return dx


def SFR_parametrize(flux, wl, ages, max_age, useTrapz = False):
    '''
    Find the k parameter in the equation SFR = k [M_sun yr^-1] L(Halpha) [(10^8 L_sun)^-1]
    
    TODO: blablabla
    
    Nh__Zt is obtained for all t in AGES differently from Nh__Z, which consists in the number
    of H-ionizing photons from MAX_AGE till today (t < MAX_AGE).
    ''' 
    from pystarlight.util.constants import L_sun, h, c, yr_sec
    
    cmInAA      = 1e-8          # cm / AA
    lambda_Ha   = 6562.8        # Angstrom
    mask_age    = ages <= max_age
    
    y = flux * wl * cmInAA * L_sun / (h * c)
    
    if useTrapz:
        import scipy.integrate as spi

        qh__Zt = np.trapz(y = y, x = wl, axis = 2) # 1 / Msol
        Nh__Zt = spi.cumtrapz(y = qh__Zt, x = ages, initial = 0, axis = 1) * yr_sec
        Nh__Z = np.trapz(y = qh__Zt[:, mask_age], x = ages[mask_age], axis = 1) * yr_sec
    else:
        qh__Zt = (y * create_dx(l)).sum(axis = 2)
        Nh__Z = (qh__Zt[:, mask_age] * create_dx(ages[mask_age])).sum(axis = 1) * yr_sec
        Nh__Zt = np.cumsum(qh__Zt * create_dx(ages), axis = 1) * yr_sec
         
    k_SFR__Z = 2.226 * lambda_Ha * L_sun * yr_sec / (Nh__Z * h * c) # M_sun / yr
    
    return qh__Zt, Nh__Zt, Nh__Z, k_SFR__Z


if __name__ == '__main__':
    for i, b in enumerate(bases):
        base        = StarlightBase(baseFile, b, hdf5 = True)
        
        max_yr      = base.ageBase[-1]
        max_yr      = 1e7
        mask        = base.l_ssp <= 912         # Angstrom
        mask_age    = base.ageBase <= max_yr
        
        f_ssp   = base.f_ssp[:,:,mask]
        l       = base.l_ssp[mask]
    
        qh__Zt, Nh__Zt, Nh__Z, k_SFR__Z = SFR_parametrize(f_ssp, l, base.ageBase, max_yr, useTrapz)
        
        print b + ':'
        for i, Z in enumerate(base.metBase):
            print '\tZ=%.4f N_H=%e k_SFR=%.2f Msun/yr' % (Z, Nh__Z[i], k_SFR__Z[i])
            
        if plot is False:
            continue
        else:
            f = plt.figure()
            gs = gridspec.GridSpec(2, 2)
            f.set_size_inches(10,10)
            ax1 = plt.subplot(gs[0, 0])
            ax2 = plt.subplot(gs[0, 1])
            ax3 = plt.subplot(gs[1, :])
    
            subpos = [0.68,  0.20, 0.45, 0.35]
            subax = add_subplot_axes(ax2, subpos)
    
            for iZ, Z in enumerate(base.metBase):
                ax1.plot(np.log10(base.ageBase), Nh__Zt[iZ, :] / 1e60,
                         c = plotConf__Z[iZ]['c'], lw = plotConf__Z[iZ]['lw'],  
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                ax2.plot(np.log10(base.ageBase), Nh__Zt[iZ, :] / Nh__Zt[iZ, -1], 
                         c = plotConf__Z[iZ]['c'], lw = plotConf__Z[iZ]['lw'], 
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                subax.plot(np.log10(base.ageBase), Nh__Zt[iZ, :] / Nh__Zt[iZ, -1], 
                           c = plotConf__Z[iZ]['c'], lw = plotConf__Z[iZ]['lw'], 
                           label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                ax3.plot(np.log10(base.ageBase), np.log10(qh__Zt[iZ, :]),
                         c = plotConf__Z[iZ]['c'], lw = plotConf__Z[iZ]['lw'], 
                         label = r'Z $=\ %.2f Z_\odot$' % (Z / base.metBase[4]))
                    
            ax3.legend(loc = 1, fontsize=14, frameon=False)
             
            ax2.axhline(y = 0.95, ls = '--', c = 'k')
            
            ax2.set_ylim([0, 1.1])
            subax.set_xlim([6.4, 7.2])
            subax.set_ylim([0.80, 1.05])
    
            subax.xaxis.set_major_locator(MultipleLocator(0.5))
            subax.xaxis.set_minor_locator(MultipleLocator(0.25))
            subax.yaxis.set_major_locator(MultipleLocator(0.1))
            subax.yaxis.set_minor_locator(MultipleLocator(0.05))
            subax.xaxis.grid(which='minor')
            subax.yaxis.grid(which='minor')
            ax3.xaxis.set_major_locator(MultipleLocator(1))
            ax3.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax3.yaxis.set_major_locator(MultipleLocator(2))
            ax3.yaxis.set_minor_locator(MultipleLocator(1))
            ax3.xaxis.grid(which='minor')
            ax3.yaxis.grid(which='minor')
    
            ax1.set_xlabel(r'$\log\ t\ [yr]$')
            ax1.set_ylabel(r'$\mathcal{N}_H(t)\ [10^{60}\ \gamma\ M_\odot{}^{-1}]$')
            ax2.set_xlabel(r'$\log\ t\ [yr]$')
            ax2.set_ylabel(r'$\mathcal{N}_H(t)/\mathcal{N}_H$')
            ax3.set_xlabel(r'$\log\ t\ [yr]$')
            ax3.set_ylabel(r'$\log\ q_H [s^{-1} M_\odot{}^{-1}]$')
                
            f.tight_layout()
            f.savefig('Nh_logt_metBase_%s.%s' % (b.replace('.', '_'), outputImgSuffix))
            plt.close(f)
