#!/usr/bin/python
from pystarlight.util.base import StarlightBase
from CALIFAUtils.plots import add_subplot_axes
from CALIFAUtils.scripts import linearInterpol, create_dx
from matplotlib.ticker import MultipleLocator
from CALIFAUtils.plots import plot_text_ax
from matplotlib import pyplot as plt
import scipy.integrate as spi
import matplotlib as mpl
import numpy as np

#with_trapz = True
with_trapz = False

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

def ltoeV(l):
    from scipy import constants as sc
    E = (sc.physical_constants['Planck constant in eV s'][0] * sc.physical_constants['speed of light in vacuum'][0]) / l
    out = [ '%.1f' % e for e in E ]
    out[0] = ''
    print l, E, out    
    return out

if __name__ == '__main__':
    bases = [ 'Padova1994.chab', 'Padova1994.salp', 'Padova2000.chab', 'Padova2000.salp' ]
    baseFile = '/Users/lacerda/LOCAL/data/Base.bc03.h5'
    base = StarlightBase(baseFile, bases[-1], hdf5 = True)
    
    l_max = base.l_ssp[-1] #no limit in wavelength
    m_l = base.l_ssp <= l_max
    l = base.l_ssp[m_l]
    #l = base.l_ssp
    
    max_age = 1.43e10 #T_univ at max
    m_t = base.ageBase <= max_age
    t = base.ageBase[m_t]
    i_Z_sun = -2 # base.metBase[-2]
    m__tl = m_t[..., np.newaxis] * m_l
    new_shape = (base.nMet, m_t.sum(), m_l.sum())
    f__Ztl = base.f_ssp[:, m__tl].reshape(new_shape)
    #f__Ztl = base.f_ssp[:, m_t, :][:, :, m_l]
    
    if with_trapz is True:
        q__Zl = np.trapz(y = f__Ztl, x = t, axis = 1)
        r__Ztl = spi.cumtrapz(y = f__Ztl, x = t, axis = 1, initial = 0)
        x_sun__tl = r__Ztl[i_Z_sun, :, :] / q__Zl[i_Z_sun, :]
    else:
        dt = create_dx(t)
        q_sun__l = (f__Ztl[i_Z_sun,:].T * dt).sum(axis = 1)
        r_sun__lt = np.cumsum(f__Ztl[i_Z_sun, :].T * dt, axis = 1)
        x_sun__tl = r_sun__lt.T / q_sun__l
        
    fractions = [ 0.99, 0.97, 0.95 ]
    alpha = [ 0.4, 0.4, 1.0 ]
    c = [ 'g', 'b', 'k' ]
    
    f = plt.figure()
    f.set_size_inches(10,8)
    ax = f.gca()
    subax_1 = add_subplot_axes(ax, [0.55,  0.44, 0.45, 0.30])
    subax_2 = add_subplot_axes(ax, [0.55,  0.1, 0.45, 0.30])
    
    for i in range(len(fractions)):
        age__l = np.empty_like(l)
        
        for j, wl in enumerate(l):
            w = np.where(x_sun__tl[:, j] >= fractions[i])
            age__l[j] = t[w][0]
        
        i_Ha = 122 # base.l_ssp[122] == 915 angs
        t_Ha = linearInterpol(l[121], l[122], age__l[121], age__l[122], 912)
        t_bol = age__l[-1] #base.l_ssp[3025] == 7000 angs
        label_vline = r'$912\ \AA\ $($%.2f$ Myr)' % (t_Ha / 1e6)
        label_hline = r't${}_{bol}\ $($%.2f$ Gyr)' % (t_bol / 1e9)
        ax.plot(l, np.log10(age__l), c = c[i], alpha = alpha[i])
        subax_1.plot(l, np.log10(age__l), c = c[i], alpha = alpha[i])
        subax_2.plot(l, np.log10(age__l), c = c[i], alpha = alpha[i])
        txt = '%.2f' % fractions[i]
        ypos = 1.05
        plot_text_ax(ax, txt, 0.44 + (i * 0.1), ypos, 15, 'bottom', 'right', c[i], alpha = alpha[i])
        #ax.axhline(y = np.log10(t_bol), c = c[i], alpha = 0.4)#, label = label_hline)
        #plot_text_ax(ax, label_vline, 0.99, 0.75, 15, 'bottom', 'right', 'b')
        #plot_text_ax(ax, label_hline, 0.99, 0.70, 15, 'bottom', 'right', 'k')
        print fractions[i], label_vline, label_hline
        del age__l

    ax.set_ylim(6.5, 10.5)
    ax.set_xlim(l[0],7000)
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    #txt = r't${}_{bol}\ =\ %.2f$' % t_bol  
    #plot_text_ax(ax, txt, 0.02, 0.98, 20, 'top', 'left', 'k')
    ax2 = ax.twiny()
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(ltoeV(ax.get_xticks() * 1e-10))

    subax_1.set_xlim(0,1500)
    subax_1.set_ylim(6.5, 10.2)
    subax_1.xaxis.set_major_locator(MultipleLocator(200))
    subax_1.xaxis.set_minor_locator(MultipleLocator(40))
    subax_1.yaxis.set_major_locator(MultipleLocator(1))
    subax_1.yaxis.set_minor_locator(MultipleLocator(0.2))

    subax_2.set_xlim(2000,3000)
    subax_2.set_ylim(8.0, 10.2)
    subax_2.xaxis.set_major_locator(MultipleLocator(200))
    subax_2.xaxis.set_minor_locator(MultipleLocator(40))
    subax_2.yaxis.set_major_locator(MultipleLocator(.5))
    subax_2.yaxis.set_minor_locator(MultipleLocator(0.1))
    
    xlim = ax.get_xlim()
    xlim_size = xlim[1] - xlim[0]
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # ypos = 1.06
    # fs = 10
    # Lya_pos = 1215. / xlim_size
    # plot_text_ax(ax, r' Ly$\alpha$', Lya_pos, ypos, fs, 'top', 'left', 'k')
    # Ha_pos = 912. / xlim_size
    # plot_text_ax(ax, r' H${}^{+}$', Ha_pos, ypos, fs, 'top', 'left', 'k')
    # He_pos = 504 / xlim_size
    # plot_text_ax(ax, r' He${}^{+}$', He_pos, ypos, fs, 'top', 'left', 'k')
    # Hee_pos = 227.9 / xlim_size
    # plot_text_ax(ax, r' He${}^{2+}$', Hee_pos, ypos, fs, 'top', 'left', 'k')
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    for tmp_ax in [ ax, subax_1, subax_2 ]:
        #Lyman limit
        tmp_ax.axvline(x = 912, ls = '--', c = 'k')
        #lyman alpha
        tmp_ax.axvline(x = 1215, ls = '--', c = 'k')
        #Helium 1e
        tmp_ax.axvline(x = 504, ls = '--', c = 'k')
        #Helium 2e
        tmp_ax.axvline(x = 227.9, ls = '--', c = 'k')
        tmp_ax.xaxis.grid(which='major')
        tmp_ax.yaxis.grid(which='major')
    
    ax.set_xlabel(r'wavelength [$\AA$]')
    ax.set_ylabel(r'$\log$ t [yr]')
    ax.legend(fontsize = 15)
    f.tight_layout()
    f.savefig('age_wl.png')