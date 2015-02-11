#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import sys
from plot_aux import plotScatterColor, plotScatterColorAxis, plot_text_ax
from califa_scripts import H5SFRData

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
    
try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)
    
H = H5SFRData(h5file)

tSF__T = H.tSF__T
tZ__U = H.tZ__U
xOkMin = H.xOkMin
tauVOkMin = H.tauVOkMin
tauVNebOkMin = H.tauVNebOkMin
tauVNebErrMax = H.tauVNebErrMax
RbinCenter__r = H.RbinCenter__r


# zones
SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
dist_zone__g = H.get_data_h5('dist_zone__g')
tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
tau_V_neb_err__g = H.get_data_h5('tau_V_neb_err__g')
McorSD__g = H.get_data_h5('McorSD__g')
alogZ_mass__Ug = H.get_data_h5('alogZ_mass__Ug')
logZ_neb_S06__g = H.get_data_h5('logZ_neb_S06__g')
tau_V__Tg = H.get_data_h5('tau_V__Tg')

SFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
SFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
tau_V__Trg = H.get_data_h5('tau_V__Trg')
tau_V_neb__rg = H.get_data_h5('tau_V_neb__rg')
alogZ_mass__Urg = H.get_data_h5('alogZ_mass__Urg')
logZ_neb_S06__rg = H.get_data_h5('logZ_neb_S06__rg')
McorSD__rg = H.get_data_h5('McorSD__rg')
morfType_GAL_zones__rg = H.get_data_h5('morfType_GAL_zones__rg')

# galaxy wide quantities replicated by zones
McorSD_GAL_zones__g = H.get_data_h5('McorSD_GAL_zones__g')
morfType_GAL_zones__g = H.get_data_h5('morfType_GAL_zones__g')
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# SKzero = np.log10(1.6e-4)
# SKslope = 1.4
# logSigmaGas = (np.log10(SFRSD_Ha__g * 1e6) - SKzero) / SKslope
# c = np.log10(0.2)
# logDGR = c + np.log10(tau_V_neb__g) - logSigmaGas
# logO_H = logZ_neb_S06__g + np.log10(4.9e-4)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#################################################################################
#################################################################################
#################################################################################

iT = 11 
tSF = tSF__T[iT]
#DGR = 10. ** (-2.35)
DGR = 10. ** (-2.21)
k = 0.2 / DGR

#radii
SigmaGas__rg = k * tau_V__Trg[iT]
#SigmaGas__rg = k * tau_V_neb__rg
r = McorSD__rg / SigmaGas__rg
f_gas__rg = 1. / (1. + r)

#x = np.log10(np.log(1 / f_gas__rg.flatten()))
x = np.ma.log10(f_gas__rg.flatten())
y = alogZ_mass__Urg[-1].flatten() 
#y = logZ_neb_S06__rg.flatten()
z = ((np.ones_like(f_gas__rg).T * RbinCenter__r).T).flatten() #RbinCenter__rg
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
xlabel = r'$\log\ f_{gas}(R)$'
ylabel = r'$\langle \log\ Z_\star \rangle_M(R)$ [$Z_\odot$]'
zlabel = r'R [HLR]'
fname = 'logfgas_alogZmass_radius_%dgals.png' % H.N_gals
ylim = [-1.5, 0.3] 
plotScatterColorAxis(f, xm, ym, zm, xlabel, ylabel, zlabel,
                     None, ylim, contour = False, run_stats = True, OLS = False)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.04))
txt = r'DGR = $10^{%.2f}$  k = %.2f' % (np.log10(DGR), k)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.grid(which = 'major')
f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize=14)
f.savefig(fname)
plt.close(f)

#x = np.log10(np.log(1 / f_gas__rg.flatten()))
x = np.ma.log10(f_gas__rg.flatten())
y = logZ_neb_S06__rg.flatten()
z = ((np.ones_like(f_gas__rg).T * RbinCenter__r).T).flatten() #RbinCenter__rg
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
xlabel = r'$\log\ f_{gas}(R)$'
ylabel = r'$\langle \log\ Z_{neb} \rangle(R)$ [$Z_\odot$]'
zlabel = r'R [HLR]'
fname = 'logfgas_logZneb_radius_%dgals.png' % H.N_gals
ylim = [-0.5, 0.2]
plotScatterColorAxis(f, xm, ym, zm, xlabel, ylabel, zlabel,
                     None, ylim, contour = False, run_stats = True, OLS = False)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
txt = r'DGR = $10^{%.2f}$  k = %.2f' % (np.log10(DGR), k)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.grid(which = 'major')
f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize=14)
f.savefig(fname)
plt.close(f)

x = np.ma.log10(f_gas__rg.flatten())
y = np.ma.log10(McorSD__rg.flatten())
z = ((np.ones_like(f_gas__rg).T * RbinCenter__r).T).flatten() #RbinCenter__rg
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
xlabel = r'$\log\ f_{gas}(R)$'
ylabel = r'$\log\ \langle \mu_\star \rangle(R)$ [$M_\odot \ pc^{-2}$]'
zlabel = r'R [HLR]'
fname = 'logfgas_logMcorSD_radius_%dgals.png' % H.N_gals
ylim = [1, 4.7]
plotScatterColorAxis(f, xm, ym, zm, xlabel, ylabel, zlabel,
                     None, ylim, contour = False, run_stats = True, OLS = False)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
txt = r'DGR = $10^{%.2f}$  k = %.2f' % (np.log10(DGR), k)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.grid(which = 'major')
f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize=14)
f.savefig(fname)
plt.close(f)


x = np.ma.log10(f_gas__rg.flatten())
y = np.ma.log10(SigmaGas__rg.flatten())
z = ((np.ones_like(f_gas__rg).T * RbinCenter__r).T).flatten() #RbinCenter__rg
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
xlabel = r'$\log\ f_{gas}(R)$'
ylabel = r'$\log\ \langle \Sigma_{gas} \rangle(R)$ [$M_\odot \ pc^{-2}$]'
zlabel = r'R [HLR]'
fname = 'logfgas_logSigmaGas_radius_%dgals.png' % H.N_gals
ylim = None
plotScatterColorAxis(f, xm, ym, zm, xlabel, ylabel, zlabel,
                     None, ylim, contour = False, run_stats = True, OLS = False)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
txt = r'DGR = $10^{%.2f}$  k = %.2f' % (np.log10(DGR), k)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.grid(which = 'major')
f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize=14)
f.savefig(fname)
plt.close(f)


x = np.ma.log10(f_gas__rg.flatten())
y = np.ma.log10(McorSD__rg.flatten() / SFRSD__Trg[iT].flatten())
z = ((np.ones_like(f_gas__rg).T * RbinCenter__r).T).flatten() #RbinCenter__rg
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
f = plt.figure()
f.set_size_inches(10, 8)
ax = f.gca()
xlabel = r'$\log\ f_{gas}(R)$'
ylabel = r'$\log\ \langle \mu_\star/\Sigma_{SFR} \rangle$ [yr]'
zlabel = r'R [HLR]'
fname = 'logfgas_McorSDSFRSD_radius_%dgals.png' % H.N_gals
ylim = None
plotScatterColorAxis(f, xm, ym, zm, xlabel, ylabel, zlabel,
                     None, ylim, contour = False, run_stats = True, OLS = False)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
txt = r'DGR = $10^{%.2f}$  k = %.2f' % (np.log10(DGR), k)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.grid(which = 'major')
f.suptitle(r'%d galaxies - tSF:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (H.N_gals, (tSF / 1e6), xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize=14)
f.savefig(fname)
plt.close(f)


#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# x = logZ_neb_S06__g
# y = logDGR
# z = np.log10(McorSD__g) 
# mask = x.mask | y.mask
# xm = x[~mask]
# ym = y[~mask]
# zm = z[~mask]
# xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
# ylabel = r'$\log$ DGR'
# zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
# fname = 'logZneb_logDGR_McorSD.png'
# xlim = None
# ylim = None
# plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname)
#
# x = alogZ_mass__Ug[0]
# x[np.isnan(x)] = np.ma.masked
# y = logDGR
# #z = np.log10(McorSD__g)
# z = dist_zone__g  
# mask = x.mask | y.mask
# xm = x[~mask]
# ym = y[~mask]
# zm = z[~mask]
# xlabel = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (tZ__U[0] / 1e9)
# ylabel = r'$\log$ DGR'
# zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
# fname = 'alogZmass0_logDGR_McorSD.png'
# xlim = None
# ylim = None
# plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname, 
#                  contour = False, run_stats = True, OLS = False)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# x = alogZ_mass__Ug[1]
# y = logDGR
# z = np.log10(McorSD__g)  
# mask = x.mask | y.mask
# xm = x[~mask]
# ym = y[~mask]
# zm = z[~mask]
# xlabel = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$]' % (tZ__U[1] / 1e9)
# ylabel = r'$\log$ DGR'
# zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
# fname = 'alogZmass1_logDGR_McorSD.png'
# xlim = None
# ylim = None
# plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname, 
#                  contour = False, run_stats = True, OLS = False)
# 
# ##
# x = logZ_neb_S06__g
# y = logDGR
# z = np.log10(McorSD_GAL_zones__g) 
# #mask = x.mask | y.mask | (tau_V_neb_err__g > 0.15)
# mask = x.mask | y.mask 
# xm = x[~mask]
# ym = y[~mask]
# zm = z[~mask]
# xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
# ylabel = r'$\log$ DGR'
# zlabel = r'$\log\ M_\star$ [$M_\odot\ pc^{-2}$]'
# fname = 'logZneb_logDGR_McorSDGAL.png'
# xlim = None
# ylim = None
# plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
