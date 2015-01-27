#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
import numpy as np
import matplotlib as mpl
import sys
from plot_aux import H5SFRData, plotScatterColor

mpl.rcParams['font.size']       = 20
mpl.rcParams['axes.labelsize']  = 20
mpl.rcParams['axes.titlesize']  = 22
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'
    
try:
    h5file = sys.argv[1]
except IndexError:
    print 'usage: %s HDF5FILE' % (sys.argv[0])
    exit(1)
    
H = H5SFRData(h5file)

tSF__T = H.tSF__T[0:20]

# zones
SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
dist_zone__g = H.get_data_h5('dist_zone__g')
tau_V_neb__g = H.get_data_h5('tau_V_neb__g')
tau_V_neb_err__g = H.get_data_h5('tau_V_neb_err__g')
McorSD__g = H.get_data_h5('McorSD__g')
alogZ_mass__g = H.get_data_h5('alogZ_mass__g')
logZ_neb_S06__g = H.get_data_h5('logZ_neb_S06__g')

# galaxy wide quantities replicated by zones
McorSD_GAL_zones__g = H.get_data_h5('McorSD_GAL_zones__g')
morfType_GAL_zones__g = H.get_data_h5('morfType_GAL_zones__g')

SKzero = np.log10(1.6e-4)
SKslope = 1.4
logSigmaGas = (np.log10(SFRSD_Ha__g * 1e6) - SKzero) / SKslope
c = 0 # np.log10(0.2)
logDGR = c + np.log10(tau_V_neb__g) - logSigmaGas
logO_H = logZ_neb_S06__g + np.log10(4.9e-4)

#################################################################################
#################################################################################
#################################################################################
 
x = logZ_neb_S06__g
y = logDGR
z = np.log10(McorSD__g) 
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
ylabel = r'$\log$ DGR'
zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
fname = 'logZneb_logDGR_McorSD.png'
xlim = None
ylim = None
plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname)

##
x = alogZ_mass__g
y = logDGR
z = np.log10(McorSD__g)  
mask = x.mask | y.mask
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
xlabel = r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]'
ylabel = r'$\log$ DGR'
zlabel = r'$\log\ \mu_\star$ [$M_\odot\ pc^{-2}$]'
fname = 'alogZmass_logDGR_McorSD.png'
xlim = None
ylim = None
plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname)

##
x = logZ_neb_S06__g
y = logDGR
z = np.log10(McorSD_GAL_zones__g) 
#mask = x.mask | y.mask | (tau_V_neb_err__g > 0.15)
mask = x.mask | y.mask 
xm = x[~mask]
ym = y[~mask]
zm = z[~mask]
xlabel = r'$\log\ Z_{neb}$ [$Z_\odot$]'
ylabel = r'$\log$ DGR'
zlabel = r'$\log\ M_\star$ [$M_\odot\ pc^{-2}$]'
fname = 'logZneb_logDGR_McorSDGAL.png'
xlim = None
ylim = None
plotScatterColor(xm, ym, zm, xlabel, ylabel, zlabel, xlim, ylim, fname)

