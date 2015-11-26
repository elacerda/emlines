#!/usr/bin/python
#
# Lacerda@Saco - 23/Jun/2014
#
import sys
import time
import numpy as np
import argparse as ap
import CALIFAUtils as C
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.plots import DrawHLRCircle
from matplotlib.pyplot import MultipleLocator
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.scripts import calc_running_stats
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import DrawHLRCircleInSDSSImage
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : None,
        'califaID' : 'K0073',
        'itSF' : 11,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--califaID', '-g',
                        metavar = 'FILE',
                        type = str,
                        default = default['califaID'])
    parser.add_argument('--itSF', '-T',
                        help = 'SF age index.',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = "FILENAME",
                        type = str,
                        default = default['output'])

    return parser.parse_args()


def newpage(NRows, NCols, K, galaxyImgFile = None):
    from os.path import isfile
    f, axArr = plt.subplots(NRows, NCols)
    f.set_size_inches((NCols * 3, NRows * 2.5))
    #turn off all the axes        
    for ax in f.axes:
        ax.set_axis_off()
    #fist plot - Galaxy Image
    if galaxyImgFile is not None and isfile(galaxyImgFile):
        ax = axON(axArr, 0, 0)
        galimg = plt.imread(galaxyImgFile)[::-1, :, :]
        plt.setp(ax.get_xticklabels(), visible = False)
        plt.setp(ax.get_yticklabels(), visible = False)
        ax.imshow(galimg, origin = 'lower')
        DrawHLRCircleInSDSSImage(ax, K.HLR_pix, K.pa, K.ba)
    return f, axArr

def axON(ax_array, row, col):
    ax = ax_array[row, col]
    ax.set_axis_on()
    return ax

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    Zsun = 0.019

    t_init_gal = time.clock()
    
    args = parser_args()
    debug = args.debug
    iT = args.itSF

    H = C.H5SFRData(args.hdf5)
    paths = C.CALIFAPaths()

    C.debug_var(debug, args = args)

    if (len(np.where(H.califaIDs == args.califaID)[0]) == 0):
        exit('<<< plot: %s: no data in HDF5 file.' % args.califaID)
        
    K = C.read_one_cube(args.califaID, EL = True, GP = True)
    
    age = H.tSF__T[iT]
    ageMyr = age / 1e6
    
    # Setup elliptical-rings geometry
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)
    
    tipos, tipo, tipo_m, tipo_p = C.get_morfologia(args.califaID)
    
    galaxyImgFile = paths.get_image_file(args.califaID)    
    
    if args.output is None:
        output_filename = '%s_dossier_v2.pdf' % args.califaID
    else:
        output_filename = args.output
    
    lines__Lyx = np.ma.empty((K.EL.N_line, K.N_y, K.N_x))
    lines__Lr = np.ma.empty((K.EL.N_line, H.NRbins))
    log_lines__Lr = np.ma.empty((K.EL.N_line, H.NRbins))
    lines_err__Lr = np.ma.empty((K.EL.N_line, H.NRbins))

    # This file will produce a big pdf with a lot of plots.
    with PdfPages(output_filename) as pdf:
        NRows = 4
        NCols = 4
        suptitle = r'%s | (b/a) pycasso: %.2f | (b/a) masterlist: %s | morphology type: %s' % (args.califaID, K.ba, np.float(K.masterListData['ba']), tipo)

        ##########################
        ######### PAGE 1 #########
        ##########################
        f, axArr = newpage(NRows, NCols, K, galaxyImgFile)
        f.suptitle(suptitle)
        row, col = 1, 0
        for line in ['4861', '5007', '6563', '6583']:
        #for i_l, line in enumerate(K.EL.lines):
            i_l = K.EL.lines.index(line) 
            if col == NCols:
                row += 1
                col = 0
            ax = axON(axArr, row, col)
            ax.set_title(r'$\log$ %s' % line)
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
            lines__Lyx[i_l] = K.zoneToYX(K.EL.flux[i_l], extensive = True, surface_density = False)
            #aux = lines__Lyx[i_l].compressed().ravel()
            #m_aux = aux > 0
            #prc = np.percentile(aux[m_aux], [2, 98])
            #print prc
            kw_imshow = dict(
                interpolation = 'nearest',
                origin = 'lower',
                aspect = 'auto',
                #vmin = prc[0],
                #vmax = prc[1],
                #vmin = -17.5,
                #vmax = -14.5,
                cmap = 'Blues',
            )
            im = ax.imshow(np.ma.log10(lines__Lyx[i_l]), **kw_imshow)
            DrawHLRCircle(ax, K, color = 'k')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size = '5%', pad = 0)
            cb = plt.colorbar(im, cax = cax)
            plt.setp(cb.ax.get_yticklabels(), fontsize = 6)
            col += 1
            ax = axON(axArr, row, col)
            #ax.set_title(r'$\log$ %s' % line)
            
            log_lines__Lr[i_l] = K.radialProfile(np.ma.log10(lines__Lyx[i_l]), H.Rbin__r, rad_scale = K.HLR_pix)
            lines__Lr[i_l] = K.radialProfile(lines__Lyx[i_l], H.Rbin__r, rad_scale = K.HLR_pix)
            lines_err__Lr[i_l] = K.radialProfile(lines__Lyx[i_l], H.Rbin__r, rad_scale = K.HLR_pix, mode = 'std')
            #ax.errorbar(H.RbinCenter__r, np.log10(lines__Lr[i_l]), np.log10(lines_err__Lr[i_l]), ls = '--', marker='o', label = r'$\log\ \langle F_{l} \rangle$')
            #ax.errorbar(H.RbinCenter__r, log_lines__r, log_lines_err__r, ls = '-', marker='o', label = r'$\langle \log\ F_{l}\rangle$')
            ax.scatter(K.zoneDistance_HLR, np.ma.log10(K.EL.flux[i_l]), marker = 'o', c = '0.5', s = 2, edgecolor = 'none', alpha = 0.4, label = '')
            ax.plot(H.RbinCenter__r, log_lines__Lr[i_l], ls = '-', marker = '.', label = r'$\langle \log\ F_{l}\rangle$')
            ax.plot(H.RbinCenter__r, np.ma.log10(lines__Lr[i_l]), ls = '--', marker = '.', label = r'$\log\ \langle F_{l} \rangle$')
            ax.set_title(r'mean flux [erg cm${}^{-2}$ s${}^{-1}$]', fontsize = 7)
            #ax.yaxis.tick_right()
            #ax.yaxis.set_label_position('right')
            #ax.yaxis.labelpad = 0.02
            ax.legend(loc = 'upper right', fontsize = 6, frameon = False)
            ax.set_xlim(0,2)
            ax.set_ylim(-17.5,-14.5)
            plt.setp(ax.get_xticklabels(), fontsize = 6)
            plt.setp(ax.get_yticklabels(), fontsize = 6)
            ax.grid()
            col += 1
                             #d_x   d_y    #s_x   #s_y
        #new_ax = f.add_axes([0.125, 0.055, 0.775, 0.03])
        #cb = plt.colorbar(im, cax = new_ax, orientation = 'horizontal')
        # here ends the first page
        #f.tight_layout()
        
        f.subplots_adjust(hspace = 0.2, wspace = 0.3)
        pdf.savefig(f)
        plt.close(f)
        ##########################
        ######### PAGE 2 #########
        ##########################
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # f, axArr = newpage(NRows, NCols, K, galaxyImgFile)
        # #tau_V
        # row, col = 0, 1
        # f.suptitle('LINES RADIAL PROFILES [HLR] | %s' % suptitle)
        # for i_l, line in enumerate(K.EL.lines):
        #     if col == NCols:
        #         row += 1
        #         col = 0
        #     ax = axON(axArr, row, col)
        #     ax.set_title(r'$\log$ %s' % line)
        #     
        #     log_lines__Lr[i_l] = K.radialProfile(np.ma.log10(lines__Lyx[i_l]), H.Rbin__r, rad_scale = K.HLR_pix)
        #     lines__Lr[i_l] = K.radialProfile(lines__Lyx[i_l], H.Rbin__r, rad_scale = K.HLR_pix)
        #     lines_err__Lr[i_l] = K.radialProfile(lines__Lyx[i_l], H.Rbin__r, rad_scale = K.HLR_pix, mode = 'std')
        #     #ax.errorbar(H.RbinCenter__r, np.log10(lines__Lr[i_l]), np.log10(lines_err__Lr[i_l]), ls = '--', marker='o', label = r'$\log\ \langle F_{l} \rangle$')
        #     #ax.errorbar(H.RbinCenter__r, log_lines__r, log_lines_err__r, ls = '-', marker='o', label = r'$\langle \log\ F_{l}\rangle$')
        #     #ax.scatter(K.zoneDistance_HLR, np.ma.log10(K.EL.flux[i_l]), marker = 'o', c = '0.9', s = 1, edgecolor = 'none', alpha = 0.1, label = '')
        #     ax.plot(H.RbinCenter__r, log_lines__Lr[i_l], ls = '-', marker = '.', label = r'$\langle \log\ F_{l}\rangle$')
        #     ax.plot(H.RbinCenter__r, np.ma.log10(lines__Lr[i_l]), ls = '--', marker = '.', label = r'$\log\ \langle F_{l} \rangle$')
        #     ax.set_ylabel(r'mean flux [erg cm${}^{-2}$ s${}^{-1}$]', fontsize = 7)
        #     ax.yaxis.labelpad = 0.02
        #     ax.legend(loc = 'upper right', fontsize = 6, frameon = False)
        #     #ax.set_ylim(-24,-19)
        #     plt.setp(ax.get_xticklabels(), fontsize = 6)
        #     plt.setp(ax.get_yticklabels(), fontsize = 6)
        #     col += 1
        # pdf.savefig(f)
        # plt.close(f)
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
                

sys.exit('bai')



    
# global
xOkMin = H.xOkMin
tauVOkMin = H.tauVOkMin
tauVNebOkMin = H.tauVNebOkMin
tauVNebErrMax = H.tauVNebErrMax
RbinCenter__r = H.RbinCenter__r
# ALL gal
##stellar
x_Y__g = H.x_Y__Tg[iT]
aSFRSD__rg = H.aSFRSD__Trg[iT]
tau_V__rg = H.tau_V__Trg[iT]
##nebular
aSFRSD_Ha__rg = H.aSFRSD_Ha__rg
tau_V_neb__rg = H.tau_V_neb__rg
# one gal
##stellar
tau_V__z = getattr(H, '%s_tau_V__Tg' % args.califaID)[iT]
atau_V__r = getattr(H, '%s_tau_V__Trg' % args.califaID)[iT]
SFRSD__z = getattr(H, '%s_SFRSD__Tg' % args.califaID)[iT]
aSFRSD__r = getattr(H, '%s_aSFRSD__Trg' % args.califaID)[iT]
x_Y__z = getattr(H, '%s_x_Y__Tg' % args.califaID)[iT]
##nebular
EW_Ha__z = getattr(H, '%s_EW_Ha__g' % args.califaID)
EW_Hb__z = getattr(H, '%s_EW_Hb__g' % args.califaID)
tau_V_neb__z = getattr(H, '%s_tau_V_neb__g' % args.califaID)
atau_V_neb__r = getattr(H, '%s_tau_V_neb__rg' % args.califaID)
tau_V_neb_err__z = getattr(H, '%s_tau_V_neb_err__g' % args.califaID)
SFRSD_Ha__z = getattr(H, '%s_SFRSD_Ha__g' % args.califaID)
aSFRSD_Ha__r = getattr(H, '%s_aSFRSD_Ha__rg' % args.califaID)

##DGR
SKzero = 1.6e-4
SKslope = 1.4
dustdim = 0.2
#DGR__g = dustdim * H.tau_V__Tg[iT] / (H.SFRSD__Tg[iT]/SKzero) ** (1./SKslope)
#DGR_Ha__g = dustdim * H.tau_V_neb__g / (H.SFRSD_Ha__g/SKzero) ** (1./SKslope)
#DGR__rg = dustdim * H.tau_V_neb__rg / (H.aSFRSD_Ha__rg/1.6e-4) ** (1./SKslope)
DGR = 10. ** (-2.21)
#DGR__rg = 10. ** (-2.21 + 1./(H.O_O3N2_M13__rg ** 2.))
#DGR__rg = 10. ** (-0.21)
k = dustdim / DGR
#k__g = dustdim / DGR__g
#k__rg = dustdim / DGR__rg
#SigmaGas__g = k * H.tau_V_neb__g
SigmaGas__g = k * H.tau_V__Tg[iT]
#SigmaGas__g = k__g * H.tau_V_neb__g
#SigmaGas__g = k__g * H.tau_V__Tg[iT]
#SigmaGas__rg = k * H.tau_V_neb__rg
SigmaGas__rg = k * H.tau_V__Trg[iT]
#SigmaGas__rg = k * H.tau_V_neb__rg
#SigmaGas__rg = k__rg * H.tau_V_neb__rg
#SigmaGas__rg = k__rg * H.tau_V__Trg[iT]
f_gas__g = 1. / (1. + (H.McorSD__Tg[iT] / SigmaGas__g))
f_gas__rg = 1. / (1. + (H.McorSD__Trg[iT] / SigmaGas__rg))
RbinCenter__rg = ((np.ones_like(f_gas__rg).T * H.RbinCenter__r).T).flatten()

#stellar
tau_V__yx = K.zoneToYX(tau_V__z, extensive = False)
SFRSD__yx = K.zoneToYX(SFRSD__z, extensive = False)
x_Y__yx = K.zoneToYX(x_Y__z, extensive = False)
#nebular
tau_V_neb__yx = K.zoneToYX(tau_V_neb__z, extensive = False)
SFRSD_Ha__yx = K.zoneToYX(SFRSD_Ha__z, extensive = False)
EW_Ha__yx = K.zoneToYX(EW_Ha__z, extensive = False)
EW_Hb__yx = K.zoneToYX(EW_Hb__z, extensive = False)
tau_V_neb_err__yx = K.zoneToYX(tau_V_neb_err__z.data, extensive = False)
F_obs_Ha__z = K.EL.flux[K.EL.lines.index('6563'), :]
F_obs_Ha__yx = K.zoneToYX(F_obs_Ha__z, extensive = True)
#mixed
deltaTau__z = tau_V_neb__z - tau_V__z 
deltaTau__yx = K.zoneToYX(deltaTau__z, extensive = False)

#DGR    
SigmaGas__r = k * atau_V__r
f_gas = 1. / (1. + (K.Mcor__z / (K.zoneArea_pc2 * SigmaGas__g)))

t_calc = time.clock()
print 'calc: elapsed time: %.2f' % (t_calc - t_init_gal)












ax = axArr[0, 1]
ax.set_axis_on()
xlabel = r'EW(H$\alpha$) [$\AA$]'
im = ax.imshow(EW_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 20, vmin = 3)
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

sigma_dev = 1
ax = axArr[0, 2]
ax.set_axis_on()
xlabel = r'$F_{obs}(H\alpha)$ [erg cm${}^{-2}$ s${}^{-1}$]'
vmax = F_obs_Ha__yx.mean() + sigma_dev * F_obs_Ha__yx.std()
vmin = 0
im = ax.imshow(F_obs_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[0, 3]
ax.set_axis_on()
xlabel = r'$\epsilon\tau_V^{neb}$'
im = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = 1)
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1, 0]
ax.set_axis_on()
x = np.ma.log10(tau_V__rg.flatten())
y = np.ma.log10(aSFRSD__rg.flatten() * 1e6)
mask = x.mask | y.mask  
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
xlabel = r'$\log\ \tau_V^{\star}(R)$'
ylabel = r'$\log\ \langle \Sigma_{SFR}^\star(t_\star, R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$'
xlim = [-1.5, xm.max()]
ylim = [-3.5, 1]
sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
nBox = 20
dxBox = (xm.max() - xm.min()) / (nBox - 1.)
X = x[~mask]
Y = y[~mask]
aux = calc_running_stats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
xbinCenter = aux[0]
xMedian = aux[1]
xMean = aux[2]
xStd = aux[3]
yMedian = aux[4]
yMean = aux[5]
yStd = aux[6]
nInBin = aux[7]
xPrc = aux[8]
yPrc = aux[9]
ax.plot(xMedian, yMedian, 'k', lw = 2)
ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
txt = '%.2f Myr' % (age / 1e6)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.98, pos_y = 0.02, fs = 14)
##########################
x = np.ma.log10(atau_V__r)
y = np.ma.log10(aSFRSD__r * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
# vmax = xr.mean() + 2. * xr.std()
# vmin = xr.mean() - 2. * xr.std()
# ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'jet_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
cb.set_label(r'R [HLR]')
##########################
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.125))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.125))
ax.grid(which = 'major')

sigma_dev = 3.

ax = axArr[1, 1]
ax.set_axis_on()
x = np.ma.log10(tau_V__yx)
y = np.ma.log10(SFRSD__yx * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
ylim = [-3.5, 1]
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# vmax = y.mean() + sigma_dev * y.std()
# vmin = y.mean() - sigma_dev * y.std()
# im = ax.imshow(ym, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin, vmax = vmax)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
txt = 'not masked: %d' % (~mask).sum()
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1, 2]
ax.set_axis_on()
x = np.ma.log10(tau_V__yx)
y = np.ma.log10(SFRSD__yx * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
xlabel = r'$\log\ \tau_V^\star$' 
ylim = [-3.5, 1]
vmin = xm.min()
im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
#im = ax.imshow(x, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1, 3]
ax.set_axis_on()
xlabel = r'$x_Y [\%]$'
vmin = xOkMin * 100.
im = ax.imshow(100. * x_Y__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(xlabel, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 0]
ax.set_axis_on()
x = np.ma.log10(tau_V_neb__rg.flatten())
y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
xlabel = r'$\log\ \tau_V^{neb}(R)$'
ylabel = r'$\log\ \langle \Sigma_{SFR}^{neb}(R)\rangle\ [M_\odot yr^{-1} kpc^{-2}]$' 
xlim = [-1.5, xm.max()]
ylim = [-3.5, 1]
sc = ax.scatter(x, y, c = 'grey', marker = 'o', s = 10., edgecolor = 'none', alpha = 0.4)
nBox = 20
dxBox = (xm.max() - xm.min()) / (nBox - 1.)
X = x[~mask]
Y = y[~mask]
aux = calc_running_stats(X, Y, dxBox = dxBox, xbinIni = X.min(), xbinFin = X.max(), xbinStep = dxBox)
xbinCenter = aux[0]
xMedian = aux[1]
xMean = aux[2]
xStd = aux[3]
yMedian = aux[4]
yMean = aux[5]
yStd = aux[6]
nInBin = aux[7]
xPrc = aux[8]
yPrc = aux[9]
ax.plot(xMedian, yMedian, 'k', lw = 2)
ax.plot(xPrc[0], yPrc[0], 'k--', lw = 2)
ax.plot(xPrc[1], yPrc[1], 'k--', lw = 2)
txt = '%.2f Myr' % (age / 1e6)
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
a_ols, b_ols, sigma_a_ols, sigma_b_ols = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.98, pos_y = 0.02, fs = 14)
##########################
x = np.ma.log10(atau_V_neb__r)
y = np.ma.log10(aSFRSD_Ha__r * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
# vmax = xr.mean() + 2. * xr.std()
# vmin = xr.mean() - 2. * xr.std()
# ax.scatter(xm, ym, c = xr, cmap = 'winter_r', marker = 'o', s = 30, edgecolor = 'black', vmax = vmax, vmin = vmin)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
sc = ax.scatter(xm, ym, c = H.RbinCenter__r, cmap = 'jet_r', marker = 'o', s = 30, edgecolor = 'black', vmax = H.RbinFin, vmin = H.RbinIni)
cb = f.colorbar(ax = ax, mappable = sc, use_gridspec = True)
cb.set_label(r'R [HLR]')
##########################
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.125))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.125))
ax.grid(which = 'major')

ax = axArr[2, 1]
ax.set_axis_on()
x = np.ma.log10(tau_V_neb__yx)
y = np.ma.log10(SFRSD_Ha__yx * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
label = r'$\log\ \Sigma_{SFR}^{neb} [M_\odot yr^{-1} kpc^{-2}]$' 
ylim = [-3.5, 1]
im = ax.imshow(y, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
txt = 'not masked: %d' % (~mask).sum()
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.set_title(label, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 2]
ax.set_axis_on()
x = np.ma.log10(tau_V_neb__yx)
y = np.ma.log10(SFRSD_Ha__yx * 1e6)
mask = x.mask | y.mask
xm = np.ma.masked_array(x, mask = mask)
ym = np.ma.masked_array(y, mask = mask)
label = r'$\log\ \tau_V^{neb}$' 
vmin = np.log10(tauVNebOkMin)
im = ax.imshow(xm, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmin = vmin)
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
ax.set_title(label, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 3]
ax.set_axis_on()
label = r'$\delta\ \tau_V$'
im = ax.imshow(deltaTau__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
mask = deltaTau__yx.mask
txt = 'not masked: %d' % (~mask).sum()
plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
ax.set_title(label, y = -0.15)
ax.grid()
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# x = np.ma.log10(tau_V__Tz[iT])
# y = np.ma.log10(SFRSD__Tz[iT] * 1e6)
# mask = x.mask | y.mask
# xm = np.ma.masked_array(x, mask = mask)
# ym = np.ma.masked_array(y, mask = mask)
# xr = xm + np.cos(a_ols) + ym * np.sin(a_ols)
# xr__yx = K.zoneToYX(xr, extensive = False)
# xlabel = r'xr'
# #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# # vmax = xr__yx.mean() + 2. * xr__yx.std()
# # vmin = xr__yx.mean() - 2. * xr__yx.std()
# # im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r', vmax = vmax, vmin = vmin)
# #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# im = ax.imshow(xr__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'winter_r')
# DrawHLRCircle(ax, K, color = 'black', lw = 2.5)
# txt = 'Nz: %d' % K.N_zone 
# plot_text_ax(ax, txt, 0.02, 0.98, 14, 'top', 'left')
# txt = 'masked: %d' % mask.sum() 
# plot_text_ax(ax, txt, 0.02, 0.92, 14, 'top', 'left')
# ax.set_title(xlabel, y=-0.15)
# ax.grid()
# f.colorbar(ax = ax, mappable = im, use_gridspec = True)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

f.suptitle(r'%s - morph:%s  b/a:%.2f  age:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (args.califaID, tipos, ba, ageMyr, xOkMin * 100., tauVOkMin, tauVNebOkMin, tauVNebErrMax), fontsize = 24)
f.subplots_adjust(left = 0.07, bottom = 0.1, right = 0.99, wspace = 0.1, top = 0.9)
f.savefig('%s_mosaic.png' % args.califaID)
plt.close(f)
t_plot = time.clock()
print 'plot: elapsed time: %.2f' % (t_plot - t_calc)
print 'total: elapsed time: %.2f' % (t_plot - t_init_gal)
