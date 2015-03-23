import sys
import numpy as np
import argparse as ap
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.plots import plot_text_ax
from CALIFAUtils.plots import calcRunningStats
from CALIFAUtils.globals import califa_work_dir
from CALIFAUtils.scripts import H5SFRData
from CALIFAUtils.scripts import read_one_cube
from CALIFAUtils.scripts import DrawHLRCircle
from CALIFAUtils.scripts import DrawHLRCircleInSDSSImage
from CALIFAUtils.scripts import get_morfologia


mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 18
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def parser_args():
    default = {
        'debug' : False,
        'hdf5' : None,
        'califaID' : 'K0073',
        'itSF' : 11,
    }

    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
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
                        help = 'age index',
                        metavar = '',
                        type = int,
                        default = default['itSF'])

    return parser.parse_args()

def print_args(args):
    for k, v in args.__dict__.iteritems():
        print k, v 

if __name__ == '__main__':
    args = parser_args()
    print_args()

    imgDir = califa_work_dir + 'images/'
    Zsun = 0.019
    
    H = H5SFRData(args.hdf5)
    tSF__T = H.tSF__T
    ageMyr = tSF__T[args.itSF] / 1e6
    
    if (len(np.where(H.califaIDs == args.califaID)[0]) == 0):
        exit('<<< plot: %s: no data.' % args.califaID)
    
    K = read_one_cube(args.califaID, EL = True)
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)

    f, axArr = plt.subplots(4, 4)
    f.set_size_inches(24, 20)
    
    for ax in f.axes:
        ax.set_axis_off()
    
    ax = axArr[0, 0]
    ax.set_axis_on()
    galaxyImgFile = imgDir + args.califaID + '.jpg'
    galimg = plt.imread(galaxyImgFile)
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    ax.imshow(galimg)
    
    ax = axArr[0, 1]
    ax.set_axis_on()
    ax.plot(H.RbinCenter__r, H.get_prop_gal(H.tau_V_neb__g, args.califaID), 'o-k')
    #ax.tick_params(axis='x', pad=30)
    ax.set_title(r'$\tau_V^{neb}(R)$')
    
    ax = axArr[0, 2]
    ax.set_axis_on()
    L_int_Ha__yx = K.zoneToYX(L_int_Ha__z, extensive = True)
    L_int_Ha__r = K.radialProfile(L_int_Ha__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, L_int_Ha__r, 'o-k')
    ax.set_title(r'$L_{H\alpha}^{int}(R)$')

    ax = axArr[0, 3]
    ax.set_axis_on()
    #Lobn__yx            = K.zoneToYX(K.Lobn__z, extensive = True)
    #logZNeb__r          = radialProfileWeighted(logZ_neb__yx, Lobn__yx, Rbin__r, K.HLR_pix, K.radialProfile)
    logZ_neb__yx = K.zoneToYX(K.EL.logZ_neb_S06__z, extensive = False)
    logZNeb__r = K.radialProfile(logZ_neb__yx, Rbin__r, rad_scale = K.HLR_pix)
    ax.plot(RbinCenter__r, logZNeb__r, 'o-k')
    #ax.set_title(r'$\langle \log\ Z_{neb}\rangle_L (HLR)$')
    ax.set_title(r'$\log\ Z_{neb}(R)$')

    ax = axArr[1, 1]
    ax.set_axis_on()
    ax.set_title(r'$\tau_V^{neb}$')
    sigma = tau_V_neb__yx.std()
    mean = tau_V_neb__yx.mean()
    vmax = mean + 2. * sigma
    vmin = mean - 2. * sigma
    im = ax.imshow(tau_V_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[2, 1]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (\tau_V^{neb})$')
    sigma = tau_V_neb_err__yx.std()
    mean = tau_V_neb_err__yx.mean()
    vmax = mean + 2. * sigma
    im = ax.imshow(tau_V_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[3, 1]
    ax.set_axis_on()
    ax.set_title(r'$F^{H\alpha}_{H\beta} / \epsilon(F^{H\alpha}_{H\beta})$')
    HaHb__yx = K.zoneToYX(K.EL.HaHb__z, extensive = True)
    err_HaHb__yx = K.zoneToYX(K.EL.HaHb_err__z, extensive = True)
    signalToNoise = np.abs(HaHb__yx) / np.abs(err_HaHb__yx) 
    sigma = signalToNoise.std()
    mean = signalToNoise.mean()
    vmax = mean + 2. * sigma
    im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)

    ax = axArr[1, 2]
    ax.set_axis_on()
    ax.set_title(r'$L_{H\alpha}^{int}$')
    sigma = L_int_Ha__yx.std()
    mean = L_int_Ha__yx.mean()
    vmax = mean + sigma
    vmin = mean - sigma
    im = ax.imshow(L_int_Ha__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[2, 2]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (L_{H\alpha}^{int})$')
    L_int_Ha_err__yx = K.zoneToYX(L_int_Ha_err__z, extensive = True)
    sigma = L_int_Ha_err__yx.std()
    mean = L_int_Ha_err__yx.mean()
    vmax = mean + sigma
    im = ax.imshow(L_int_Ha_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[3, 2]
    ax.set_axis_on()
    ax.set_title(r'$L_{H\alpha}^{int} / \epsilon(L_{H\alpha}^{int})$')
    signalToNoise = np.abs(L_int_Ha__yx) / np.abs(L_int_Ha_err__yx) 
    sigma = signalToNoise.std()
    mean = signalToNoise.mean()
    vmax = mean + sigma
    im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)

    ax = axArr[1, 3]
    ax.set_axis_on()
    ax.set_title(r'$\log\ Z_{neb}$')
    sigma = logZ_neb__yx.std()
    mean = logZ_neb__yx.mean()
    vmax = mean + sigma
    vmin = mean - sigma
    im = ax.imshow(logZ_neb__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax, vmin = vmin)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[2, 3]
    ax.set_axis_on()
    ax.set_title(r'$\epsilon (log\ Z_{neb})$')
    logZ_neb_err__yx = K.zoneToYX(K.EL.logZ_neb_S06_err__z, extensive = False)
    sigma = logZ_neb_err__yx.std()
    mean = logZ_neb_err__yx.mean()
    vmax = mean + sigma
    im = ax.imshow(logZ_neb_err__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
    ax = axArr[3, 3]
    ax.set_axis_on()
    ax.set_title(r'$\log\ Z_{neb} / \epsilon(log\ Z_{neb})$')
    signalToNoise = np.abs(logZ_neb__yx) / np.abs(logZ_neb_err__yx) 
    sigma = signalToNoise.std()
    mean = signalToNoise.mean()
    vmax = mean + sigma
    im = ax.imshow(signalToNoise, origin = 'lower', interpolation = 'nearest', aspect = 'auto', vmax = vmax)
    f.colorbar(ax = ax, mappable = im, use_gridspec = False)
    
#    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
    f.savefig(K.califaID + '_' + versionSuffix + '_nebular.png')
    plt.close(f)