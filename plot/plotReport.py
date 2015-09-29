#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014 
#                 - 10/Sept/2015
#
from matplotlib.ticker import MultipleLocator
import sys
import CALIFAUtils as C
from scipy import stats as st
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
#from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.plots import plot_linreg_params
from CALIFAUtils.plots import plot_text_ax
#from CALIFAUtils.plots import plot_zbins
#from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'scatter' : False,
        'hdf5' : None,
        'output' : '',
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--scatter', '-s',
                        action = 'store_true',
                        default = default['scatter'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
    parser.add_argument('--maskradius', '-R',
                        help = 'initial RDisc value in HLR',
                        metavar = 'NUM',
                        type = float,
                        default = default['maskradius'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default['output'])

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    H = C.H5SFRData(args.hdf5)
    
    minR = 0
    fnamesuffix = '.pdf'
    
    if args.maskradius is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
        maxR = H.Rbin__r[-1]
        maxR = 3
        maskRadiusOk__g = (H.zone_dist_HLR__g >= minR) & (H.zone_dist_HLR__g <= maxR) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * ((H.RbinCenter__r >= minR) & (H.RbinCenter__r <= maxR))).T
        fnamesuffix = '_maskradius.pdf'
        
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1]
        fnamesuffix = '_%s%s' % (gals_txt, fnamesuffix)
        #integrated
        l_gals, _ = C.sort_gals(args.slice_gals)
        l_gals = sorted(l_gals)
        gals_slice__integr = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
        for g in l_gals:
            i = H.califaIDs_all.tolist().index(g)
            gals_slice__integr[i] = True
            
    txt_suptitle = r'$\Longrightarrow$ %s  NGals:%d  $x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f ' % (gals_txt, N_gals, H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     txt_suptitle = r'$\Longrightarrow$  NGals:%d' % (N_gals)
#     
#     maskmtype_ok = {
#         'list_E' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_S0' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_S0a' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sa' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sab' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sb' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sbc' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sc' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Scd' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sd' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sdm' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Sm' : dict(gals = [], zones = None, radius = None, integrated = None),
#         'list_Irr' : dict(gals = [], zones = None, radius = None, integrated = None),
#     }
#     
#     path = '/Users/lacerda/CALIFA/list_v20_q050.d15a/'
#     for k, v in maskmtype_ok.iteritems():
#         filegals = path + k + '.txt'
#         v['gals'], _ = C.sort_gals(filegals, order = 1)
#         v['zones'] = H.get_mask_zones_list(v['gals'])
#         v['radius'] = H.get_mask_radius_list(v['gals'])
#         v['integrated'] = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
#         for g in v['gals']:
#             i = H.califaIDs_all.tolist().index(g)
#             v['integrated'][i] = True
#     
#     # new lists using the old lists
#     for nl in [ 'S0_S0a', 'Sa_Sab', 'Sc_Scd_Sd' ]:
#         mtypes = nl.split('_')
#         k1 = 'list_' + mtypes[0]
#         k2 = 'list_' + mtypes[1]
#         nk = 'list_' + nl
#         maskmtype_ok[nk] = {}
#         for k in maskmtype_ok[k1].iterkeys():
#             if k != 'gals': 
#                 m1 = maskmtype_ok[k1][k]
#                 m2 = maskmtype_ok[k2][k]
#                 nv = m1 | m2 
#             else:
#                 l1 = maskmtype_ok[k1][k].tolist()
#                 l2 = maskmtype_ok[k2][k].tolist()
#                 nv = sorted(l1 + l2)  
#             maskmtype_ok[nk][k] = nv
#             
#     mtype_bins = [ 'E', 'S0_S0a', 'Sa_Sab', 'Sb', 'Sbc', 'Sc_Scd_Sd' ]
#     mtype_colors = {
#         'E' : 'brown',
#         'S0_S0a' : 'red',
#         'Sa_Sab' :'green',
#         'Sb' : '#00D0C9',
#         'Sbc' : '#0076C9',
#         'Sc_Scd_Sd' : 'blue',
#     } 
#             
#     SFRSD = {}
#     SFRSD['radius'] = {}
#     SFRSD['zones'] = {}
#     SFRSD['integrated'] = {}
#     ### ZONES ###
#     n_xm = []
#     n_ym = []
#     y = H.SFRSD_Ha__g
#     for iT, tSF in enumerate(H.tSF__T):
#         x = H.SFRSD__Tg[iT]
#         xm, ym = C.ma_mask_xyz(x, y, mask = ~maskRadiusOk__g)
#         n_xm.append(xm)
#         n_ym.append(ym)
#     SFRSD['zones']['x'] = n_xm
#     SFRSD['zones']['y'] = n_ym
#     ### RADIUS ###
#     n_xm = []
#     n_ym = []
#     y = H.aSFRSD_Ha__rg
#     for iT, tSF in enumerate(H.tSF__T):
#         x = H.aSFRSD__Trg[iT]
#         xm, ym = C.ma_mask_xyz(x, y, mask = ~maskRadiusOk__rg)
#         n_xm.append(xm)
#         n_ym.append(ym)
#     SFRSD['radius']['x'] = n_xm
#     SFRSD['radius']['y'] = n_ym
#     ### integrated ###
#     n_xm = []
#     n_ym = []
#     y = H.integrated_SFRSD_Ha__g
#     for iT, tSF in enumerate(H.tSF__T):
#         x = H.integrated_SFRSD__Tg[iT]
#         xm, ym = C.ma_mask_xyz(x, y)
#         n_xm.append(xm)
#         n_ym.append(ym)
#     SFRSD['integrated']['x'] = n_xm
#     SFRSD['integrated']['y'] = n_ym
#     
#     default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)
#     series_types = ['radius', 'zones', 'integrated' ] 
# 
#     if args.dryrun:
#         sys.exit('dryrun')
#     
#     SFRSD__mtypes = {}
#     for mt in mtype_bins:
#         key = 'list_' + mt
#         SFRSD__mtypes[key] = {}
#         SFRSD__mtypes[key]['radius'] = {}
#         SFRSD__mtypes[key]['zones'] = {}
#         SFRSD__mtypes[key]['integrated'] = {}
#         for k in series_types:
#             n_rs = []
#             n_R = []
#             print mt, k, maskmtype_ok[key][k].sum()
#             for iT, tSF in enumerate(H.tSF__T):
#                 x = np.ma.log10(SFRSD[k]['x'][iT])
#                 y = np.ma.log10(SFRSD[k]['y'][iT])
#                 # The radius mask is already inside SFRSD
#                 mask = ~(maskmtype_ok[key][k])
#                 xm, ym = C.ma_mask_xyz(x, y, mask = mask)
#                 rs = runstats(xm.compressed(), ym.compressed(), nBox = 10, **default_rs_kwargs)
#                 n_rs.append(rs)
#                 if k == 'zones':
#                     R = np.ma.masked_array(H.zone_dist_HLR__g, mask = xm.mask).compressed()
#                 elif k == 'radius':
#                     R = np.ma.masked_array(H.Rtoplot(), mask = xm.mask).compressed()
#                 else:
#                     R = None
#                 n_R.append(R)
#             SFRSD__mtypes[key][k]['R'] = n_R
#             SFRSD__mtypes[key][k]['rs'] = n_rs
# 
#     if args.dryrun:
#         sys.exit('dryrun')
#     
#     if len(args.output) != 0:
#         output = args.output + fnamesuffix
#     else:
#         output = 'SFRSD' + fnamesuffix
#     
#     default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
# 
#     tSF_txt = [ '32 Myr', '300 Myr', '1.5 Gyr', '14.2 Gyr' ]
#     with PdfPages(output) as pdf:
#         series_key = 'radius'
#         print output
#         ########################
#         ########################
#         ########################
#         rs_kwargs = default_rs_kwargs.copy()
#         sc_kwargs = default_sc_kwargs.copy()
#         NRows = 1
#         NCols = H.N_T
#         page_size_inches = (NCols * 3.5, NRows * 4)
#         grid_shape = (NRows, NCols)
#         f = plt.figure()
#         f.set_size_inches(page_size_inches)
#         #f.suptitle(txt_suptitle, fontsize = 10)
#         last_iT = H.N_T - 1
#         for iT, tSF in enumerate(H.tSF__T):
#             ax = plt.subplot2grid(grid_shape, loc = (0, iT))
#             ax.set_axis_on()
#             for i_mt, mt in enumerate(mtype_bins):
#                 k = 'list_' + mt
#                 v = SFRSD__mtypes[k][series_key]
#                 rs = v['rs'][iT]
#                 #ax.scatter(xm, ym, c = mtype_colors[k[5:]], **sc_kwargs)
#                 label = '%s' % k.replace('_', '-')[5:]
#                 ax.plot(rs.xS + 9, rs.yS + 9, c = mtype_colors[k[5:]], label = label)
#                 ax.set_xlabel(r'$\log\ \Sigma_{SFR}^\star (t_\star, R)\ [M_\odot Gyr^{-1} pc^{-2}]$')
#                 ax.set_ylabel('$\log\ \Sigma_{SFR}^{Neb} (R)\ [M_\odot Gyr^{-1} pc^{-2}]$')
#                 ax.grid(True)
#                 ax.set_xticks([ -1, 0, 1, 2, 3])
#                 ax.set_yticks([ -1, 0, 1, 2, 3])
#                 ax.set_xlim([-1, 3])
#                 ax.set_ylim([-1, 3])
#             if iT == last_iT:
#                 ax.legend(bbox_to_anchor = (1.6, 0.9), fontsize = 12, frameon = False)
#             ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
#             plot_text_ax(ax, txt = tSF_txt[iT], xpos = 0.02, ypos = 0.98, fontsize = 12, va = 'top', ha = 'left', c = 'k')
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         # cb = f.colorbar(sc)
#         # cb.set_label('R [HLR]')
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         f.subplots_adjust(bottom = 0.15, wspace = 0.3, right = 0.9)
#         #f.savefig('DGR_%s_%s' % (prop_key, output))
#         pdf.savefig(f)
#         plt.close(f)
# 
#         rs_kwargs = default_rs_kwargs.copy()
#         sc_kwargs = default_sc_kwargs.copy()
#         NRows = H.N_T
#         NCols = len(mtype_bins)
#         grid_shape = (NRows, NCols)
#         f = plt.figure()
#         f.set_size_inches((page_size_inches[0], NRows * 2.5))
#         #f.suptitle(txt_suptitle, fontsize = 10)
#         for iT, tSF in enumerate(H.tSF__T):
#             for i_mt, mt in enumerate(mtype_bins):
#                 ax = plt.subplot2grid(grid_shape, loc = (iT, i_mt))
#                 ax.set_axis_on()
#                 k = 'list_' + mt
#                 v = SFRSD__mtypes[k][series_key]
#                 rs = v['rs'][iT]
#                 label = '%s' % k.replace('_', '-')[5:]
#                 c = mtype_colors[k[5:]]
#                 ax.plot(rs.xS + 9, rs.yS + 9, c = mtype_colors[k[5:]], label = label)
#                 sc = ax.scatter(rs.x + 9, rs.y + 9, c = v['R'][iT], vmax = 3, vmin = 0, cmap = mpl.cm.spectral_r, **sc_kwargs)
#                 ax.plot(rs.xPrc[0] + 9, rs.yPrc[0] + 9, '--', c = c, label = None)
#                 ax.plot(rs.xPrc[1] + 9, rs.yPrc[1] + 9, '--', c = c, label = None)
#                 ax.plot(rs.xPrc[2] + 9, rs.yPrc[2] + 9, '--', c = c, label = None)
#                 ax.plot(rs.xPrc[3] + 9, rs.yPrc[3] + 9, '--', c = c, label = None)
#                 ax.grid(True)
#                 ax.set_xticks([ -1, 0, 1, 2, 3])
#                 ax.set_yticks([ -1, 0, 1, 2, 3])
#                 ax.set_xlim([-1, 3])
#                 ax.set_ylim([-1, 3])
#                 ax.plot(ax.get_xlim(), ax.get_xlim(), lw = .5, ls = "-", c = 'k')
#                 if i_mt > 0:
#                     plt.setp(ax.get_yticklabels(), visible = False)
#                 txt = '%s (%d)' % (label, len(rs.x))
#                 plot_text_ax(ax, txt = txt, xpos = 0.02, ypos = 0.98, fontsize = 12, va = 'top', ha = 'left', c = c)
#                 plot_text_ax(ax, txt = tSF_txt[iT], xpos = 0.98, ypos = 0.02, fontsize = 12, va = 'bottom', ha = 'right', c = 'k')
#                 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#                 #     if (iT == H.N_T - 1) and i_mt == 2:
#                 #         ax.set_xlabel(r'$\log\ \Sigma_{SFR}^\star (t_\star, R)\ [M_\odot Gyr^{-1} pc^{-2}]$', fontsize = 16)
#                 # else:
#                 #     if iT == 2:
#                 #         ax.set_ylabel('$\log\ \Sigma_{SFR}^{Neb} (R)\ [M_\odot Gyr^{-1} pc^{-2}]$', fontsize = 16)
#                 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         f.text(0.5, 0.04, r'$\log\ \Sigma_{SFR}^\star (t_\star, R)\ [M_\odot Gyr^{-1} pc^{-2}]$', ha='center', fontsize = 16)
#         f.text(0.04, 0.5, '$\log\ \Sigma_{SFR}^{Neb} (R)\ [M_\odot Gyr^{-1} pc^{-2}]$', va='center', rotation='vertical', fontsize = 16)
#         #f.tight_layout()        
#         f.subplots_adjust(left = 0.1, bottom = 0.1, hspace = 0.0, wspace = 0.0, right = 0.85)
#         cbar_ax = f.add_axes([0.88, 0.1, 0.02, 0.8])
#         cb = f.colorbar(sc, cax=cbar_ax)
#         cb.set_label('R [HLR]')
#         #f.savefig('DGR_%s_%s' % (prop_key, output))
#         pdf.savefig(f)
#         plt.close(f)
#         
# #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# #         series_key = 'zones'
# #         ########################
# #         ########################
# #         ########################
# #         rs_kwargs = default_rs_kwargs.copy()
# #         sc_kwargs = default_sc_kwargs.copy()
# #         NRows = 1
# #         NCols = H.N_T
# #         page_size_inches = (NCols * 3.5, NRows * 4)
# #         grid_shape = (NRows, NCols)
# #         f = plt.figure()
# #         f.set_size_inches(page_size_inches)
# #         f.suptitle(txt_suptitle, fontsize = 10)
# #         last_iT = H.N_T - 1
# #         for iT, tSF in enumerate(H.tSF__T):
# #             ax = plt.subplot2grid(grid_shape, loc = (0, iT))
# #             ax.set_axis_on()
# #             for i_mt, mt in enumerate(mtype_bins):
# #                 k = 'list_' + mt
# #                 v = SFRSD__mtypes[k][series_key]
# #                 rs = v['rs'][iT]
# #                 #ax.scatter(xm, ym, c = mtype_colors[k[5:]], **sc_kwargs)
# #                 label = '%s' % k.replace('_', '-')[5:]
# #                 ax.plot(rs.xS + 9, rs.yS + 9, c = mtype_colors[k[5:]], label = label)
# #                 ax.set_xlabel(r'$\log\ \Sigma_{SFR}^\star (t_\star)\ [M_\odot Gyr^{-1} pc^{-2}]$')
# #                 ax.set_ylabel('$\log\ \Sigma_{SFR}^{Neb}\ [M_\odot Gyr^{-1} pc^{-2}]$')
# #                 ax.grid(True)
# #                 ax.set_xticks([ -1, 0, 1, 2, 3])
# #                 ax.set_yticks([ -1, 0, 1, 2, 3])
# #                 ax.set_xlim([-1, 3])
# #                 ax.set_ylim([-1, 3])
# #             if iT == last_iT:
# #                 ax.legend(bbox_to_anchor = (1.6, 0.9), fontsize = 12, frameon = False)
# #             ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
# #             plot_text_ax(ax, txt = tSF_txt[iT], xpos = 0.02, ypos = 0.98, fontsize = 12, va = 'top', ha = 'left', c = 'k')
# #         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# #         # cb = f.colorbar(sc)
# #         # cb.set_label('R [HLR]')
# #         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# #         f.subplots_adjust(bottom = 0.15, wspace = 0.3, right = 0.9)
# #         #f.savefig('DGR_%s_%s' % (prop_key, output))
# #         pdf.savefig(f)
# #         plt.close(f)
# # 
# #         rs_kwargs = default_rs_kwargs.copy()
# #         sc_kwargs = default_sc_kwargs.copy()
# #         NRows = H.N_T
# #         NCols = len(mtype_bins)
# #         grid_shape = (NRows, NCols)
# #         f = plt.figure()
# #         f.set_size_inches((page_size_inches[0], NRows * 2.5))
# #         #f.suptitle(txt_suptitle, fontsize = 10)
# #         for iT, tSF in enumerate(H.tSF__T):
# #             for i_mt, mt in enumerate(mtype_bins):
# #                 ax = plt.subplot2grid(grid_shape, loc = (iT, i_mt))
# #                 ax.set_axis_on()
# #                 k = 'list_' + mt
# #                 v = SFRSD__mtypes[k][series_key]
# #                 rs = v['rs'][iT]
# #                 label = '%s' % k.replace('_', '-')[5:]
# #                 c = mtype_colors[k[5:]]
# #                 ax.plot(rs.xS + 9, rs.yS + 9, c = mtype_colors[k[5:]], label = label)
# #                 #ax.scatter(rs.x + 9, rs.y + 9, c = v['R'][iT], vmax = 3, vmin = 0, **sc_kwargs)
# #                 ax.plot(rs.xPrc[0] + 9, rs.yPrc[0] + 9, '--', c = c, label = None)
# #                 ax.plot(rs.xPrc[1] + 9, rs.yPrc[1] + 9, '--', c = c, label = None)
# #                 ax.plot(rs.xPrc[2] + 9, rs.yPrc[2] + 9, '--', c = c, label = None)
# #                 ax.plot(rs.xPrc[3] + 9, rs.yPrc[3] + 9, '--', c = c, label = None)
# #                 ax.grid(True)
# #                 ax.set_xticks([ -1, 0, 1, 2, 3])
# #                 ax.set_yticks([ -1, 0, 1, 2, 3])
# #                 ax.set_xlim([-1, 3])
# #                 ax.set_ylim([-1, 3])
# #                 ax.plot(ax.get_xlim(), ax.get_xlim(), lw = .5, ls = "-", c = 'k')
# #                 if i_mt > 0:
# #                     plt.setp(ax.get_yticklabels(), visible = False)
# #                 txt = '%s (%d)' % (label, len(rs.x))
# #                 plot_text_ax(ax, txt = txt, xpos = 0.02, ypos = 0.98, fontsize = 12, va = 'top', ha = 'left', c = c)
# #                 plot_text_ax(ax, txt = tSF_txt[iT], xpos = 0.98, ypos = 0.02, fontsize = 12, va = 'bottom', ha = 'right', c = 'k')
# #                 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# #                 #     if (iT == H.N_T - 1) and i_mt == 2:
# #                 #         ax.set_xlabel(r'$\log\ \Sigma_{SFR}^\star (t_\star, R)\ [M_\odot Gyr^{-1} pc^{-2}]$', fontsize = 16)
# #                 # else:
# #                 #     if iT == 2:
# #                 #         ax.set_ylabel('$\log\ \Sigma_{SFR}^{Neb} (R)\ [M_\odot Gyr^{-1} pc^{-2}]$', fontsize = 16)
# #                 #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# #         f.text(0.5, 0.04, r'$\log\ \Sigma_{SFR}^\star (t_\star)\ [M_\odot Gyr^{-1} pc^{-2}]$', ha='center', fontsize = 16)
# #         f.text(0.04, 0.5, '$\log\ \Sigma_{SFR}^{Neb}\ [M_\odot Gyr^{-1} pc^{-2}]$', va='center', rotation='vertical', fontsize = 16)
# #         #f.tight_layout()        
# #         f.subplots_adjust(left = 0.1, bottom = 0.1, hspace = 0.0, wspace = 0.0, right = 0.95)
# #         #f.savefig('DGR_%s_%s' % (prop_key, output))
# #         pdf.savefig(f)
# #         plt.close(f)    
# #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    
    
    
     
     
     
     
     
     
     
     
     
     
     
     
    #tSF__T              = H.tSF__T[0:20]
    tSF__T = H.tSF__T
    RRange = H.RRange
    RColor = H.RColor
    SFR__Tg = H.get_data_h5('SFR__Tg')
    SFR_Ha__g = H.get_data_h5('SFR_Ha__g')
    aSFRSD__Trg = H.get_data_h5('aSFRSD__Trg')
    aSFRSD_Ha__rg = H.get_data_h5('aSFRSD_Ha__rg')
    SFRSD__Tg = H.get_data_h5('SFRSD__Tg')
    SFRSD_Ha__g = H.get_data_h5('SFRSD_Ha__g')
    dist_zone__g = H.get_data_h5('dist_zone__g')
     
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    f.suptitle(txt_suptitle, fontsize = 11)
    filename = 'SFR_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFR__Tg[iT])
            y = np.ma.log10(SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-6, 0]
            yran = [-6, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    xlabel = r'$\log\ t_\star$ [yr]'  
    x = np.log10(tSF__T)
    plot_linreg_params(a, x, xlabel,
                       'a', 'SFR_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFR_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age%s' % fnamesuffix)
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age%s' % fnamesuffix, 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical') 
    filename = 'SFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)  
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__g & gals_slice__g))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'SFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'SFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'SFRSD_Rs_age%s' % fnamesuffix, 1., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')   
    filename = 'aSFRSD_linregress_report%s' % fnamesuffix
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(aSFRSD__Trg[iT] * 1e6)
            y = np.ma.log10(aSFRSD_Ha__rg * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-3.5, 1]
            yran = [-3.5, 1]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.4)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # ax = axArr[i, j]
            # age = tSF__T[iT]
            # n_mask = n_tot = 0
            # for iR, RUp in enumerate(RRange):
            #     if iR == 0:
            #         RMask = RbinCenter__r <= RUp
            #         txt = 'R <= %.1f HLR' % RUp
            #     else:
            #         RDown = RRange[iR - 1]
            #         RMask = (RbinCenter__r > RDown) & (RbinCenter__r <= RUp)
            #         txt = '%.1f < R <= %.1f HLR' % (RDown, RUp)
            #     c = RColor[iR] 
            #     x = np.ma.log10(aSFRSD__Trg[iT, RMask, :].flatten() * 1e6)
            #     y = np.ma.log10(aSFRSD_Ha__rg[RMask, :].flatten() * 1e6)
            #     mask = x.mask | y.mask
            #     xm = x[~mask]
            #     ym = y[~mask]
            #     n_mask += mask.sum()
            #     n_tot += len(x)
            #     if i == 0 and j == 0:
            #         if iR > 0:
            #             jump = 0.25 * (iR - 1.)
            #         pos_x = iR + jump
            #         #pos_y = pos_y_ini - (iR * pos_step)
            #         pos_y = 1.8
            #         textbox = dict(alpha = 0.)
            #         ax.text(pos_x, 1.1, txt,
            #                 fontsize = Rfontsize, color = c,
            #                 transform = ax.transAxes,
            #                 va = 'top', ha = 'left',
            #                 bbox = textbox)
            #     scat = ax.scatter(xm, ym, c = c, marker = 'o', s = 1., edgecolor = 'none', alpha = 1.)
            # C.debug_var(args.debug, masked = n_mask, not_masked = n_tot - n_mask, total = n_tot)
            # #print 'SigmaSFR x SigmaSFR_Ha Age: %.2f Myr: masked %d points of %d (Total: %d)' % (age / 1e6, n_mask, n_tot, n_tot - n_mask)
            # #ax.legend(loc = 'lower left', fontsize = 12, frameon = False)
            # age = tSF__T[iT]
            # xran = [-3.5, 1.]
            # yran = [-3.5, 1.]
            # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # # binsx = np.linspace(-4.5, 1., 51)
            # # binsy = np.linspace(min(ym),max(ym), 51)
            # # density_contour(xm, ym, binsx, binsy, ax=ax)
            # #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # x = np.ma.log10(aSFRSD__Trg[iT].flatten() * 1e6)
            # y = np.ma.log10(aSFRSD_Ha__rg.flatten() * 1e6)
            # xm, ym = C.ma_mask_xyz(x, y)
            # a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            # b2[iT] = (ym - xm).mean()
            # Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            # Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            # Y2 = xm + b2[iT]
            # Yrms = (ym - Y2).std()
            # ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            # if b2[iT] >= 0:
            #     txt = r'y = x + %.2f (rms:%.2f)' %  (b2[iT], Yrms)
            # else:
            #     txt = r'y = x - %.2f (rms:%.2f)' %  (-1. * b2[iT], Yrms)
            # C.debug_var(args.debug, y_hold_x = txt)
            # plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            # txt = '%.2f Myr' % (age / 1e6)
            # plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            # txt = '%.4f' % (Rs[iT])
            # plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            # txt = '%.4f' % (Rp[iT])
            # plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            # ax.set_xlim(xran)
            # ax.set_ylim(yran)
            # ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            # ax.xaxis.set_major_locator(MultipleLocator(1))
            # ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            # ax.yaxis.set_major_locator(MultipleLocator(1))
            # ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            # if i == NRows - 1 and j == 0:
            #     plt.setp(ax.get_xticklabels(), visible = True)
            #     plt.setp(ax.get_yticklabels(), visible = True)
            # C.debug_var(args.debug, age = age)
            # C.debug_var(args.debug, Rs = Rs[iT])
            # C.debug_var(args.debug, Rp = Rp[iT])
            # iT += 1
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
             
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'aSFRSD_linregress_slope_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep2_age%s' % fnamesuffix, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'aSFRSD_Rs_age%s' % fnamesuffix, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuffix, 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    pos_y_ini = 0.38
    pos_step = 0.09
    Rfontsize = 12
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27)
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star(R)}{\Sigma_{SFR}^\star(@1HLR)}$'
    ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}(R)}{\Sigma_{SFR}^{neb}(@1HLR)}$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    fnamesuftmp = '_norm%s' % fnamesuffix   
    filename = 'aSFRSD_linregress_report%s' % fnamesuftmp
    C.debug_var(args.debug, filename = filename)
    NAxes = len(f.axes)
    iT = 0
    jump = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)          
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j]
            age = tSF__T[iT]
            SFRSD_norm_GAL__g = (H.aSFRSD__Trg[iT][10, :] + H.aSFRSD__Trg[iT][9, :] / 2.)
            SFRSD_Ha_norm_GAL__g = (H.aSFRSD_Ha__rg[10, :] + H.aSFRSD_Ha__rg[9, :] / 2.)
            aSFRSD_norm__rg = H.aSFRSD__Trg[iT] / SFRSD_norm_GAL__g
            aSFRSD_Ha_norm__rg = H.aSFRSD_Ha__rg / SFRSD_Ha_norm_GAL__g
            xran = [-1.5, 1.5]
            yran = [-1.5, 1.5]
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # binsx = np.linspace(-4.5, 1., 51)
            # binsy = np.linspace(min(ym),max(ym), 51)
            # density_contour(xm, ym, binsx, binsy, ax=ax)
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            x = np.ma.log10(aSFRSD_norm__rg)
            y = np.ma.log10(aSFRSD_Ha_norm__rg)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~(maskRadiusOk__rg & gals_slice__rg))
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 0.3, edgecolor = 'none', alpha = 0.6)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'aSFRSD_linregress_slope_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'aSFRSD_linregress_intercep2_age%s' % fnamesuftmp, 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'aSFRSD_Rs_age%s' % fnamesuftmp, 1., 16) 
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age%s' % fnamesuftmp, 1., 16) 
    ############# Integrated #############
    ############# Integrated #############
    ############# Integrated #############
    integrated_SFR__Tg = H.integrated_SFR__Tg
    integrated_SFR_Ha__g = H.integrated_SFR_Ha__g
    integrated_SFRSD__Tg = H.integrated_SFRSD__Tg
    integrated_SFRSD_Ha__g = H.integrated_SFRSD_Ha__g
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(300)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ SFR_\star^{int}(t_\star)\ [M_\odot yr^{-1}]$' 
    ylabel = r'$\log\ SFR_{neb}^{int}\ [M_\odot yr^{-1}]$'
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFR_linregress_report.png'
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFR__Tg[iT])
            y = np.ma.log10(integrated_SFR_Ha__g)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFR x SFR_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())        
            xran = [-5, 2]
            yran = [-5, 2]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'integrated_SFR_linregress_slope_age.png', 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFR_linregress_intercep_age.png', 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFR_linregress_intercep2_age.png', 0., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rs', 'SFR_Rs_age.png', 1., 16) 
     
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # plot_linreg_params(sigma, x, xlabel, 
    #                    r'$\sigma$', 'SFR_linregress_sigma_age.png')
    # plot_linreg_params(r**2., x, xlabel, 
    #                    r'$r^2$', 'SFR_linregress_sqrcorr_age.png', 1., 16) 
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'integrated_SFR_Rs_age.png', 1., 16) 
    ###############################
    ###############################
    ###############################
    NRows = 4
    NCols = 5
    f, axArr = plt.subplots(NRows, NCols)
    f.set_dpi(96)
    f.set_size_inches(11.69, 8.27) 
    f.suptitle(txt_suptitle, fontsize = 11)
    plt.setp([a.get_xticklabels() for a in f.axes], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes], visible = False)
    xlabel = r'$\log\ \Sigma_{SFR}^\star(int, t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(int)\ [M_\odot yr^{-1} kpc^{-2}]$' 
    f.text(0.5, 0.04, xlabel, ha = 'center', va = 'center')
    f.text(0.06, 0.5, ylabel, ha = 'center', va = 'center', rotation = 'vertical')
    filename = 'integrated_SFRSD_linregress_report.png'
    C.debug_var(args.debug, filename = filename)   
    iT = 0
    a = np.ones_like(tSF__T)
    b = np.ones_like(tSF__T)
    b2 = np.ones_like(tSF__T)
    Rs = np.empty_like(tSF__T)
    Rp = np.empty_like(tSF__T)  
    for i in range(0, NRows):
        for j in range(0, NCols):
            ax = axArr[i, j] 
            x = np.ma.log10(integrated_SFRSD__Tg[iT] * 1e6)
            y = np.ma.log10(integrated_SFRSD_Ha__g * 1e6)
            xm, ym = C.ma_mask_xyz(x, y, mask = ~gals_slice__integr)
            age = tSF__T[iT]
            C.debug_var(args.debug, masked = xm.mask.sum(), not_masked = len(x) - xm.mask.sum(), total = len(x))
            #print 'integrated SFRSD x SFRSD_Ha Age: %.2f Myr: masked %d points of %d (total: %d)' % (age / 1e6, xm.mask.sum(), len(x), len(x) - xm.mask.sum())
            xran = [-5, 0]
            yran = [-5, 0]
            scat = ax.scatter(xm, ym, c = 'black', marker = 'o', s = 10, edgecolor = 'none', alpha = 0.8)
            a[iT], b[iT], sigma_a, sigma_b = plotOLSbisectorAxis(ax, xm, ym, pos_x = 0.96, y_pos = 0.05, fs = 8)
            b2[iT] = (ym - xm).mean()
            Rs[iT], _ = st.spearmanr(xm.compressed(), ym.compressed())
            Rp[iT], _ = st.pearsonr(xm.compressed(), ym.compressed())        
            Y2 = xm + b2[iT]
            Yrms = (ym - Y2).std()
            ax.plot(xm, Y2, c = 'b', ls = '--', lw = 0.5)
            if b2[iT] >= 0:
                txt = r'y = x + %.2f (rms:%.2f)' % (b2[iT], Yrms)
            else:
                txt = r'y = x - %.2f (rms:%.2f)' % (-1. * b2[iT], Yrms)
            C.debug_var(args.debug, y_hold_x = txt)
            plot_text_ax(ax, txt, 0.96, 0.09, 8, 'bottom', 'right', color = 'b')
            txt = '%.2f Myr' % (age / 1e6)
            plot_text_ax(ax, txt, 0.05, 0.96, 8, 'top', 'left')
            txt = '%.4f' % (Rs[iT])
            plot_text_ax(ax, txt, 0.05, 0.89, 8, 'top', 'left')
            txt = '%.4f' % (Rp[iT])
            plot_text_ax(ax, txt, 0.05, 0.84, 8, 'top', 'left')
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.plot(ax.get_xlim(), ax.get_xlim(), ls = "--", c = ".3")
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.xaxis.set_minor_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(1))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            if i == NRows - 1 and j == 0:
                plt.setp(ax.get_xticklabels(), visible = True)
                plt.setp(ax.get_yticklabels(), visible = True)
            C.debug_var(args.debug, age = age)
            C.debug_var(args.debug, Rs = Rs[iT])
            C.debug_var(args.debug, Rp = Rp[iT])
            iT += 1
    f.subplots_adjust(wspace = 0, hspace = 0, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig(filename)
    plt.close(f)
    x = np.log10(tSF__T)
    xlabel = r'$\log\ t_\star$ [yr]'
    plot_linreg_params(a, x, xlabel,
                       'a', 'integrated_SFRSD_linregress_slope_age.png', 1., 16) 
    plot_linreg_params(b, x, xlabel,
                       r'b [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFRSD_linregress_intercep_age.png', 0., 16)
    plot_linreg_params(b2, x, xlabel,
                       r'b2 [$M_\odot\ yr^{-1}\ kpc^{-2}$]', 'integrated_SFRSD_linregress_intercep2_age.png', 0., 16)
    plot_linreg_params(Rs, x, xlabel,
                       'Rs', 'integrated_SFRSD_Rs_age.png', 1., 16)
    plot_linreg_params(Rp, x, xlabel,
                       'Rp', 'SFR_Rs_age.png', 1., 16) 
  
