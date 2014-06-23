#!/usr/bin/python
from astropy.io import fits
import numpy as np

def read_rbs_fits(fitsFile = None):
    flux__z = {}
    err_flux__z = {}
    fwhm__z = {}
    err_fwhm__z = {}
    ew__z = {}

    if fitsFile:
        hdu = fits.open(fitsFile)

        # Building the data arrays
        for i in range(len(hdu))[2:]:
            flux__z[hdu[i].name] = np.array(object = hdu[i].data['flux'], dtype = hdu[i].data.dtype['flux'])
            err_flux__z[hdu[i].name] = np.array(object = hdu[i].data['eflux'], dtype = hdu[i].data.dtype['eflux'])
            fwhm__z[hdu[i].name] = np.array(object = hdu[i].data['fwhm'], dtype = hdu[i].data.dtype['fwhm'])
            err_fwhm__z[hdu[i].name] = np.array(object = hdu[i].data['efwhm'], dtype = hdu[i].data.dtype['efwhm'])
            ew__z[hdu[i].name] = np.array(object = hdu[i].data['ew'], dtype = hdu[i].data.dtype['EW'])

        del hdu

    return flux__z, err_flux__z, fwhm__z, err_fwhm__z, ew__z