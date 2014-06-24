#!/usr/bin/python
import pyfits
import numpy as np

def read_rgb_fits(fitsFile = None):
    flux__z = {}
    err_flux__z = {}
    fwhm__z = {}
    err_fwhm__z = {}
    ew__z = {}

    if fitsFile:
        hdu = pyfits.open(fitsFile)

        # Building the data arrays
        for i in range(len(hdu))[2:]:
            arrFlux = np.array(object = hdu[i].data['flux'][1:], dtype = hdu[i].data.dtype['flux'])
            mask = arrFlux > 0 
            flux__z[hdu[i].name] = np.ma.masked_array(arrFlux, mask = ~mask)
            
            arr = np.array(object = hdu[i].data['eflux'][1:], dtype = hdu[i].data.dtype['eflux'])
            err_flux__z[hdu[i].name] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[i].data['fwhm'][1:], dtype = hdu[i].data.dtype['fwhm'])
            fwhm__z[hdu[i].name] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[i].data['efwhm'][1:], dtype = hdu[i].data.dtype['efwhm'])
            err_fwhm__z[hdu[i].name] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[i].data['ew'][1:], dtype = hdu[i].data.dtype['EW'])
            ew__z[hdu[i].name] = np.ma.masked_array(arr, mask = ~mask)

        del hdu

    return flux__z, err_flux__z, fwhm__z, err_fwhm__z, ew__z
