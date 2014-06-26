#!/usr/bin/python
import pyfits
import numpy as np
import sys
import os
import ast
import re

def read_rgb_fits(fitsFile = None, lines = None):
    __funcname__ = sys._getframe().f_code.co_name
    
    if not fitsFile:
        return False 
    
    califaID = re.search('K[0-9][0-9][0-9][0-9]', fitsFile).group(0)
    
    if os.path.isfile(fitsFile): 
        try:
            hdu = pyfits.open(fitsFile)
        except IOError:
            print __funcname__ + ': cannot read file ' + fitsFile
            return False
        
        N_zones = (hdu['6563'].header['NAXIS2']) - 1
        str_lines = hdu[0].header['FLINES']

        # Building the data arrays
        if lines:
            N_lines = len(lines)
            
            for line in lines:
                if str_lines.find(line) < 0:
                    print self.__name__ + ': ' + califaID + ': line ' + line + ' not found.'
                    N_lines = N_lines - 1
        else:
            lines = ast.literal_eval(str_lines)
            N_lines = len(hdu) - 2
        
        testData = hdu[lines[0]].data['flux'][1:] > 0
        
        if testData.astype(int).sum() == 0:
            print __funcname__ + ': ' + califaID + ': test data: all zeros'
            return False
        
        flux__lz        = np.ma.zeros((N_lines, N_zones))
        err_flux__lz    = np.ma.zeros((N_lines, N_zones))
        fwhm__lz        = np.ma.zeros((N_lines, N_zones))
        err_fwhm__lz    = np.ma.zeros((N_lines, N_zones))
        ew__lz          = np.ma.zeros((N_lines, N_zones))
        
        for i_l, line in enumerate(lines):
            arrFlux = np.array(object = hdu[line].data['flux'][1:], dtype = hdu[line].data.dtype['flux'])
            mask = arrFlux > 0
            
            flux__lz[i_l] = np.ma.masked_array(arrFlux, mask = ~mask)
            
            arr = np.array(object = hdu[line].data['eflux'][1:], dtype = hdu[line].data.dtype['eflux'])
            err_flux__lz[i_l] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[line].data['fwhm'][1:], dtype = hdu[line].data.dtype['fwhm'])
            fwhm__lz[i_l] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[line].data['efwhm'][1:], dtype = hdu[line].data.dtype['efwhm'])
            err_fwhm__lz[i_l] = np.ma.masked_array(arr, mask = ~mask)
            
            arr = np.array(object = hdu[line].data['ew'][1:], dtype = hdu[line].data.dtype['EW'])
            ew__lz[i_l] = np.ma.masked_array(arr, mask = ~mask)
            
        del hdu
    else:
        print __funcname__ + ': ' + fitsFile + ': file does not exist'
        return False

    return flux__lz, err_flux__lz, fwhm__lz, err_fwhm__lz, ew__lz, lines
