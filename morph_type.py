#!/usr/bin/python
import pyfits
import numpy as np

def morph_number(hubtyp, hubsubtyp):
    spirals_subtyp = {
        '0'     : 8,
        '0a'    : 9,
        'a'     : 10,
        'ab'    : 11, 
        'b'     : 12, 
        'bc'    : 13, 
        'c'     : 14, 
        'cd'    : 15, 
        'd'     : 16, 
        'dm'    : 17, 
        'm'     : 18,
    }
    
    if hubtyp == 'E':
        morphNumber = np.int(hubsubtyp)
    elif hubtyp == 'S':
        morphNumber = spirals_subtyp[hubsubtyp]
    elif hubtyp == 'I':
        morphNumber = 19
    else:
        morphNumber = -1
        
    return morphNumber

def get_morph(califaName = None, califaID = -1):
    morphFile = '/Users/lacerda/CALIFA/morph_eye_class.fits'
    hdu = pyfits.open(morphFile)
    
    gal = None
    
    if califaName:
        califaID = np.int(califaName[1:])
        
    if califaID >= 0:
        gal = hdu[1].data[np.where(hdu[1].data['CALIFAID'] == califaID)]
        
    del hdu
     
    return gal