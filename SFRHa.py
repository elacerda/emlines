#!/usr/bin/python
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# 
# Calculate the k for SFR(Halpha) = k . L(Halpha)
# 
#     Lacerda@Saco - 9/Jul/2014
#     
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
import numpy as np
from pystarlight.util.base import StarlightBase
from pystarlight.util.constants import L_sun, h, c, yr_sec

useTrapz = True

bases = [ 'Padova1994.chab', 'Padova1994.salp', 'Padova2000.chab', 'Padova2000.salp' ]

baseFile    = '/Users/lacerda/LOCAL/data/Base.bc03.h5'

for i, b in enumerate(bases):
    base        = StarlightBase(baseFile, b, hdf5 = True)
    
    cmInAA  = 1e-8      # cm / AA
    
    lambda_Ha   = 6562.8                    # Angstrom
    #max_yr      = base.ageBase[-1]
    max_yr      = 100e6
    mask        = base.l_ssp <= 912         # Angstrom
    mask_age    = base.ageBase <= max_yr
    
    f_ssp   = base.f_ssp[:,:,mask]
    l       = base.l_ssp[mask]
    
    if useTrapz:
        qh = np.trapz(y = f_ssp * l * cmInAA * L_sun / (h * c), 
                      x = l, axis = 2) # 1 / Msol
         
        Nh = np.trapz(y = qh[:, mask_age], x = base.ageBase[mask_age], axis = 1) * yr_sec 
         
        k_SFR = 2.226 * lambda_Ha * L_sun * yr_sec / (Nh * h * c)
    else:
        def create_dx(x):
            dx          = np.zeros_like(x)
            dx[1:]      = (x[1:] - x[:-1])/2.   # dl/2 from right neighbor
            dx[:-1]     += dx[1:]               # dl/2 from left neighbor
            dx[0]       = 2 * dx[0]
            dx[-1]      = 2 * dx[-1]
            return dx
        
        age         = base.ageBase[mask_age]    # years
         
        dl = create_dx(l)
        qh__Zt = (f_ssp * l * dl * cmInAA * L_sun / (h * c)).sum(axis = 2)
         
        d_age = create_dx(age)
        Nh = (qh__Zt[:, mask_age] * d_age).sum(axis = 1) * yr_sec
         
        k_SFR = 2.226 * lambda_Ha * L_sun * yr_sec / (Nh * h * c) # M_sun / yr
    
    print b + ':'
    for i, Z in enumerate(base.metBase):
        print '\tZ=%.4f N_H=%e k_SFR=%.2f Msun/yr' % (Z, Nh[i], k_SFR[i])
