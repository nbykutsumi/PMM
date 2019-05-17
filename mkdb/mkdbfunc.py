import numpy as np
from numpy import *
import sys

def mk_nonlin_comb(a3tb):
    ny,nx,ntb = a3tb.shape 
    dtype = a3tb.dtype
    ncomb = ntb + (ntb+1)*ntb/2 + (ntb-1)
    a3out = np.zeros([ny,nx,ncomb],dtype=dtype)

    #-- Tb :ntb ----
    a3out[:,:,:ntb] = a3tb
    kt = ntb-1

    #-- Tbi * Tbj : (ntb+1)*(ntb)/2 --- 
    for i in range(ntb):
        for j in range(i,ntb):
            kt = kt + 1
            a2x = a3tb[:,:,i] * a3tb[:,:,j]
            a3out[:,:,kt] = a2x

    #-- (Tbi-Tbi+1)/(Tbi+Tbi+1) : (ntb-1) ---
    for i in range(ntb-1):
        kt  = kt + 1
        a2x = (a3tb[:,:,i]-a3tb[:,:,i+1])/(a3tb[:,:,i]+a3tb[:,:,i+1])
        a3out[:,:,kt] = a2x

    return a3out
