import h5py
from numpy import *
import numpy as np

srcPath = './GPMCOR_DPR_1602200752_0925_011244_L2S_DD2_06A_ASS_test0000.h5'
with h5py.File(srcPath, 'r') as h:
    lkey1 = h.iterkeys()
    for key1 in lkey1:
        print ''
        print key1
        lkey2 = h[key1].iterkeys()
        for key2 in lkey2:
            print ''
            print '/%s/%s'%(key1,key2)

            try:
                lkey3 = h[key1+'/' + key2].iterkeys()
                for key3 in lkey3:
                    print h[key1+'/' + key2 + '/' + key3]
            except:
                print h[key1+'/' + key2]

