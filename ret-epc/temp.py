# %%
import matplotlib
matplotlib.use('Agg')
%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import glob
import os
import glob

sensor='MHS'
sate = 'NOAA18'
lmrmspath = np.sort(glob.glob('/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match/%s.%s/????/??/??/mrms*'%(sensor, sate)))
epcbasedir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/us.v01'

a1mrms = []
a1good = []
a1epc  = []
for mrmspath in lmrmspath:
    mrmsdir = os.path.dirname(mrmspath)
    Year,Mon,Day = map(int, mrmspath.split('/')[-4:-1])
    oid = int(mrmspath.split('.')[-3])
    iscan, escan = map(int, mrmspath.split('.')[-2].split('-'))

    #-- goodfrac of mrms -----
    goodpath = glob.glob(mrmsdir + '/goodfrac.%06d.*'%(oid))[0]

    #-- EPC ------------------
    epcdir = epcbasedir + '/%s.%s/%04d/%02d/%02d'%(sensor,sate,Year,Mon,Day)
    epcpath= glob.glob(epcdir + '/nsurfNScmb.%06d.*.npy'%(oid))[0]


    a2mrms = np.load(mrmspath)
    a2good = np.load(goodpath)
    a2epc  = np.load(epcpath)

    a2mask1 = ma.masked_less(a2mrms, 0).mask
    a2mask2 = ma.masked_less(a2epc, 0).mask
    a2mask3 = ma.masked_less(a2good, 0.5).mask
    a2mask  = a2mask1 + a2mask2 + a2mask3

    a1mrms.extend( ma.masked_where(a2mask, a2mrms).compressed())
    a1epc .extend( ma.masked_where(a2mask, a2epc ).compressed())

plt.scatter(a1mrms, a1epc)
plt.ylim([0,10])
plt.xlim([0,10])
plt.show()
print len(a1mrms)
print len(a1epc)
#for mrmspath in lmrms:
#    oid = mrmspath.split('.')[-3]
#    srcdir = os.path.dirname(mrmspath)
#    goodpath = glob.glob(srcdir + '/goodfrac.%s.*'%(oid))[0]
#    a2mrms = np.load(mrmspath)
#    a2good = np.load(goodpath)
#    print oid,a2mrms.max(), a2good.max()
#

# %%
