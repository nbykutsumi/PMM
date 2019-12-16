import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import os, sys
import glob
import h5py
import numpy as np
from datetime import datetime, timedelta
import socket

iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,30)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'

else:
    print 'check myhost'
    sys.exit()
#*******************************
a1dh = np.array([])
for DTime in lDTime:
    print DTime
    Year,Mon,Day = DTime.timetuple()[:3]
    #-- Read DPR-Ku  ----------------------------------------
    dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
    dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = dprDir + '/2B.GPM.DPRGMI.*.??????.V???.HDF5'
    lcmbPath = glob.glob(ssearch)

    if len(lcmbPath)==0:
        continue

    fig = plt.figure(figsize=(5,5))
    ax  = fig.add_axes([0.2,0.2,0.7,0.7])
    for cmbPath in lcmbPath:
        oid = int(cmbPath.split('.')[-3])
        with h5py.File(cmbPath, 'r') as h:
            a2prec  = h['NS/surfPrecipTotRate'][:]
            a3prwat = h['NS/precipTotWaterCont'][:]
            a3prwat = a3prwat[:,:,::-1]  # top to bottom --> bottom to top
 
        dprbaseDir = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06'
        dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)
        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No Ku file for oid=',oid
            continue

        with h5py.File(dprPath, 'r') as h:
            a2stop = h['NS/PRE/heightStormTop'][:]

        
        a3prtop = np.ones(a3prwat.shape) * np.arange(88)
        a2bintop = ma.masked_where(a3prwat<=0, a3prtop).argmax(axis=2)
        a2prtop = ma.masked_equal(a2bintop,0)*0.25 + 0.125
        a2stop  = ma.masked_less_equal(a2stop,0)*0.001

        a2mask  = ma.masked_less_equal(a2bintop, 0).mask
        a2mask  = a2mask + a2stop.mask
        a2mask  = a2mask + ma.masked_less_equal(a2prec, 0.1).mask

        a1prtop = ma.masked_where(a2mask, a2prtop).compressed()
        a1stop  = ma.masked_where(a2mask, a2stop ).compressed()

        a1dhTmp    = a1prtop - a1stop
        a1dh = np.concatenate([a1dh, a1dhTmp])



plt.hist( a1dh, bins=np.arange(-8+0.25,8,0.5) )
plt.savefig('/home/utsumi/temp/ret/hist.dh.org.png')
                                                          

#plt.xlim([0,22])
#ax.set_xlabel('Storm top by Ku product')
#ax.set_ylabel('Top height of condensed water')
#plt.savefig('/home/utsumi/temp/ret/scatter.stop.prtop.png')
#
#sys.exit()
