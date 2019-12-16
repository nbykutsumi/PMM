import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket

#*******************************
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
iYM  = [2014,6]
eYM  = [2014,6]
lYM  = util.ret_lYM(iYM,eYM)
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)


a1out = np.array([])
for Year,Mon in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
    
    for DTime in lDTime:
        print DTime
        Year,Mon,Day = DTime.timetuple()[:3]
        if [Year,Mon,Day] in lskipdates:
            continue
        srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)

        ssearch = srcDir + '/%s.??????.npy'%('profrad')
        lprofPath = np.sort(glob.glob(ssearch))
    
        if len(lprofPath)==0:
            print 'No files'
            print ssearch
            sys.exit() 
        for profPath in lprofPath:
            oid = int(profPath.split('.')[-2])    
            #a1ph   = np.load(srcDir + '/stoprad.%06d.npy'%(oid)) 
            a1ph   = np.load(srcDir + '/heightStormToprad.%06d.npy'%(oid)) 
            #a1freez= np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid)) 
            a2prof = np.load(srcDir + '/profrad.%06d.npy'%(oid)) # bottom to top

            a2bin  = np.ones(a2prof.shape).astype(int16) * np.arange(36).astype(int16)
            a2bin  = ma.masked_where(a2prof <=0.01, a2bin).argmax(axis=1)
            a1wh   = ma.masked_less_equal(a2bin,0) * 0.5 + 0.25
            a1ph   = ma.masked_less(a1ph,0) * 0.001

            a1dh   = (ma.masked_equal(a1wh,0) - ma.masked_equal(a1ph,0)).flatten()
            a1out = np.concatenate([a1out, a1dh])

            if a1dh.shape[0]==0:
                print profPath
                continue
            print a1dh.max()


plt.hist( a1out, bins=np.arange(-8+0.25,8,0.5) )
plt.savefig('/home/utsumi/temp/ret/hist.dh.png')
