import numpy as np
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import random
import socket
import scipy.stats as stats
import pickle
import string

season = 'JJADJF'
#season = 6
lrettype = ['epc']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#metric= 'rmse'
#metric= 'cc'
#metric = 'dwatNorm'
#metric = 'dconvfrac'
#metric = 'dstop'
nsample = 1
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    epcbaseDir  = '/tank/utsumi/PMM/retepc'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'

else:
    print 'check myhost'
    sys.exit()

nz    = 25

#** Functions ***************************
def ret_lym(season):
    if season=='JJADJF':
        lYM = util.ret_lYM([2014,6],[2014,8]) + util.ret_lYM([2014,12],[2015,2])

    elif season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
    elif type(season)==int:
        if season <6:
            lYM = [[2015,season]]
        else:
            lYM = [[2014,season]]
    elif season =='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    elif season =='DJ':
        lYM = util.ret_lYM([2014,12],[2015,1])

    else:
        print 'check season',season
        sys.exit()
    return lYM

#*****************************************
lYM = ret_lym(season)
for rettype in lrettype:

    #-- Initialize ------
    aq = np.array([])
    #--------------------


    for Year,Mon in lYM:
    
        eDay = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,1)
        eDTime = datetime(Year,Mon,eDay)
        dDTime = timedelta(days=1)
        lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
        #lDTime = lDTime[29:]  # test 
        for DTime in lDTime:
            Year,Mon,Day = DTime.timetuple()[:3]
    
            pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)

            ssearch= pairDir + '/Latitude.*.npy'
            llatPath= sort(glob.glob(ssearch))
            #-- random sampling -----
            if len(llatPath)>=nsample:
                llatPath = random.sample(llatPath, nsample)
            else:
                pass
            #------------------------
            
            for latPath in llatPath:
                print latPath
                oid = int(latPath.split('.')[-2])

                a1precradTmp= np.load(pairDir+ '/precrad.%06d.npy'%(oid))
                a2profradTmp= np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                a1stopradTmp=np.load(pairDir + '/stoprad.%06d.npy'%(oid))
                a1elevTmp   = np.load(pairDir + '/surfaceElevationrad.%06d.npy'%(oid))

                #--- Screen missing surface precipitation ---
                a1flagP = ma.masked_greater_equal(a1precradTmp,0).mask

                #--- Screen high mountains ------------------
                a1flagE  = ma.masked_less(a1elevTmp, 1000).mask

                #--- Screen storm top ----
                a1flagStop= ma.masked_inside(a1stopradTmp,2000+0.01,12500-0.01).mask
                ##---------------------
                a1flag = a1flagP * a1flagE * a1flagStop

                a1precradTmp = a1precradTmp[a1flag]
                a2profradTmp = a2profradTmp[a1flag,:] 
                a1stopradTmp = a1stopradTmp[a1flag]

                #--- Mask missing data --
                a2profradTmp = ma.masked_less(a2profradTmp,0)
                if a1precradTmp.shape[0]==0:
                    continue
    

                #--- Find bin number for storm top ----
                a1idx = np.arange(a2profradTmp.shape[0]).astype('int16')
                a1stopbin = ((a1stopradTmp-500*4) / 500).astype('int16')
                a1stopq   = a2profradTmp[a1idx, a1stopbin]
               
                aq = np.concatenate([aq, a1stopq])

                #sys.exit()
    #-----------------
    print 'mean=',aq.mean()
    print 'std=', aq.std()
