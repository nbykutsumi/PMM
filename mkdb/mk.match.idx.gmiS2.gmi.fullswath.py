from numpy import *
#import myfunc.util as util
from f_match_fov import *
import sys, os, glob
from datetime import datetime, timedelta
import numpy as np
import h5py


#iDTime = datetime(2016,4,18)
#eDTime = datetime(2016,4,18)
iDTime = datetime(2018,1,1)
eDTime = datetime(2018,1,1)

def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

dDTime = timedelta(days=1)
lDTime = ret_lDTime(iDTime, eDTime, dDTime)

#lDTime = [DTime for DTime in lDTime if not ((datetime(2017,9,26)<=DTime)&(datetime(2017,9,29))]
#lDTime = [DTime for DTime in lDTime if not (datetime(2017,9,26)<=DTime)&(DTime<=datetime(2017,9,29))]


ix0 = 0   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 220  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
wx  = ex0-ix0 +1
oid = -9999

baseDirGMI = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05'
pmwPath    = ''  # Dont't change.
obaseDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX'

dydxDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx'
dyPath0  = dydxDir + '/dy.000.npy'
dxPath0  = dydxDir + '/dx.000.npy'
dyPath180= dydxDir + '/dy.180.npy'
dxPath180= dydxDir + '/dx.180.npy'

#*******************************************************
# Check standard input
#-------------------------------------------------------
argvs = sys.argv
if len(argvs) ==1:
    print '*****************************'
    print 'No standard input'
    print '*****************************'
elif len(argvs) >2:
    print 'Too many standard input'
    print 'Python [prog] [parameter-string]'
    sys.exit()
else:
    largvs = argvs[1].split()
    dargv = {}
    for argvs in largvs:
        key,param = argvs.split('=')
        dargv[key] = param
    Year       = int(dargv['Year'])
    Mon        = int(dargv['Mon'] )
    Day        = int(dargv['Day'] )
    pmwPath    = dargv['pmwPath']   # L1C data path
    dyPath0    = dargv['dyPath0']
    dxPath0    = dargv['dxPath0']
    dyPath180  = dargv['dyPath180']
    dxPath180  = dargv['dyPath180']
    obaseDir   = dargv['obaseDir'] 


    lDTime = [datetime(Year,Mon,Day)]
#*******************************************************
# Functions
#-------------------------------------------------------
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass

#*******************************************************

dy = 4 # searching pixel number
dx = 4 # searching pixel number
thdist = 10  # km

#*******************************************************
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]   
    print DTime 

    if pmwPath =='':
        srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
        #ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
        ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.HDF5'
        lsrcPathGMI = sort(glob.glob(ssearchGMI))

        if len(lsrcPathGMI)==0:
            print 'No GMI file',Year,Mon,Day
            print ssearchGMI
            continue
    else:
        lsrcPathGMI = [pmwPath]

    for srcPathGMI in lsrcPathGMI:
        print srcPathGMI
        oid       = srcPathGMI.split('.')[-3] 

        #if oid != '001573':continue   # test

        with h5py.File(srcPathGMI, 'r') as h:
            Lat0 = h['S1/Latitude'][:]
            Lon0 = h['S1/Longitude'][:]
    
            Lat1 = h['S2/Latitude'][:]
            Lon1 = h['S2/Longitude'][:]
            
            a1scori = h['S1/SCstatus/SCorientation'][:]

        #-- Read dy and dx file --
        a2dy0    = np.load(dyPath0)
        a2dx0    = np.load(dxPath0)
        a2dy180  = np.load(dyPath180)
        a2dx180  = np.load(dxPath180)

        #- adjust swath length -
        nyorg    = Lon0.shape[0]
        nytmp    = a2dy0.shape[0]
        dny      = abs(nyorg-nytmp)
        if nytmp==nyorg:
            pass
        elif dny >= 5:
            print 'too large dny',dny
            sys.exit()

        elif nytmp <nyorg:
            a2dy0 = concatenate([a2dy0[:dny],a2dy0], axis=0)
            a2dx0 = concatenate([a2dx0[:dny],a2dx0], axis=0)
            a2dy180 = concatenate([a2dy180[:dny],a2dy180], axis=0)
            a2dx180 = concatenate([a2dx180[:dny],a2dx180], axis=0)

        elif nytmp >nyorg:
            a2dy0   = a2dy0[dny:]
            a2dx0   = a2dx0[dny:]
            a2dy180 = a2dy180[dny:]
            a2dx180 = a2dx180[dny:]
        else:
            print 'Something is wrong'
            sys.exit()
    
        #------------------------- 
        LonSub0 = Lon0[:,ix0:ex0+1]
        LatSub0 = Lat0[:,ix0:ex0+1]
        a2dy0   = a2dy0[:,ix0:ex0+1]
        a2dx0   = a2dx0[:,ix0:ex0+1]
        a2dy180 = a2dy180[:,ix0:ex0+1]
        a2dx180 = a2dx180[:,ix0:ex0+1]


        X1,X2,X3,X4,Y1,Y2,Y3,Y4 = f_match_fov.match_s1s2(LonSub0.T, LatSub0.T, Lon1.T, Lat1.T, a1scori, a2dx0.T, a2dy0.T, a2dx180.T, a2dy180.T, dx, dy, thdist)
        #X1,X2,X3,X4,Y1,Y2,Y3,Y4 = f_match_fov.match_gmi_dpr(LonSub0.T, LatSub0.T, Lon1.T, Lat1.T)
        
        print 'calc done'
        X1 = X1.T
        #X2 = X2.T
        #X3 = X3.T
        #X4 = X4.T
        Y1 = Y1.T
        #Y2 = Y2.T
        #Y3 = Y3.T
        #Y4 = Y4.T

        # xfort, yfort --> xpy, ypy
        X1 = (ma.masked_less(X1,0)-1).data
        #X2 = (ma.masked_less(X2,0)-1).data
        #X3 = (ma.masked_less(X3,0)-1).data
        #X4 = (ma.masked_less(X4,0)-1).data

        Y1 = (ma.masked_less(Y1,0)-1).data
        #Y2 = (ma.masked_less(Y2,0)-1).data
        #Y3 = (ma.masked_less(Y3,0)-1).data
        #Y4 = (ma.masked_less(Y4,0)-1).data
        
        outDir = obaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        #outDir = obaseDir + '/temp/%04d/%02d/%02d'%(Year,Mon,Day)
        mk_dir(outDir)
        
        nameGMI= os.path.basename(srcPathGMI)
        #for i in [1,2,3,4]:
        for i in [1]:
            if i==1:
                X=X1
                Y=Y1
            elif i==2:
                X=X2
                Y=Y2
            elif i==3:
                X=X3
                Y=Y3
            elif i==4:
                X=X4
                Y=Y4

            xPath= outDir + '/Xpy.%d.%s.npy'%(i,oid)
            yPath= outDir + '/Ypy.%d.%s.npy'%(i,oid)

            np.save(xPath, X.astype(int16))
            np.save(yPath, Y.astype(int16))
        print xPath
       



