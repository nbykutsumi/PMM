import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l2_dpr as l2_dpr
import myfunc.IO.GPM.l1_gmi as l1_gmi
from f_match_fov import *
import sys, os, glob
from datetime import datetime, timedelta
import numpy as np

iDTime = datetime(2014,9,16)
eDTime = datetime(2014,9,16)

dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
lDTime = [DTime for DTime in lDTime if not (datetime(2017,9,26)<=DTime)&(DTime<=datetime(2017,9,29))]

radar  = 'Ku'
gmi = l1_gmi.L1_GMI()
dpr = l2_dpr.L2_DPR()
mwscan= 'S1'

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
wx  = ex0-ix0 +1
verGMI = '05'
verDPR = '06'
subverGMI = 'A'
subverDPR = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
fullverDPR = '%s%s'%(verDPR,subverDPR)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
#baseDirDPR = '/work/hk01/PMM/NASA/GPM.Ku/2A/V%s'%(verDPR)
baseDirDPR = '/work/hk01/PMM/NASA/GPM.DPRGMI/2B/V%s'%(verDPR)
obaseDir   = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.IDX'%(fullverGMI, mwscan, ix0, ex0, radar, fullverDPR)


for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]   
    print DTime 
    
    srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
    #ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
    ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.HDF5'
    lsrcPathGMI = sort(glob.glob(ssearchGMI))

    if len(lsrcPathGMI)==0:
        print 'No GMI file',Year,Mon,Day
        print ssearchGMI
        continue

    for srcPathGMI in lsrcPathGMI:
        oid       = srcPathGMI.split('.')[-3] 

        # test --
        #loid = [3019,3296,3556,3609,3694,4186,4817,4832,5130]
        #if int(oid) not in loid: continue   # test

        #print srcPathGMI


        srcDirDPR = baseDirDPR + '/%04d/%02d/%02d'%(Year,Mon,Day)

        ##-- Ku case ---
        #ssearch   = srcDirDPR  + '/2A.GPM.Ku.V8-20180723.20170102-S072749-E090022.016175.V06A.HDF5'
        #ssearch   = srcDirDPR  + '/2A.GPM.Ku.*.%s.V06A.HDF5'%(oid)

        #-- CMB case ---
        ssearch   = srcDirDPR  + '/2B.GPM.DPRGMI.CORRA2018.20140501-S060650-E073923.000974.V06A.HDF5'
        ssearch   = srcDirDPR  + '/2B.GPM.DPRGMI.*.%s.V06A.HDF5'%(oid)

        #--------------


        lsrcPathDPR = sort(glob.glob(ssearch))
        if len(lsrcPathDPR)==0:
            print 'No DPR file for',Year,Mon,Day
            print ssearch
            sys.exit()
        elif len(lsrcPathDPR)>1:
            print 'too many files for',Year,Mon,Day
            print ssearch
            sys.exit()

        srcPathDPR = lsrcPathDPR[0] 
        
        Lat0 = gmi.load_var_granule(srcPathGMI, '%s/Latitude'%(mwscan))
        Lon0 = gmi.load_var_granule(srcPathGMI, '%s/Longitude'%(mwscan))
        dtime0= gmi.load_dtime_granule(srcPathGMI, mwscan)
        
        Lat1 = dpr.load_var_granule(srcPathDPR, 'NS/Latitude')
        Lon1 = dpr.load_var_granule(srcPathDPR, 'NS/Longitude')
        dtime1= dpr.load_dtime_granule(srcPathDPR, 'NS')
        
        scori = dpr.load_var_granule(srcPathDPR, 'NS/scanStatus/SCorientation')
        
        #print 'GMI.shape=',Lat0.shape
        #print 'GPR.shape=',Lat1.shape
        
        LonSub0 = Lon0[:,ix0:ex0+1]
        LatSub0 = Lat0[:,ix0:ex0+1]
        ny,nx   = LonSub0.shape
        print ny,nx

        if Lat1.shape[0]==0: continue

        X1,X2,X3,X4,Y1,Y2,Y3,Y4 = f_match_fov.match_gmi_dpr(LonSub0.T, LatSub0.T, Lon1.T, Lat1.T)
        #sys.exit()        
        
        print 'calc done'
        X1 = X1.T
        X2 = X2.T
        X3 = X3.T
        X4 = X4.T
        Y1 = Y1.T
        Y2 = Y2.T
        Y3 = Y3.T
        Y4 = Y4.T

        # xfort, yfort --> xpy, ypy
        X1 = (ma.masked_less(X1,0)-1).data
        X2 = (ma.masked_less(X2,0)-1).data
        X3 = (ma.masked_less(X3,0)-1).data
        X4 = (ma.masked_less(X4,0)-1).data

        Y1 = (ma.masked_less(Y1,0)-1).data
        Y2 = (ma.masked_less(Y2,0)-1).data
        Y3 = (ma.masked_less(Y3,0)-1).data
        Y4 = (ma.masked_less(Y4,0)-1).data

        
        outDir = obaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        
        nameGMI= os.path.basename(srcPathGMI)
        for i in [1,2,3,4]:
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
       



