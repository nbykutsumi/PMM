import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l1_gmi as l1_gmi
from f_match_fov import *
import sys, os, glob
from datetime import datetime, timedelta
import numpy as np

#iDTime = datetime(2016,12,31)
iDTime = datetime(2017,7,2)
#eDTime = datetime(2018,1,1)
eDTime = datetime(2017,7,2)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)


#lDTime = [DTime for DTime in lDTime if not ((datetime(2017,9,26)<=DTime)&(datetime(2017,9,29))]
lDTime = [DTime for DTime in lDTime if not (datetime(2017,9,26)<=DTime)&(DTime<=datetime(2017,9,29))]


radar  = 'Ku'
gmi = l1_gmi.L1_GMI()
mwscan= 'S1'

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
wx  = ex0-ix0 +1
verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
#obaseDir   = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.IDX'%(fullverGMI, mwscan, ix0, ex0, radar, fullverDPR)
obaseDir   = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.S2.IDX'%(fullverGMI, mwscan, ix0, ex0)



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
        sys.exit()

    for srcPathGMI in lsrcPathGMI:
        #print srcPathGMI
        oid       = srcPathGMI.split('.')[-3] 

        if oid != '018986':continue   # test

        Lat0 = gmi.load_var_granule(srcPathGMI, '%s/Latitude'%('S1'))
        Lon0 = gmi.load_var_granule(srcPathGMI, '%s/Longitude'%('S1'))
        dtime0= gmi.load_dtime_granule(srcPathGMI, mwscan)

        Lat1 = gmi.load_var_granule(srcPathGMI, '%s/Latitude'%('S2'))
        Lon1 = gmi.load_var_granule(srcPathGMI, '%s/Longitude'%('S2'))
        dtime1= gmi.load_dtime_granule(srcPathGMI, mwscan)
        
        scori = gmi.load_var_granule(srcPathGMI, 'S1/SCstatus/SCorientation')
        
        #print 'S1.shape=',Lat0.shape
        #print 'S2.shape=',Lat1.shape
        
        LonSub0 = Lon0[:,ix0:ex0+1]
        LatSub0 = Lat0[:,ix0:ex0+1]
        ny,nx   = LonSub0.shape
        
        X1,X2,X3,X4,Y1,Y2,Y3,Y4 = f_match_fov.match_gmi_dpr(LonSub0.T, LatSub0.T, Lon1.T, Lat1.T)
        
        #print 'calc done'
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

        #y,x = 627, 111
        #print 'S2 1st:y,x=',Y1[y,x-83], X1[y,x-83]
        #print 'S2 2nd:y,x=',Y2[y,x-83], X2[y,x-83]
        #print 'S2 3rd:y,x=',Y3[y,x-83], X3[y,x-83]
        #print 'S2 4th:y,x=',Y4[y,x-83], X4[y,x-83]

        
        #outDir = obaseDir + '/temp/%04d/%02d/%02d'%(Year,Mon,Day)
        outDir = obaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        
        nameGMI= os.path.basename(srcPathGMI)
        #for i in [1]:
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
       



