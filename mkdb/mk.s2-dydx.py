from numpy import *
import h5py
from datetime import datetime, timedelta
import glob
from f_match_fov import *
import myfunc.util as util
import os, sys
import numpy as np

lori = [0,180]
#lori = [180]

drev = {0:'016608',180:'016173'}  # ori=0: 016608: 2017-1-30,  ori=180: 016173: 2017-1-2
''' both orbits have swath length of 2946 '''
dictDTime= {0:datetime(2017,1,30), 180:datetime(2017,1,2)}


nx = 221
verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
Year = 2017
Mon  = 1

dx = 35
dy = 20
dycurve= 20
thdist = 20
miss   = -9999

#********************************************
for ori in lori:
    rev   = drev[ori]
    DTime = dictDTime[ori]
    Year,Mon,Day = DTime.timetuple()[:3]

    baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
    srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
    #ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
    ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.%s.*.HDF5'%(rev)
    print ssearchGMI
    lsrcPathGMI = glob.glob(ssearchGMI)
    srcPathGMI  = lsrcPathGMI[0]
    print srcPathGMI

    #-- Read HDF file --
    with h5py.File(srcPathGMI) as h:
        a2lat1 = h['S1/Latitude'][:]
        a2lon1 = h['S1/Longitude'][:]
        a2lat2 = h['S2/Latitude'][:]
        a2lon2 = h['S2/Longitude'][:]

        a1ori  = h['S1/SCstatus/SCorientation'][:]
    #wx =20
    #cx = 60
    #a2lat1 = a2lat1[:,cx-wx:cx+wx+1]
    #a2lon1 = a2lon1[:,cx-wx:cx+wx+1]
    #a2lat2 = a2lat2[:,cx-wx:cx+wx+1] 
    #a2lon2 = a2lon2[:,cx-wx:cx+wx+1]

    #print a2lat1.shape
    #print a1ori.min(), a1ori.max()
    #sys.exit()

    X1,X2,X3,X4,Y1,Y2,Y3,Y4 = f_match_fov.match_gmi_gen(a2lon1.T, a2lat1.T, a2lon2.T, a2lat2.T, dx, dy, dycurve, thdist)

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


    obaseDir= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
    outDir = obaseDir + '/S1.ABp000-220.GMI.S2.dydx'
    util.mk_dir(outDir)

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


        ny,nx = X.shape
        a2xorg, a2yorg = np.meshgrid(range(nx), range(ny))

        a2dx = (ma.masked_less(X,0) - a2xorg).filled(miss)
        a2dy = (ma.masked_less(Y,0) - a2yorg).filled(miss)

        xPath= outDir + '/dx.%03d.npy'%(ori)
        yPath= outDir + '/dy.%03d.npy'%(ori)

        np.save(xPath, a2dx.astype(int16))
        np.save(yPath, a2dy.astype(int16))
        print xPath



