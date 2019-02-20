import os, sys
import glob
import myfunc.util as util
from datetime import datetime, timedelta

rootDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'

#lvarNameFull = ['GMI.Latitude','GMI.Longitude','GMI.Tc','GMI.TcS2','GMI.epc-s1','GMI.epcid-s1','GMI.surfacePrecipitation','GMI.surfaceTypeIndex','Ku.V06A.flagAnvil','Ku.V06A.heightStormTop','Ku.V06A.heightZeroDeg','Ku.V06A.typePrecip']
#lvarNameFull = ['GMI.Tc']
#lvarNameFull = ['Ku.V06A.typePrecip']
lvarNameFull = ['Ku.V06A.flagAnvil','Ku.V06A.heightStormTop','Ku.V06A.heightZeroDeg','Ku.V06A.typePrecip']
#lvarNameFull = ['GMI.Latitude','GMI.Longitude','GMI.Tc','GMI.TcS2','GMI.epc-s1','GMI.epcid-s1','GMI.surfacePrecipitation','GMI.surfaceTypeIndex']
#lvarNameFull = ['GMI.Latitude','GMI.Longitude','GMI.Tc','GMI.TcS2','GMI.surfacePrecipitation','GMI.surfaceTypeIndex']

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,12,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)

for varNameFull in lvarNameFull:
    sensor  = varNameFull.split('.')[0]
    varName = varNameFull.split('.')[-1]

    srcPath0 = None
    gnum0 = 0
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        srcDir   = rootDir + '/S1.ABp103-117.%s/%04d/%02d/%02d'%(varNameFull,Year,Mon,Day)
        if sensor == 'Ku':
            ssearch  = srcDir + '/%s.1.*.npy'%(varName)
        elif sensor=='GMI':
            ssearch  = srcDir + '/%s.*.npy'%(varName)
        else:
            print 'check sensor',sensor

        #print ssearch
        lsrcPath = sorted(glob.glob(ssearch))
        #print lsrcPath 
        if len(lsrcPath)==0:
            print '-----------------------------'
            print varNameFull, DTime
            print 'No Directory or Files'
            print srcDir

        for srcPath in lsrcPath:
            fileName = os.path.basename(srcPath)
            gnum = int(fileName.split('.')[-2])
            #print fileName
            if gnum != gnum0+1:
                print '-----------------------------'
                print varNameFull
                print 'prev='
                print srcPath0
                print 'next='
                print srcPath
                print DTime, 'prevous=',gnum0,'next=',gnum

            gnum0 = gnum
            srcPath0 = srcPath
