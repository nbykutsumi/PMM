from numpy import *
import numpy as np
import myfunc.util as util
from datetime import datetime, timedelta 
import glob
import os, sys
import socket

myhost = socket.gethostname()
if myhost =="shui":
    tankDir = '/tank'
elif myhost =="well":
    tankDir = '/home/utsumi/mnt/lab_tank'


iDTime = datetime(2015,3,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

nz   = 50
vres = 250 # m
tmpbottom = -1000 # m
tmpnz     = 2*nz + 4  # 104
miss = -9999.

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    #baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s'%(expr)
    baseDir = tankDir + '/utsumi/PMM/retepc/%s'%(expr)
    #matchbaseDir = #'/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo' 
    matchbaseDir = tankDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo'    

    srcDir  = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    elevDir = matchbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
   
    ssearch = srcDir + '/prwatprofNS-rs.??????.y-9999--9999.nrec%05d.npy'%(DB_MAXREC)
    lsrcPath = glob.glob(ssearch)
    lsrcPath = np.sort(lsrcPath)

    for srcPath in lsrcPath:
        oid  = int(srcPath.split('.')[-4])

        a3in = np.load(srcPath)
        a2elev= np.load(elevDir + '/gtopo.%06d.npy'%(oid))

        a2surfbin = ((a2elev-tmpbottom) / vres).astype('int16')
        a2surfbin = tmpnz - 1 - a2surfbin  # top = 0, bottom = tmpnz-1
        #print a2surfbin.min(), a2surfbin.max(), a2elev.min()
        #print a3in.shape

        nl,nx = a3in.shape[:2]

        a3tmp = np.ones([nl,nx,tmpnz], float32)* miss

        X,Y = np.meshgrid( np.arange(nx), np.arange(nl) )
        a1x = X.flatten().astype('int32')
        a1y = Y.flatten().astype('int32')
        a1surfbin = a2surfbin.flatten()

        for iz in range(nz):
            a1k = a1surfbin - nz+1 + iz  # if iz=nz-1 --> a1k=a1surfbin
            a3tmp[a1y,a1x,a1k] = a3in[:,:,iz].flatten()

        a3out = a3tmp[:,:,nz:2*nz].astype('float32')

        ifname = os.path.basename(srcPath)
        ofname = 'prwatprofNS.' + '.'.join(ifname.split('.')[1:])
        outPath= srcDir + '/' + ofname
        np.save(outPath, a3out)

        print outPath
        print a3out.shape


