from numpy import *
import myfunc.util as util
import os, sys
import glob
import h5py
import numpy as np
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

iDTime = datetime(2015,2,1)
eDTime = datetime(2015,2,28)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

thpr = 0.1
miss_out= -9999.
#lvar  = ['S1/Latitude', 'S1/Longitude'] 
#lvar  = ['S1/Latitude', 'S1/Longitude']+ ['S1/qualityFlag','S1/surfaceTypeIndex'] 
lvar  = ['S1/Latitude', 'S1/Longitude','S1/surfaceTypeIndex'] 


#------------------------------------------------
def ave_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''
    #-- Average 9 grids --
    nydpr,nxdpr = a2in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a2datTmp    = empty([9,len(a1y)], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a1datTmp= a2in[a1yTmp.filled(0),a1xTmp.filled(0)]

        a2datTmp[itmp,:] = a1datTmp

    #------------
    a1datTmp = ma.masked_equal(a2datTmp,miss).mean(axis=0)
    a1datTmp[a1dprmask] = miss

    return a1datTmp


#------------------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    
    gprofbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/2A/V05'
    gprofDir = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = gprofDir + '/2A.GPM.GMI.GPROF2017v1.*.??????.V05A.HDF5'
    lgprofPath = sort(glob.glob(ssearch))
    
    for gprofPath in lgprofPath: 
        oid = int(gprofPath.split('/')[-1].split('.')[-3])
        #print gprofPath
 
        #-- Read and Save profile database of GPROF (Only once) ----
        if not os.path.exists(gprofPath):
            continue

        #-- Read DPR -----------------------------------------------
        dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
        dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No DPR file for oid=',oid
            continue

        with h5py.File(dprPath, 'r') as h:
            a2sfcprecd = h['NS/surfPrecipTotRate'][:]
    
        #-- Read GMI-DPR matching index file ----------------------
        xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
        xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
        yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)

        if not os.path.exists(xPath):
            continue    
        a1x = np.load(xPath).flatten()
        a1y = np.load(yPath).flatten()

        #-- Read PMW precip ---
        pmwDir = tankbaseDir + '/utsumi/PMM/retepc/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
        pmwPath = pmwDir + '/nsurfMScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid,DB_MAXREC)
        a2sfcprecp = np.load(pmwPath)[:,83:137+1]

        #-- Reshape PMW --
        a1sfcprecp=a2sfcprecp.flatten()
    
        #-- Extract matching pixels from DPR array ---------------- 
        a1mask1 = ma.masked_less(a1x,0).mask
        a1mask2 = ma.masked_less(a1y,0).mask
        a1mask  = a1mask1 + a1mask2
    
        a1x = ma.masked_where(a1mask, a1x).filled(0)
        a1y = ma.masked_where(a1mask, a1y).filled(0)
    
        a2sfcprecd = ma.masked_less(a2sfcprecd,0)
        a1sfcprecd = ave_9grids_2d(a2sfcprecd, a1y, a1x, miss=-9999.).filled(-9999.)

        #a1sfcprecd[a1mask] = miss_out

        #-- Screen no precipitation cases -----------------------
        a1flag1 = ma.masked_greater(a1sfcprecd, thpr).mask
        a1flag2 = ma.masked_greater(a1sfcprecp, thpr).mask
        a1flag  = a1flag1 + a1flag2 
       
        # screen a1sfcprecd==-9999. --
        a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
        a1flag4 = ma.masked_not_equal(a1sfcprecp, -9999.).mask

        a1flag  = a1flag * a1flag3 * a1flag4
        #-- Start var loop ----------------------------------------
        for var in lvar:
            varName = var.split('/')[-1]
            with h5py.File(gprofPath,'r') as h: 
                a2var = h[var][:,83:137+1]  # (Y,X)

            a1var = a2var.flatten() 
            a1var = a1var[a1flag] 
   
            outbaseDir = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s'%(expr)
            outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            util.mk_dir(outDir)
            np.save(outDir + '/%s.%06d.npy'%(varName, oid), a1var)
            print outDir, oid
