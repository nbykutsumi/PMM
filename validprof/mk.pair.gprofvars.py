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

iDTime = datetime(2014,6,4)
eDTime = datetime(2014,6,4)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

thpr = 0.1
miss_out= -9999.
#lvar  = ['S1/Latitude', 'S1/Longitude'] 
lvar  = ['S1/Latitude', 'S1/Longitude']+ ['S1/qualityFlag','S1/surfaceTypeIndex'] 
#------------------------------------------------

def ave_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
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
        print gprofPath
 
        #-- Read and Save profile database of GPROF (Only once) ----
        if not os.path.exists(gprofPath):
            continue

        with h5py.File(gprofPath,'r') as h:

            '''    
            #---- hgtTopLayer -----
            [ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,
            6. ,  6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 11. , 12. ,
            13. , 14. , 15. , 16. , 17. , 18. ], dtype=float32)
            '''
    
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


        #-- Read GPROF --------------------------------------------
        with h5py.File(gprofPath,'r') as h: 
            a2sfcprecg = h['S1/surfacePrecipitation'][:,83:137+1]  # (Y,X)
 
        #-- Reshape GPROF --
        a1sfcprecg=a2sfcprecg.flatten()
    
        #-- Extract matching pixels from DPR array ---------------- 
        a1mask1 = ma.masked_less(a1x,0).mask
        a1mask2 = ma.masked_less(a1y,0).mask
        a1mask  = a1mask1 + a1mask2
    
        nyg,nxg = a2sfcprecg.shape
        #*** surface precip (DPR) *****
        a2sfcprecd = ma.masked_less(a2sfcprecd,0)
        a1sfcprecd = ave_9grids_2d(a2sfcprecd, a1y, a1x, miss=-9999.).filled(-9999.)

        a1sfcprecd[a1mask] = miss_out

        #-- Screen no precipitation cases -----------------------
        a1flag1 = ma.masked_greater(a1sfcprecd, thpr).mask
        a1flag2 = ma.masked_greater(a1sfcprecg, thpr).mask
        a1flag  = a1flag1 + a1flag2 
       
        # screen a1sfcprecd==-9999. --
        a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
        a1flag  = a1flag * a1flag3
        #print oid
        #print 'mask1',~a1mask1.sum(),a1mask1.shape
        #print 'mask2',~a1mask2.sum(),a1mask2.shape
        #print 'mask ',~a1mask.sum(), a1mask.shape
        #print 'a1flag1',a1flag1.sum(), a1flag1.shape
        #print 'a1flag2',a1flag2.sum(), a1flag2.shape
        #print 'a1flag3',a1flag3.sum(), a1flag3.shape
        #sys.exit()

        #-- Start var loop ----------------------------------------
        for var in lvar:
            print var
            varName = var.split('/')[-1]
            with h5py.File(gprofPath,'r') as h: 
                a2var = h[var][:,83:137+1]  # (Y,X)

            a1var = a2var.flatten() 
            a1var = a1var[a1flag] 
   
            outbaseDir = tankbaseDir + '/utsumi/PMM/validprof/pair/gprof'
            outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            util.mk_dir(outDir)
            np.save(outDir + '/%s.%06d.npy'%(varName, oid), a1var)
            print outDir
