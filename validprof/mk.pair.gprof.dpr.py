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
iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

miss_out= -9999.
#------------------------------------------------
def ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale, lspecies=[0,2,3,4]):
    nh = 28
    ny,nx = a2tIndex.shape
    a4out = empty([len(lspecies),ny, nx, nh], dtype='float32')

    for i,species in enumerate(lspecies):
        a1profScale= a3profScale[:,:,species].flatten()
        a1profNum  = a3profNum[:,:,species].flatten()
        a1tIndex   = a2tIndex.flatten()

        #-- Handle non-precipitation pixels --
        a1flag = ma.masked_equal(a1profNum, 0).mask
        a1profNum[a1flag] = 1
        a1tIndex[a1flag] = 1

        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:, a1tIndex-1, species]
        a2prof[a1flag,:] = 0.0 
        a4out[i] = a2prof.reshape(ny,nx,nh)

    return a4out

def prof250mTo500m(a3prof, miss_out):
    ny,nx,nz = a3prof.shape
    #a3out = array([ma.masked_less(a3prof[:,:,i:i+2],0).mean(2) for i in range(0,nz,2)])
    a3out = ((ma.masked_less(a3prof[:,:,0::2],0) + ma.masked_less(a3prof[:,:,1::2],0) )*0.5).filled(miss_out)
    #a3out = ma.masked_less( concatenate([a3prof[:,:,0::2].reshape(ny,nx,-1,1), a3prof[:,:,1::2].reshape(ny,nx,-1,1)], axis=3), 0).mean(axis=3).filled(miss_out)

    return a3out


def gprofLayerconversion(a3prof): 
    ny,nx,nz = a3prof.shape
    a3outTop = zeros([ny,nx, (nz-20)*2],float32)
    a3outTop[:,:,0::2] = a3prof[:,:,20:]
    a3outTop[:,:,1::2] = a3prof[:,:,20:]
    return concatenate([a3prof[:,:,:20], a3outTop], axis=2)
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
            a4clusterProf= h['GprofDHeadr/clusterProfiles'][:]  # (profNumber, nlev, nT, nspecies) = (80, 28, 12, 5)
            a1hgtTopLayer= h['GprofDHeadr/hgtTopLayer'][:]  # (28,)
            species    = h['GprofDHeadr/speciesDescription'][:]
            species    = [''.join( map(chr, line) ) for line in species]
    
            '''
            #--- Species ---
            ['Rain Water Content\x00\x05\x00',
             'Cloud Water Content\x00!',
             'Ice Water Content\x00\x00\x00\x00',
             'Snow Water Content\x00  ',
             'Grauple/Hail Content\x00']
    
            #---- hgtTopLayer -----
            #[ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,
            6. ,  6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 11. , 12. ,
            13. , 14. , 15. , 16. , 17. , 18. ], dtype=float32)
            '''
    
        #-- Read GPROF ---------------------------------------------
        with h5py.File(gprofPath,'r') as h: 
            a2qFlag    = h['S1/qualityFlag'][:,83:137+1]  # (Y,X)
            a2sfcprecg = h['S1/surfacePrecipitation'][:,83:137+1]  # (Y,X)
            a2mlPrecip = h['S1/mostLikelyPrecipitation'][:,83:137+1] # (Y,X) 
            a2tIndex   = h['S1/profileTemp2mIndex'][:,83:137+1] # (Y,X)  zero=missing value? 
            a3profNum  = h['S1/profileNumber'][:,83:137+1,:] # (Y,X, nspecies)
            a3profScale= h['S1/profileScale'][:,83:137+1,:]  # (Y,X, nspecies)
    
        a4profg = ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale)
        a3profg = ma.masked_less(a4profg,0).sum(axis=0).filled(miss_out)
    
        a3profg = gprofLayerconversion(a3profg)
    
        #-- Reshape GPROF --
        a1sfcprecg=a2sfcprecg.flatten()
        a2profg = a3profg.reshape(-1,36)
    
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
            a3profd = h['NS/precipTotWaterCont'][:,:,16:]     # g/m3  (0-18 g/m3), 250m layers, missing=-9999.9,  Cut-off first 16 layers (22km-18km)
    
    
        a3profd = prof250mTo500m(a3profd, miss_out=miss_out)[:,:,::-1] # convert to 500m layers up to 18km (36-layers): From bottom to top. 
    
    
        #-- Read GMI-DPR matching index file ----------------------
        xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
        xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
        yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)


        if not os.path.exists(xPath):
            continue    
        a1x = np.load(xPath).flatten()
        a1y = np.load(yPath).flatten()
    
        #-- Extract matching pixels from DPR array ---------------- 
        a1mask1 = ma.masked_less(a1x,0).mask
        a1mask2 = ma.masked_less(a1y,0).mask
        a1mask  = a1mask1 + a1mask2
    
        a1x = ma.masked_where(a1mask, a1x).filled(0)
        a1y = ma.masked_where(a1mask, a1y).filled(0)
    
        nyg,nxg = a2sfcprecg.shape
        a1sfcprecd = a2sfcprecd[a1y,a1x]
        a2profd    = a3profd[a1y,a1x,:]
    
        a1sfcprecd[a1mask] = miss_out
        a2profd[a1mask,:] = miss_out
        #-- Screen no precipitation cases -----------------------
        a1flag1 = ma.masked_greater(a1sfcprecd, 0).mask
        a1flag2 = ma.masked_greater(a1sfcprecg, 0).mask
        a1flag  = a1flag1 + a1flag2 
       
        # screen a1sfcprecd==-9999. --
        a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
        a1flag  = a1flag * a1flag3
 
        a1sfcprecd = a1sfcprecd[a1flag] 
        a1sfcprecg = a1sfcprecg[a1flag] 
        a2profd    = a2profd[a1flag]   # Bottom to top
        a2profg    = a2profg[a1flag]   # Bottom to top
   
        outDir     = tankbaseDir + '/utsumi/validprof/pair.gprof/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        np.save(outDir + '/profpmw.%06d.npy'%(oid), a2profg)
        np.save(outDir + '/profrad.%06d.npy'%(oid), a2profd)
        np.save(outDir + '/precpmw.%06d.npy'%(oid), a1sfcprecg)    
        np.save(outDir + '/precrad.%06d.npy'%(oid), a1sfcprecd)
        print outDir 
