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
eDTime = datetime(2014,6,1)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

thpr   = 0.1
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

#***************************************

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    gprofbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/2A/V05'
    gprofDir = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = gprofDir + '/2A.GPM.GMI.GPROF2017v1.*.??????.V05A.HDF5'
    lgprofPath = sort(glob.glob(ssearch))
    if len(lgprofPath)==0:
        print 'no flies'
        print ssearch
    
    for gprofPath in lgprofPath: 
        oid = int(gprofPath.split('/')[-1].split('.')[-3])
        #print gprofPath
 
        #-- Read and Save profile database of GPROF (Only once) ----
        if not os.path.exists(gprofPath):
            print gprofPath
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


        sys.exit()
