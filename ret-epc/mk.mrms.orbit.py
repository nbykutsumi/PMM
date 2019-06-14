from numpy import *
import numpy as np
import h5py
import myfunc.util as util
import glob
import gzip
from collections import deque


iYM = [2014,10]
eYM = [2014,10]
lYM = util.ret_lYM(iYM,eYM)

#*****************
# Read orbit list
#*****************
dorbit = {}
for (Year,Mon) in lYM:
    listDir = '/work/hk01/utsumi/PMM/US/obtlist'
    listPath= listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines=f.readlines(); f.close()
    for line in lines:
        line = map(int,line.split(','))
        yyyy,mm,dd,oidTmp,iy,ey = line
        dorbit[oidTmp] = [yyyy,mm,dd,iy,ey]


#*****************
# Start oid loop
#*****************
#loid = dorbit.keys()
loid = [3556]
for oid in loid:
    Year,Mon,Day,iy,ey = dorbit[oid]

    #*****************
    # Read GMI L1C
    #*****************
    gmiDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
    #gmiPath= gmiDir + '/1C.GPM.GMI.XCAL2016-C.20140505-S193910-E211142.001045.V05A.HDF5'
    ssearch = gmiDir + '/1C.GPM.GMI.XCAL2016-C.%04d%02d%02d-S??????-E??????.%06d.V???.HDF5'%(Year,Mon,Day,oid)
    gmiPath = glob.glob(ssearch)[0]
    with h5py.File(gmiPath,'r') as h:
        a2latgmi = h['S1/Latitude'][:]
        a2longmi = h['S1/Longitude'][:]
    
    iy,ey    = dorbit[oid][3:4+1]
    a2latgmi = a2latgmi[iy:ey+1]
    a2longmi = a2longmi[iy:ey+1]
    nygmi,nxgmi = a2latgmi.shape
    
    a1latgmi = a2latgmi.flatten()
    a1longmi = a2longmi.flatten()

    #*****************
    # Initialize output    
    #*****************
    a1out = ones(len(a1latgmi),float32)*(-9999.)

    #*****************
    # Make MRMS list
    #*****************
    mrmsDir = '/work/hk01/PMM/MRMS/match-GMI'
    #mrmsPath= mrmsDir + '/GMI-995-V04A.MRMS-20140502.152400.matched-130W_55W_20N_55N.extract.dat.gz'
    ssearch = mrmsDir + '/GMI-%d-V04A.MRMS-*.matched-130W_55W_20N_55N.extract.dat.gz'%(oid)
    lmrmsPath = glob.glob(ssearch)
    if len(lmrmsPath)==0:
        print 'No MRMS',Year,Mon,Day,oid

    #*****************
    # Read MRMS
    #*****************
    for mrmsPath in lmrmsPath:    
        with gzip.open(mrmsPath, 'rt') as f:
            lines = f.readlines()
            a1lon = deque([])
            a1lat = deque([])
            a1dat = deque([])
            for line in lines[1:]:
                line = map(float,line.split())
                a1lon.append(line[0])
                a1lat.append(line[1])
                a1dat.append(line[2])
        #*****************
        # Project MRMS over GMI track 
        #*****************
        for i in range(len(a1lon)):
            lat = a1lat[i]
            lon = a1lon[i]
            dat = a1dat[i]
            #lat = 46.291878  # test
            #lon = -129.8957  # test
        
            RADEARTH = 6371
            DTR      = 0.017453
            adist = RADEARTH*np.arccos(np.cos(DTR*a1longmi-DTR*lon)*np.cos(DTR*a1latgmi)*np.cos(DTR*lat) + np.sin(DTR*a1latgmi)*np.sin(DTR*lat))

            j = adist.argmin()
            a1out[j] = dat
        
    #*****************
    # Save
    #*****************
    a2out = a1out.reshape(nygmi,nxgmi).astype(float32)
    outDir= '/work/hk01/PMM/MRMS/match-GMI-orbit'
    util.mk_dir(outDir)
    outPath= outDir + '/GMI.MRMS.130W_55W_20N_55N.%04d%02d%02d.%06d.%05d-%05d.npy'%(Year,Mon,Day,oid,iy,ey)
    np.save(outPath,a2out)
    print outPath
   
     
