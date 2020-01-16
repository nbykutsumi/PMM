import numpy as np
import pylab as pl
#import matplotlib.gridspec as gridspec
from glob import glob
import numpy.ma as ma
import sys,os, glob
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
import socket
import h5py

stopexpr = 'best01'
outexpr  = stopexpr+'cr'
#stopact  = 'LTQZ'
stopact  = 'HTQZ'

#iDTime = datetime(2015,1,1)
#eDTime = datetime(2015,5,31)
#iDTime = datetime(2017,1,1)
#eDTime = datetime(2017,12,31)
iDTime = datetime(2014,12,16)
eDTime = datetime(2015,5,31)


lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,11,25],[2014,12,9],[2014,12,10]]

lisurf = np.arange(1,14+1)
miss  = -9999.
#Year = 2017
#lDTime  = np.load('/home/utsumi/bin/PMM/stop/ldtime-%04d-test.npy'%(Year), allow_pickle=True)
#lDTime  = np.sort(lDTime)
##n = len(lDTime)
##lDTime = lDTime[:36]
##lDTime = lDTime[36:]

myhost = socket.gethostname()
if myhost =='shui':
    pmmbaseDir ='/tank/utsumi/PMM'
    hdfbaseDir ='/work/hk02/PMM/NASA'
    stopbaseDir= '/tank/utsumi/PMM/stop'

    figDir = '/home.rainbow/utsumi/public_html/tempfig/stop'
elif myhost =='well':
    pmmbaseDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM'
    hdfbaseDir ='/home/utsumi/mnt/lab_work/hk02/PMM/NASA'
    stopbaseDir='/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'

else:
    print 'check myhost',myhost



#******************************************
# Functions
#******************************************


#******************************************
#**************************************************************
# Read parameters for correction
#--------------------------------------------------------------
paraPath = '/home/utsumi/bin/PMM/stop/para-corr-lin.csv'
f=open(paraPath,'r'); lines=f.readlines(); f.close()
dslope = {}
disept = {}
for line in lines[1:]:
    line = line.split(',')
    isurf= int(line[0])
    dslope[isurf] = float(line[2])
    disept[isurf] = float(line[3])
  
#**************************************************************
# Main loop start
#******************************************
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    if [Year,Mon,Day] in lskipdates: continue

    #** Read HDF GMI Tb ****
    ibaseDir= stopbaseDir + '/orbit/%s-%s-ssn%s'%(stopexpr,stopact, 0)
    iDir = ibaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)


    ssearch = iDir + '/stop.??????.npy'

    lsrcPath= glob.glob(ssearch)
    if len(lsrcPath)==0: continue

    lsrcPath = np.sort(lsrcPath)

    for srcPath in lsrcPath:
        oid = int(srcPath.split('.')[-2])

        #******************************************
        # Original stop file
        #******************************************
        a2in = np.load(srcPath)
    
        #** Read HDF GPROF surface type ****
        gprofDir = hdfbaseDir + '/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
        gprofPath= glob.glob(gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid))[0]
    
        with h5py.File(gprofPath,'r') as h:
            a2surf = h['/S1/surfaceTypeIndex'][:]
       
  
        if a2surf.shape[0]==0:
            print 'Blank array for'
            print gprofPath
            continue 
        #******************************************
        a2out = a2in.copy() 
        #******************************************
        # Screen surf type
        #******************************************
        for isurf in lisurf:
            a2out = ((ma.masked_where(a2surf !=isurf, a2out) - disept[isurf])/dslope[isurf]).data

        #**********************************
        # Replace invalid
        #**********************************
        a2out = ma.masked_where(a2in==miss, a2out).filled(miss)
 
        #******************************************
        # Save
        #******************************************
        outbaseDir= stopbaseDir + '/orbit/%s-%s-ssn%s'%(outexpr,stopact, 0)
        outDir = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        outPath= outDir + '/stop.%06d.npy'%(oid)
        np.save(outPath, a2out)
        print outPath
