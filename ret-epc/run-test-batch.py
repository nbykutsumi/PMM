import sys, os, shutil, socket
import myfunc.util as util
import glob
import subprocess
from datetime import datetime, timedelta
import numpy as np
import h5py
import shutil

iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,1)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)

batchsize = 3
#batchsize = 1

DB_MAXREC = 10000
DB_MINREC = 1000

#** Constants ******
#expr = 'glb.wprof.org'
#expr = 'glb.wprof.batch'
#expr = 'glb.wprof.rnr'
#expr = 'glb.nprof'
#expr = 'glb.v02.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
expr = 'test.batch2'

#prog = 'ret-myepc-29bins.py'
prog = 'ret-testbatch.py'

sensor  = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    gmibaseDir  = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/tank/utsumi/PMM/MATCH.GMI.V05A'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    outbaseDir = '/tank/utsumi/PMM/retepc/%s'%(expr) 
elif myhost =="well":
    #gmibaseDir  = '/media/disk2/share/data/PMM/NASA/GPM.GMI/1C/V05'
    gmibaseDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    #outbaseDir = '/media/disk2/share/PMM/retepc/%s'%(expr)
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s'%(expr)


#*******************
# Copy program
#*******************
copyDir = './progtemp'
util.mk_dir(copyDir)
progtime = datetime.now().strftime('%Y-%m-%d-%H:%M-%S-%f')
progcopy = copyDir + '/%s.%s'%(prog,progtime)
shutil.copy(prog, progcopy)
print progcopy
#*******************

icount = 0
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]    

    lgmiPathAll = np.sort(glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon,Day)))
    if len(lgmiPathAll)==0:
        continue

    #** Make batch list ***   
    llgmiPath = [lgmiPathAll[i*batchsize:(i+1)*batchsize] for i in range(int(len(lgmiPathAll)/batchsize) +1)]

    ##-- test --
    #llgmiPath = [['/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/07/01/1C.GPM.GMI.XCAL2016-C.20140701-S121951-E135224.001927.V05A.HDF5']]
    ##----------
    for lgmiPath in llgmiPath:
        icount = icount+1

        #** Loop in each batch ***
        a3tb1 = []
        a3tb2 = []
        a2lat = []
        a2lon = []
        a2x   = []
        a2y   = []
        a2t2m = []
        a2tqv = []
        a2elev= []
        loid  = []
        lny   = [] 
        for gmiPath in lgmiPath:
            oid = int(gmiPath.split('.')[-3])
            print 'oid=',oid
            #if oid <=5161: continue  # test

            loid.append(oid)
            #------------
            srcPath = glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
            s2xPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
            s2yPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
            t2mPath = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
            #tqvPathTmp = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid))[0]
            tqvPath = ''
            #elevPathTmp = glob.glob(matchbaseDir + '/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid))[0]
            elevPath = ''
        
            print srcPath
    
            #-- Append HDF file contents-----
            with h5py.File(srcPath,'r') as h: 
                a3tb1.append( h['/S1/Tc'][:])
                a3tb2.append( h['/S2/Tc'][:])
                a2lat.append( h['/S1/Latitude'][:])
                a2lon.append( h['/S1/Longitude'][:])

            #-- Append npy files -----
            a2x.append(np.load(s2xPath))
            a2y.append(np.load(s2yPath))
            a2t2m.append(np.load(t2mPath))
            if tqvPath !='':
                a2tqv.append(np.load(tqvPath))
            if elevPath !='':
                a2elev.append(np.load(elevPath))
            #-- Append NREC -----------------
            print s2xPath
            print np.load(s2xPath).shape[0], a3tb1[-1].shape[0]
            lny.append(a3tb1[-1].shape[0])
            #--------------------------------

        if len(loid)==0: continue

        #-- Concatenate data ----
        a3tb1 = np.concatenate(a3tb1, axis=0)
        a3tb2 = np.concatenate(a3tb2, axis=0)
        a2lat = np.concatenate(a2lat, axis=0)
        a2lon = np.concatenate(a2lon, axis=0)
        a2x   = np.concatenate(a2x, axis=0)
        a2y   = np.concatenate(a2y, axis=0)
        a2t2m = np.concatenate(a2t2m,axis=0)
        if tqvPath !='':
            a2tqv  = np.concatenate(a2tqv, axis=0) 
        if elevPath !='':
            a2elev = np.concatenate(a2elev,axis=0)

        #----------------------- 
        dargv = {}
        outDirTmp  = outbaseDir + '/%04d/%02d/%02d/temp'%(Year,Mon,Day)
        util.mk_dir(outDirTmp)
    

        #-- Batch file names ----
        oid = loid[0]
        srcPathTmp = outDirTmp + '/tb.%06d.HDF5'%(oid)
        s2xPathTmp = outDirTmp + '/Xpy.%06d.npy'%(oid)
        s2yPathTmp = outDirTmp + '/Ypy.%06d.npy'%(oid)
        t2mPathTmp = outDirTmp + '/t2m.%06d.npy'%(oid)
        if tqvPath !='':
            tqvPathTmp = outDirTmp + '/tqv.%06d.npy'%(oid)
        else:
            tqvPathTmp = ''
        if elevPath !='':
            elevPathTmp= outDirTmp + '/elev.%06d.npy'%(oid)
        else:
            elevPathTmp = ''
        #-- Save HDF file -------
        with h5py.File(srcPathTmp, 'w') as h:
            h.create_group('S1')
            h.create_group('S2')
            h.create_dataset('S1/Tc', data=a3tb1)
            h.create_dataset('S2/Tc', data=a3tb2)
            h.create_dataset('S1/Latitude', data=a2lat)
            h.create_dataset('S1/Longitude',data=a2lon)
            h.flush()

        #-- Save npy files ------
        np.save(s2xPathTmp, a2x)
        np.save(s2yPathTmp, a2y)
        np.save(t2mPathTmp, a2t2m)
        if tqvPath !='':
            np.save(tqvPathTmp, a2tqv)
        if elevPath !='':
            np.save(elevPathTmp, a2elev)

        #***** Set parameter dictionary *************
        dargv['sensor']=sensor
        dargv['coefDir']=coefDir
        dargv['dbDir']=dbDir
    
        #** Save program ******
        progDir = outbaseDir + '/prog'
        util.mk_dir(progDir)
        stime    = datetime.now().strftime('%Y.%m.%d_%H:%M:%S')
    
        iprog    = os.path.basename(__file__)
        oprog    = progDir + '/%s.%s'%(iprog, stime)
        shutil.copyfile(iprog, oprog)
        print oprog  
    
        iprog    = os.path.basename(prog)
        oprog    = progDir + '/%s.%s'%(iprog, stime)
    
        if icount==1:
            shutil.copyfile(iprog, oprog)
            print oprog  
    
        #**********************
        dargv['oid'] = oid
        dargv['clat'] = -9999
        dargv['clon'] = -9999
        dargv['dlatlon'] = -9999
        dargv['iscan'] = -9999
        dargv['escan'] = -9999
        dargv['dscan'] = -9999
        #------------
        dargv['NEM'] = 12
        dargv['NTBREG'] = 13
        dargv['NEM_USE'] = 3
        dargv['NPCHIST'] = 29
        dargv['NLEV_DPR'] = 50    # extract this number of layers
        dargv['NLEV_PRECIP'] = 50
        dargv['thwtmin'] = 0.1
        dargv['miss'] = -9999.
        dargv['miss_int32'] = np.int32(-9999)  
        #------------
        dargv['DB_MAXREC'] = DB_MAXREC
        dargv['DB_MINREC'] = DB_MINREC
        dargv['DB_USE_MINREC'] = 2
        dargv['NDB_EXPAND'] = 20
        dargv['DB_RAINFRAC'] = 0.0001 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
        dargv['MAX_T2M_DIFF'] = 10
        dargv['MAX_TQV_DIFF'] = 10
        dargv['MAX_RMA_0'] = 0.05 # Maximum Ratio of missing amount (0>=mm/h) acceptable for rain / no-rain classification # -9999. --> No screening
        dargv['flag_top_var'] = 0  # 0: No top-ranked vars. 1: Retrieve top-ranked vars.
        #dargv['flag_top_var'] = 1  # 0: No top-ranked vars. 1: Retrieve top-ranked vars.
        dargv['outDir'] = outDirTmp
    
        #------------
        oid = dargv['oid']
        dargv['srcPath'] = srcPathTmp
        dargv['s2xPath'] = s2xPathTmp
        dargv['s2yPath'] = s2yPathTmp
        dargv['t2mPath'] = t2mPathTmp
        dargv['tqvPath'] = tqvPathTmp
        dargv['elevPath']= elevPathTmp
    
        #-------------
        sargv = ['%s=%s'%(key, dargv[key]) for key in dargv.keys()]
        sargv = ' '.join(sargv)
        
        
        lcmd = ['python', progcopy, sargv]
        print lcmd
        subprocess.call(lcmd)

        #***************************
        # Split output files 
        #***************************
        iscan = dargv['iscan']
        escan = dargv['escan']
        DB_MAXREC = dargv['DB_MAXREC']

        
        ssearch = outDirTmp + '/*.nrec%d.npy'%(DB_MAXREC)
        ltmpPath = glob.glob(ssearch)
        for tmpPath in ltmpPath:
            varName = os.path.basename(tmpPath).split('.')[0]

            print varName
            atmp = np.load(tmpPath)
            y0 = 0
            for i in range(len(loid)):
                oid = loid[i]
                ny  = lny[i]
                aout= atmp[y0:y0+ny]
                y0  = y0+ny
                stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC)

                outDir  = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                outPath = outDir + '/%s.%s.npy'%(varName, stamp)
                np.save(outPath, aout)
                print outPath

        
        #***************************
        # Clean temporaty files
        #***************************
        ssearch = outDirTmp + '/*'
        ltmpPath = glob.glob(ssearch)
        for tmpPath in ltmpPath:
            print 'REMOVE    ',os.path.exists(tmpPath), tmpPath

            os.remove(tmpPath)
        os.rmdir(outDirTmp)
        print ''

os.remove(progcopy)
print 'remove progcopy',progcopy

