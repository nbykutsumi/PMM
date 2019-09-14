import sys, os, shutil, socket
import myfunc.util as util
import glob
import subprocess
from datetime import datetime, timedelta
import numpy as np

iDTime = datetime(2014,7,1)
eDTime = datetime(2014,7,1)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)



#** Constants ******
#expr = 'glb.wprof.org'
expr = 'test'
#prog = 'ret-myepc-29bins.py'
prog = 'ret-test.py'
sensor  = 'GMI'

myhost = socket.gethostname()
if myhost =="shui":
    gmibaseDir  = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/tank/utsumi/PMM/MATCH.GMI.V05A'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    outbaseDir = '/tank/utsumi/PMM/retepc/%s'%(expr) 
elif myhost =="well":
    #gmibaseDir  = '/media/disk2/share/data/PMM/NASA/GPM.GMI/1C/V05'
    gmibaseDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    #outbaseDir = '/media/disk2/share/PMM/retepc/%s'%(expr)
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s'%(expr)


#*******************

icount = 0
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]    

    lgmiPath = np.sort(glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon,Day)))
    if len(lgmiPath)==0:
        continue

    for gmiPath in lgmiPath:
        icount = icount+1
        oid = int(gmiPath.split('.')[-3])
        print 'oid=',oid
        if oid != 1927: continue  # test

        dargv = {}
        outDir  = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)

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

        #------------
        #oid     = 3556
        #clat    = 34    # SE.US case. oid = 003556
        #clon    = -86   # 2014/10/14  05:42:03 UTC
        #
        #dlatlon = 3  # used to search the domain center
        #dscan   = 90
        ##dscan   = 55
        ##dscan   = 5
        
        dargv['oid'] = oid
        dargv['clat'] = 50.8
        dargv['clon'] = 28.1
        dargv['dlatlon'] = 20
        dargv['iscan'] = -9999
        dargv['escan'] = -9999
        dargv['dscan'] = 10

        #dargv['oid'] = oid
        #dargv['clat'] = -9999
        #dargv['clon'] = -9999
        #dargv['dlatlon'] = -9999
        #dargv['iscan'] = -9999
        #dargv['escan'] = -9999
        #dargv['dscan'] = -9999
        #------------
        dargv['NEM'] = 12
        dargv['NTBREG'] = 13
        dargv['NEM_USE'] = 3
        dargv['NPCHIST'] = 29
        dargv['NLEV_DPR'] = 50    # extract this number of layers
        dargv['NLEV_PRECIP'] = 50
        dargv['thwtmin'] = 0.1
        dargv['miss'] = -9999.
        #------------
        dargv['DB_MAXREC'] = 20000
        dargv['DB_MINREC'] = 5000
        #dargv['DB_USE_MINREC'] = 200
        dargv['DB_USE_MINREC'] = 2  # test
        dargv['NDB_EXPAND'] = 20
        #dargv['DB_RAINFRAC'] = 0.01 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
        dargv['DB_RAINFRAC'] = 0.0001 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
        dargv['MAX_T2M_DIFF'] = 10
        dargv['MAX_TQV_DIFF'] = 10
        
        dargv['outDir'] = outDir

        #print matchbaseDir + '/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid)   # test
        #------------
        oid = dargv['oid']
        dargv['srcPath'] = glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
        dargv['s2xPath'] = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
        dargv['s2yPath'] = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
        dargv['t2mPath'] = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
        dargv['tqvPath'] = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['tqvPath'] = ''
        dargv['elevPath']= glob.glob(matchbaseDir + '/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['elevPath']= ''

        #dargv['srcPath'] = glob.glob('/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
        #dargv['s2xPath'] = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['s2yPath'] = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['t2mPath']  = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['tqvPath']  = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid))[0]
        #dargv['elevPath']= glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid))[0]

        
        #-------------
        sargv = ['%s=%s'%(key, dargv[key]) for key in dargv.keys()]
        sargv = ' '.join(sargv)
        
        
        lcmd = ['python', prog, sargv]
        print lcmd
        subprocess.call(lcmd)
