import os, sys
import subprocess, glob
from datetime import datetime, timedelta

#*******************
# Function
#*******************
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass
#------------
def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]
#*******************

iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,1)
dDTime = timedelta(days=1)
lDTime = ret_lDTime(iDTime,eDTime,dDTime)


DB_MAXREC = 10000
sensor = 'GMI'
gmibaseDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05'
matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
tankDir = '/home/utsumi/mnt/lab_tank'


packdir = '/home/utsumi/temp/test-ret'

##***********************
## Coefficient dir
##***********************
#idir = coefDir
#odir = packdir + '/EPCDB/EPC_COEF/GMI'
#cmd = 'rsync -avr %s %s'%(idir, odir)
#print cmd
#subprocess.call(cmd.split(' '))

##***********************
## Orbit-wise data
##***********************

lgmiPathAll = []
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    lgmiPathTmp = sorted(glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon,Day)))
    lgmiPathAll = lgmiPathAll + lgmiPathTmp

for gmiPath in lgmiPathAll[:1]:
    oid = int(gmiPath.split('.')[-3])
    Year,Mon,Day = map(int, os.path.dirname(gmiPath).split('/')[-3:])
    print 'oid=',oid
    #if oid <=2780: continue  # test

    #------------
    srcPath = glob.glob(gmibaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
    s2xPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
    s2yPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
    t2mPath = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
    rnrPath = glob.glob(tankDir + '/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000/%04d/%02d/%02d/nsurfNScmb.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid))[0]


    osrcPath = packdir + '/GMI.L1C/' + '/'.join(srcPath.split('/')[-4:])
    os2xPath = packdir + '/MATCH.GMI.V05A/' + '/'.join(s2xPath.split('/')[-5:])
    os2yPath = packdir + '/MATCH.GMI.V05A/' + '/'.join(s2yPath.split('/')[-5:])
    ot2mPath = packdir + '/MATCH.GMI.V05A/' + '/'.join(t2mPath.split('/')[-5:])
    ornrPath = packdir + '/RNR/' + '/'.join(rnrPath.split('/')[-4:])

    #for [ipath,opath] in [[srcPath,osrcPath],[s2xPath,os2xPath],[os2yPath,os2yPath],[t2mPath,ot2mPath],[rnrPath,ornrPath]]:
    for [ipath,opath] in [[s2xPath,os2xPath],[os2yPath,os2yPath],[t2mPath,ot2mPath],[rnrPath,ornrPath]]:
        mk_dir(os.path.dirname(opath))
        cmd = 'rsync -av %s %s'%(ipath, opath)
        lcmd= cmd.split(' ')
        subprocess.call(lcmd) 
        print 'copy'
        print opath


