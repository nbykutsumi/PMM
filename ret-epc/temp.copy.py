# %%
import myfunc.util as util
import os, sys, shutil
from datetime import datetime, timedelta

iDTime = datetime(2014,11,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,12,9],[2014,12,10]]
useorblist = True

DB_MAXREC = 10000
DB_MINREC = 1000
#expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

tankDir = '/home/utsumi/mnt/lab_tank'
workbaseDir = '/home/utsumi/mnt/lab_work'
epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)

#*******************************
miss_out= -9999.

#*******************
# oid list
#*******************
orblistPath= tankDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.1000obts.txt'
#orblistPath= tankDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.2obts.txt'

f=open(orblistPath,'r'); lines=f.readlines(); f.close()
for gmiPath in lines:
    gmiPath = gmiPath.strip()
    oid = int(os.path.basename(gmiPath).split('.')[-3])
    Year,Mon,Day= map(int, os.path.dirname(gmiPath).split('/')[-3:])
    if [Year,Mon,Day] in lskipdates:
        continue
    DTimeTmp = datetime(Year,Mon,Day)
    if DTimeTmp < iDTime: continue
    if DTimeTmp > eDTime: continue

    srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf02.minrec1000.maxrec10000/%04d/%02d/%02d'%(Year,Mon,Day)
    ipath = srcdir + '/nsurfConvNScmb.%06d.y-9999--9999.nrec10000.npy'%(oid)

    srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.relsurf01.minrec1000.maxrec10000/%04d/%02d/%02d'%(Year,Mon,Day)
    opath = srcdir + '/nsurfConvNScmb.%06d.y-9999--9999.nrec10000.npy'%(oid)

    print ipath 
    print opath
    print ''
    shutil.copyfile(ipath,opath)
    if not os.path.exists(ipath):
        sys.exit()

#------------------------------------------------



# %%
