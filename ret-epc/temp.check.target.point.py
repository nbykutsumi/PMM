from numpy import *
import numpy as np
import os, sys
import myfunc.util as util
import socket, glob
import EPCDB
import h5py
import epcfunc

NEM = 12
DB_MAXREC = 10000
sensor = 'GMI'
myhost = socket.gethostname()
if myhost == 'shui':
    #srcDir = '/home/utsumi/temp/out'
    #srcDir = '/home/utsumi/temp/out/my'
    srcbaseDir = '/tank/utsumi/PMM/retepc'
    tankbaseDir  = '/tank'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    workbaseDir  = '/work'
    mrmsDir  = '/work/hk01/PMM/MRMS/match-GMI-orbit'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    #srcbaseDir = '/media/disk2/share/PMM/retepc/glb.wprof'
    #srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.wprof'
    srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    workbaseDir  = '/home/utsumi/mnt/lab_work'
    mrmsDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/MRMS/match-GMI-orbit'
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    #figDir   = '/home/utsumi/temp/ret'
    figDir   = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/ret'

else:
    print 'check hostname',myhost
    sys.exit()
#**************************************************************
def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return array(lout)

#**************************************************************
# Read parameters
#--------------------------------------------------------------
#-- Read PC coefficient file --
coefPath = coefDir + '/coef_pc.txt'
a2coef   = read_table(coefPath)
a2coef   = a2coef[:,1:]
#-- Read EPC range files --
rangePath = coefDir + '/PC_MIN_MAX_29.txt'
a2pc_edge = read_table(rangePath)

# expand the lowest and highest ranges
a2pc_edge[:,0]   = a2pc_edge[:,0] - 1.e6
a2pc_edge[:,-1]  = a2pc_edge[:,-1]+ 1.e6


#-- Read PC ave and std file --
pcavePath  = coefDir + '/ave_pc.txt'
#pcavePath  = '/home/utsumi/bin/ENSPR/ave_pc.txt'
a2pc_avestd= read_table(pcavePath)
a1pc_ave   = a2pc_avestd[:,1]
a1pc_std   = a2pc_avestd[:,2]



db = EPCDB.EPCDB()

oid  = 1686
Year,Mon,Day = 2014,6,16
iy, ey = 1746, 1766
clat    = 52.0
clon    = 270.0 -360.
#DB_MAXREC = 20000
DB_MAXREC = 10000
#DB_MINREC = 5000
DB_MINREC = 1000

y,x = 9,92
expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
srcDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)

idxdbPath = srcDir + '/top-idxdbNS.%06d.y%04d-%04d.nrec%d.npy'%(oid, iy, ey, DB_MAXREC)
irecPath  = srcDir + '/top-irecNS.%06d.y%04d-%04d.nrec%d.npy'%(oid, iy, ey, DB_MAXREC)
#precnsPath = srcDir + '/nsurfNS.%06d.y%04d-%04d.nrec%d.npy'%(oid, iy, ey, DB_MAXREC)
a2idxdb = np.load(idxdbPath)
a2irec  = np.load(irecPath)
#a2precns= np.load(precnsPath)


#********************************************
#-- Read Obs Tb data --
tbDir   = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)

ssearch = tbDir + '/1C.GPM.GMI.XCAL*.%06d.V05A.HDF5'%(oid)
tbPath  = glob.glob(ssearch)[0]


with h5py.File(tbPath, 'r') as h:
    a3tb1    = h['/S1/Tc'][:]
    a3tb2org = h['/S2/Tc'][:]
    #a2lat    = h['/S1/Latitude'][:]
    #a2lon    = h['/S1/Longitude'][:]
    #a2lat2= h['/S2/Latitude'][:]
    #a2lon2= h['/S2/Longitude'][:]

#********************************************
#-- Matchup and Joint S1 and S2 Tb --
s2xPath= tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid)
s2yPath= tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid)

a1x2  = np.load(s2xPath).flatten()
a1y2  = np.load(s2yPath).flatten()

a1mask= ma.masked_less(a1x2,0).mask
a1x2  = ma.masked_less(a1x2,0).filled(0)
a1y2  = ma.masked_less(a1y2,0).filled(0)

nytmp, nxtmp, ztmp = a3tb2org.shape
a2tb2 = a3tb2org[a1y2, a1x2]
a2tb2[a1mask] = -9999.
a3tb2 = a2tb2.reshape(nytmp,nxtmp,-1)
a3tbobs = concatenate([a3tb1, a3tb2],axis=2)
a3tbobs = a3tbobs[iy:ey,:,:]
a1tbobs = a3tbobs[y,x,:]
print 'tb obs',a1tbobs
#****************************************************
# Convert Tb to EPC
#----------------------------------------------------
print 'calc epc'
a1epcobs = epcfunc.mk_epc_12pc(a1tbobs.reshape(1,1,13), a2coef).reshape(-1,)
print 'calc epc done'
#print a3epc
#****************************************************
# Find EPC bin numbers
#----------------------------------------------------
#NPCHIST = 29
#print 'calc idx'
##a2idx_db = epcfunc.mk_epc_id_25bins(a3epc, a2pc_edge)
#a2idx_db = epcfunc.mk_epc_id_nbins(a3epc, a2pc_edge, NPCHIST)
#print 'calc idx done'


#********************************************
#-- Read Obs environmental variables  --
#********************************************
matchbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A'
t2mPath = matchbaseDir + '/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid)
tqvPath = matchbaseDir + '/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid)


t2mobs = np.load(t2mPath)[iy+y,x]
tqvobs = np.load(tqvPath)[iy+y,x]
print 'obs epc',a1epcobs
print 'obs t2mv=%.1f'%t2mobs, ' tqv=%.4f'%tqvobs

#********************************************
idx_db = int(a2idxdb[y,x])
irec   = int(a2irec[y,x])
print idx_db

dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
db.set_idx_db(dbDir, idx_db)

a2tbdb = db.get_var('tb')
a2epcdb= db.get_var('pc_emis')
a1nsurfdb = db.get_var('precip_nsfc_NS')
a1t2mdb = db.get_var('t2m')
a1tqvdb = db.get_var('tqv')
#print a2tbdb[irec]
#print a1nsurfdb[irec]


#--- sort database based on nsurf ---
a1k = np.arange(a2tbdb.shape[0]).astype(int32)
a2dat = np.concatenate([a1nsurfdb.reshape(-1,1), a1k.reshape(-1,1)],axis=1)
a2dat = np.array(sorted(a2dat, key=lambda x: (x[0]), reverse=True))

a1nsurfdbsort = a2dat[:,0]
a1ksort = a2dat[:,1].astype(int32)
print a1ksort
for i,k in enumerate(a1ksort[:10]):
    a1epcdb = a2epcdb[k]
    rmsd = np.sqrt(np.square((a1epcdb - a1epcobs)/a1pc_std).sum()/NEM)
    print ''
    print i, 'rmsd=','%.4f'%rmsd,'nsurfNS=',a1nsurfdb[k]
    print 't2m_db=%.4f'%a1t2mdb[k], ' tqv_db=%.4f'%(a1tqvdb[k])
    print 'EPC',a2epcdb[k]

    


