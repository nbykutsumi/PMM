from numpy import *
import os, sys
import myfunc.util as util
import EPCDB
import socket
import epcfunc

DB_MAXREC = 10000
sensor = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    gmibaseDir  = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/tank/utsumi/PMM/MATCH.GMI.V05A'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    #outbaseDir = '/tank/utsumi/PMM/retepc/%s'%(expr)
elif myhost =="well":
    #gmibaseDir  = '/media/disk2/share/data/PMM/NASA/GPM.GMI/1C/V05'
    gmibaseDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    #outbaseDir = '/media/disk2/share/PMM/retepc/%s'%(expr)
    #outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s'%(expr)

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

#--------------------------------------------------------------
NPCHIST = 29

db = EPCDB.EPCDB()

idx_db = 6628
db.set_idx_db(dbDir, idx_db)
a2tb = db.get_var('tb')
ny,nz = a2tb.shape


#****************************************************
# Convert Tb to EPC
#----------------------------------------------------
print 'calc epc'
a3epc = epcfunc.mk_epc_12pc(a2tb.reshape(ny,1,nz), a2coef)
print 'calc epc done'
print a3epc



#****************************************************
# Find EPC bin numbers
#----------------------------------------------------
print 'calc idx'
#a2idx_db = epcfunc.mk_epc_id_25bins(a3epc, a2pc_edge)
a2idx_db = epcfunc.mk_epc_id_nbins(a3epc, a2pc_edge, NPCHIST)
print 'calc idx done'
print a2idx_db
print a2idx_db.min(), a2idx_db.max()


