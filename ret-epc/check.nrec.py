from numpy import *
import os, sys
import socket
import myfunc.util as util

myhost = socket.gethostname()
if myhost =="shui":
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    outbaseDir = '/home/utsumi/temp/ret'
elif myhost =="well":
    #gmibaseDir  = '/media/disk2/share/data/PMM/NASA/GPM.GMI/1C/V05'
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    #outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/hometemp/ret'
    outbaseDir = '/home/utsumi/temp/ret'

def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return array(lout)

def read_nrain(idx_db):
    '''
     /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
     /* second six= same for when T2m > 278K */
    '''
    #srcPath = dbDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    nrainDir= dbDir + '/nrain'
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    a1nrain = read_table(srcPath,type=int32)[0]
    a1nrain_cold = a1nrain[:6]
    a1nrain_warm = a1nrain[6:]
    return a1nrain_warm, a1nrain_cold


ndb = 29*29*29  # 24388
sout = ''
for idx_db in range(ndb):
    nrainDir= dbDir + '/nrain'
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    if not os.path.exists(srcPath):
        continue

    print idx_db
    _,a1nrain = read_nrain(idx_db)
    stmp = '%d,'%(idx_db) + ','.join(map(str,a1nrain)) + '\n'
    sout = sout + stmp

outPath = outbaseDir + '/nrain.csv'
util.mk_dir(outbaseDir)
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
