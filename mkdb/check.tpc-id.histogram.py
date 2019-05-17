from numpy import *
import numpy as np
import myfunc.util as util
import random

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
NPC_USE = 3
NBINS   = 22

#pctype = 'tpc-3pc'
pctype = 'tpcid-s1-%dpc-%dbin'%(NPC_USE,NBINS)
#pctype = 'epc'
#-- Read orbit list  ---
if pctype.split('-')[0] =='tpcid':
    listDir = '/work/hk01/utsumi/PMM/TPCDB/list'
elif pctype.split('-')[0]=='epcid':
    listDir = '/work/hk01/utsumi/PMM/EPCDB/list'
else:
    print 'check pctype'
    sys.exit()

lorbit  = []
for (Year,Mon) in lYM:
    listPath = listDir + '/list.1C.V05.%04d%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines = f.readlines(); f.close()
    for line in lines:
        line = map(int, line.split(','))
        #print line
        lorbit.append(line)

#-- Sampling fron orbit list ---
random.seed(0)
#nsample = 4
nsample = int(len(lorbit))
lorbitTmp = random.sample(lorbit, nsample)

#if pctype in ['tpc-3pc','epc']:
#    lbnd = arange(0,25*25*25+1)-0.5
#elif pctype == 'tpc-4pc':
#    lbnd = arange(0,10**4+1)-0.5

lbnd = arange(0,NBINS**NPC_USE+1)-0.5

a1hist = zeros(len(lbnd)-1,int32)
for orbinfo in lorbit:
    oid,Year,Mon,Day,itime,etime = orbinfo
    #srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.tpcid-s1/2017/06/05'


    srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.%s/%04d/%02d/%02d'%(pctype,Year,Mon,Day)
    srcPath= srcDir + '/%s.%06d.npy'%(pctype,oid)
    #if pctype=='tpc-3pc':
    #    srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.tpcid-s1/%04d/%02d/%02d'%(Year,Mon,Day)
    #    srcPath= srcDir + '/tpcid-s1.%06d.npy'%(oid)
    #elif pctype=='tpc-4pc':
    #    srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.tpcid-s1-4pc/%04d/%02d/%02d'%(Year,Mon,Day)
    #    srcPath= srcDir + '/tpcid-s1-4pc.%06d.npy'%(oid)

    #elif pctype=='epc':
    #    srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.epcid-s1/%04d/%02d/%02d'%(Year,Mon,Day)
    #    srcPath= srcDir + '/epcid-s1.%06d.npy'%(oid)



    a2id = np.load(srcPath)
     
    a1histTmp,atmp = np.histogram(a2id, bins=lbnd)
    a1hist = a1hist + a1histTmp

outPath = '/home/utsumi/temp/temp.%s.hist.csv'%(pctype)
sout = util.array2csv(a1hist.reshape(-1,1))
f=open(outPath,'w'); f.write(sout); f.close()
print outPath


