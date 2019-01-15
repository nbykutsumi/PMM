import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l1_gmi as l1_gmi

import glob
from datetime import datetime, timedelta
import numpy as np
import sys
from f_match_fov import *

NEM     = 11
NTBREG  = 9
NEM_USE = 4

gmi    = l1_gmi.L1_GMI()

iDTime = datetime(2017,1,2)
eDTime = datetime(2017,1,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)

verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
mwscan = 'S1'

outrootDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s'%(fullverGMI)

var = 'S1/Tc'

#-- Functions ------
def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return array(lout)

def mk_epc_9ch(a3tb):
    coefDir  = '/home/utsumi/bin/ENSPR'
    coefPath = coefDir + '/coef_pc.txt'
    a2coef   = read_table(coefPath)
    a2coef   = a2coef[:,1:]
    
    ny,nx,ntmp = a3tb.shape
    a3epc = zeros([ny,nx,NEM])
    
    # +a0
    kt  = 0
    for iem in range(NEM):
        a3epc[:,:,iem] = a2coef[kt,iem]
    
    for i in range(NTBREG):   # NTBRERG=9
        print 'iTBREG=',i,'out of',NTBREG
        kt =kt+1
    
        # + bi*Tbi
        for iem in range(NEM): # NEM = 11
            b = a2coef[kt, iem]
            a3epc[:,:,iem] = a3epc[:,:,iem] + b*a3tb[:,:,i]
    
            ##-- test --
            #if iem==emtemp:
            #    print 'i,b, Tb',i,b,a3tb[iy,ix,i]
            ##----------
    
        # + cij * Tbi * Tbj
        for j in range(i,NTBREG):
            kt = kt+1
            #print i,j
            a2tmp = a3tb[:,:,i] * a3tb[:,:,j]
    
            for iem in range(NEM):
                c = a2coef[kt, iem]
                a3epc[:,:,iem] = a3epc[:,:,iem] + c*a3tb[:,:,i]*a3tb[:,:,j]
    
                ##-- test --
                #if iem==emtemp:
                #    print 'i,j,c, Tbi*Tbj',i,j,c,a3tb[iy,ix,i]*a3tb[iy,ix,j]
                ##----------
    
    
    # + d * (V-H)/(V+H)
    kt = kt+1
    for iem in range(NEM):
        d = a2coef[kt,iem]
        a3epc[:,:,iem] = a3epc[:,:,iem] \
             +d*(a3tb[:,:,0]-a3tb[:,:,1])/(a3tb[:,:,0]+a3tb[:,:,1])
    
    kt = kt+1
    for iem in range(NEM):
        d = a2coef[kt,iem]
        a3epc[:,:,iem] = a3epc[:,:,iem] \
             +d*(a3tb[:,:,2]-a3tb[:,:,3])/(a3tb[:,:,2]+a3tb[:,:,3])
    
    kt = kt+1
    for iem in range(NEM):
        d = a2coef[kt,iem]
        a3epc[:,:,iem] = a3epc[:,:,iem] \
             +d*(a3tb[:,:,5]-a3tb[:,:,6])/(a3tb[:,:,5]+a3tb[:,:,6])
    
    return a3epc

def mk_epc_id(a3epc):
    ny,nx,nz  = a3epc.shape
    a2id_db  = zeros([ny,nx],int32)
    for iem in range(NEM_USE):
        a1bin = a2pc_range[iem]
        a2idTmp = np.digitize(a3epc[:,:,iem], a1bin, right=False) - 1
    
        a2idTmp = ma.masked_outside(a2idTmp,0,9)
        a2id_db = a2id_db + a2idTmp*pow(10, iem)

    a2id_db = a2id_db.filled(-9999)
    return a2id_db



#-- Read EPC range files --
coefDir   = '/home/utsumi/bin/ENSPR'
rangePath = coefDir + '/PC_MIN_MAX_10_no_overlap.txt'
a2tmp     = read_table(rangePath)
a2tmp     = a2tmp[:,3:]

nytmp, nxtmp = a2tmp.shape
a2pc_range     = zeros([nytmp,nxtmp/2+1]) # 11*(10+1)
a2pc_range_min = zeros([nytmp,nxtmp/2])   # 11*10
a2pc_range_max = zeros([nytmp,nxtmp/2])   # 11*10
for ytmp in range(nytmp):
    for xtmp in range(0,nxtmp,2):
        a2pc_range_min[ytmp,xtmp/2] = a2tmp[ytmp,xtmp]
        a2pc_range_max[ytmp,xtmp/2] = a2tmp[ytmp,xtmp+1]

a2pc_range[:,:-1] = a2pc_range_min
a2pc_range[:,-1]  = a2pc_range_max[:,-1]


#-----------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
    srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.HDF5'

    lsrcPathGMI = glob.glob(ssearchGMI)
    
    if len(lsrcPathGMI)==0:
        print 'No GMI file',Year,Mon,Day
        print ssearchGMI
        sys.exit()
    
    for srcPathGMI in lsrcPathGMI:
        oid  = srcPathGMI.split('.')[-3]
        datoutTmp = gmi.load_var_granule(srcPathGMI,var)
        a3tb = datoutTmp[:,cx-w:cx+w+1,:]
        print 'make epc'
        a3epc= mk_epc_9ch(a3tb)
        print 'make epc done'
        a2epcid = mk_epc_id(a3epc)

        #-- save epc --
        outvarName    = 'epc-s1'
        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, outvarName)

        outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        outPath    = outDir + '/%s.%s.npy'%(outvarName, oid)

        util.mk_dir(outDir)
        np.save(outPath, a3epc.astype(float32))
        print outPath

        #-- save epcid --
        outvarName    = 'epcid-s1'
        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, outvarName)

        outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        outPath    = outDir + '/%s.%s.npy'%(outvarName, oid)

        util.mk_dir(outDir)
        np.save(outPath, a2epcid.astype(int16))
        print outPath

    
