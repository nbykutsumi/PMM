from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l1_gmi as l1_gmi
import mkdbfunc
import glob
from datetime import datetime, timedelta
import numpy as np
import sys

NPC_USE = 3
#NPC_USE = 4
if NPC_USE==3:
    #NBINS=25
    NBINS=22
elif NPC_USE==4:
    NBINS=10
else:
    print 'check NBINS',NPC_USE
    sys.exit()

percentiletype = 'all'
#percentiletype = 'precip'

gmi    = l1_gmi.L1_GMI()

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,12,31)
#eDTime = datetime(2017,1,1)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

cx  = 110  # GMI center angle bin (py-idx) # GMI total angle bin=221
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

#------------------------------------------
def mk_tpc_id(a3tpc, NBINS):
    ny,nx,nz  = a3tpc.shape
    a2id_db  = zeros([ny,nx],int32)
    for iem in range(NPC_USE):
        a1bin = a2pc_range[iem]
        a2idTmp = np.digitize(a3tpc[:,:,iem], a1bin, right=False) - 1
    
        a2idTmp = ma.masked_outside(a2idTmp,0,NBINS-1)
        a2id_db = a2id_db + a2idTmp*pow(NBINS, NPC_USE-1-iem)

    a2id_db = a2id_db.filled(-9999)
    return a2id_db


def mk_tpc_id_25bins(a3tpc):
    ny,nx,nz  = a3tpc.shape
    a2id_db  = zeros([ny,nx],int32)
    for iem in range(NPC_USE):
        a1bin = a2pc_range[iem]
        a2idTmp = np.digitize(a3tpc[:,:,iem], a1bin, right=False) - 1
    
        a2idTmp = ma.masked_outside(a2idTmp,0,24)
        a2id_db = a2id_db + a2idTmp*pow(25, NPC_USE-1-iem)

    a2id_db = a2id_db.filled(-9999)
    return a2id_db

def mk_tpc_id_10bins(a3tpc):
    ny,nx,nz  = a3tpc.shape
    a2id_db  = zeros([ny,nx],int32)
    for iem in range(NPC_USE):
        a1bin = a2pc_range[iem]
        a2idTmp = np.digitize(a3tpc[:,:,iem], a1bin, right=False) - 1
    
        a2idTmp = ma.masked_outside(a2idTmp,0,9)
        a2id_db = a2id_db + a2idTmp*pow(10, NPC_USE-1-iem)

    a2id_db = a2id_db.filled(-9999)
    return a2id_db


#-- Read EPC range files --
coefDir   = '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
rangePath = coefDir + '/PC_MIN_MAX_%d_no_overlap_%scases.txt'%(NBINS, percentiletype)
a2tmp     = read_table(rangePath)
a2tmp     = a2tmp[:,3:]

nytmp, nxtmp = a2tmp.shape
a2pc_range     = zeros([nytmp,nxtmp/2+1]) # 12*(25+1)
a2pc_range_min = zeros([nytmp,nxtmp/2])   # 12*25
a2pc_range_max = zeros([nytmp,nxtmp/2])   # 12*25
for ytmp in range(nytmp):
    for xtmp in range(0,nxtmp,2):
        a2pc_range_min[ytmp,xtmp/2] = a2tmp[ytmp,xtmp]
        a2pc_range_max[ytmp,xtmp/2] = a2tmp[ytmp,xtmp+1]

a2pc_range[:,:-1] = a2pc_range_min
a2pc_range[:,-1]  = a2pc_range_max[:,-1]
a2pc_range[:,0]   = a2pc_range[:,0] - 1.e6
a2pc_range[:,-1]  = a2pc_range[:,-1] + 1.e6

#-- Read EPC conversion Coefficient file --
coefDir   = '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
#coefPath = coefDir + '/coef_pc.txt'
#a2coef   = read_table(coefPath)
#a2coef   = a2coef[:,1:]
coefPath  = coefDir + '/egvec.nltb.npy'
a2coef    = np.load(coefPath)[:NPC_USE,:]  # (n-th, nComb)

#-- Read Mean and STD file --
coefDir   = '/work/hk01/utsumi/PMM/TPCDB/PC_COEF'
meanPath  = coefDir + '/mean.nltb.npy'
stdPath   = coefDir + '/std.nltb.npy'
a1mean = np.load(meanPath)
a1std  = np.load(stdPath)
#-----------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    '''
    # Do not directly read GMI HDF file.
    # It may cause mismuch between Tc(S1) and Tc(S2)
    '''

    tcbaseDir  = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, 'Tc')
    tc2baseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, 'TcS2')

    tcDir      = tcbaseDir  + '/%04d/%02d/%02d'%(Year,Mon,Day)
    tc2Dir     = tc2baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)

    ssearchTc  = tcDir + '/Tc.*.npy'
    lsrcPathTc = sort(glob.glob(ssearchTc))
    
    if len(lsrcPathTc)==0:
        print 'No Tc file',Year,Mon,Day
        print ssearchTc
        sys.exit()
    
    for srcPathTc in lsrcPathTc:
        oid  = srcPathTc.split('.')[-2]
        srcPathTcS2 = tc2Dir + '/TcS2.1.%s.npy'%(oid)

        a3tb_s1 = np.load(srcPathTc)
        a3tb_s2 = np.load(srcPathTcS2)

        a3tb = concatenate([a3tb_s1, a3tb_s2], axis=2)


        #--- Make non-linear combination of Tb --
        nlen = a3tb.shape[0]
        a3comb= mkdbfunc.mk_nonlin_comb(a3tb)

        #--- Convert to PC of nl-Tb -------------
        a3comb= (a3comb - a1mean)/a1std
        a3tpc = np.dot(a3comb, a2coef.T)

        #--- Make PC index ---------------------
        #a2tpcid = mk_tpc_id_25bins(a3tpc)
        a2tpcid = mk_tpc_id(a3tpc, NBINS)

        #--- Mask where Tb is invalid ----------
        a2masktb= ma.masked_outside(a3tb,50, 350).mask.any(axis=2)
        a2tpcid = ma.masked_where(a2masktb, a2tpcid).filled(-9999)

        #-- save tpcid --
        outvarName = 'tpcid-s1-%dpc-%dbin'%(NPC_USE,NBINS)
        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, mwscan, cx-w, cx+w, outvarName)

        outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        outPath    = outDir + '/tpcid-s1.%s.npy'%(oid)

        util.mk_dir(outDir)
        np.save(outPath, a2tpcid.astype(int16))
        print outPath

    
