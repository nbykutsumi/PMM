from numpy import *
import calendar
import myfunc.util as util
import numpy as np
from collections import deque
import sys, os, glob
import gc
import h5py

prod     = '1C'
verGMI   = '05'
subverGMI= 'A'
fullverGMI='%s%s'%(verGMI, subverGMI)

iYM = [2017,12]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15624]]  # 25*25*25 = 15625
lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15000],[15000,17500],[17500,20000],[20000,22500],[22500,25000]]  # 29*29*29 = 24389


cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)
worg= 221  # GMI total angle bins

lvar = ['t2m']

'''
int8 : -128 ~ +127
int16: -32768 ~ +32767
int32: -2147483648 ~ +2147483647
'''
dattype={
 't2m' :'float32'
}

dnvect ={
 't2m' :1
}

#-- if vector size=11 --
maxmem = 3*1000*1000*1000 # N GB
dmaxrec= {}
for var in lvar:
    nvect = dnvect[var]
    dmaxrec[var] = int(maxmem / (4*nvect*10000))

#-----------------------

listDir  = '/work/hk01/utsumi/PMM/EPCDB/list'
extractidDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, 'S1', cx-w, cx+w, 'epcid-s1')

gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'


outbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)

#-- Functions --
def csv2list(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split(',')
        lout.append(line)
    return lout
#---------------
for Year,Mon in lYM:
    for  varName in lvar:
        for epcid_range in lepcid_range:
            epcid_min, epcid_max = epcid_range
            maxrec   = dmaxrec[var]
            grp      = '/'.join(var.split('/')[:-1])
            listPath = listDir + '/list.shuffle.%s.V%s.%04d%02d.csv'%(prod,verGMI,Year,Mon)
            lobt = csv2list(listPath)
        
            dstack = {}
            dnum   = {}
            for epcid in range(epcid_min,epcid_max): 
                dstack[epcid] = deque([])
                dnum  [epcid] = 0
    
            ##-- test --
            #lobt = lobt[:4]
            ##----------    
            for (obtnum, Year,Mon,Day,time0, time1) in lobt:   
                Year,Mon,Day = map(int, [Year,Mon,Day])
                print varName, epcid_range,obtnum,Year,Mon,Day

                #-- Read EPC-id --
                extractDir = extractidDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                ssearch = extractDir + '/epcid-s1.%s.npy'%(obtnum)
                epcidPath = glob.glob(ssearch)[0]
                a2epcid   = np.load(epcidPath)
                a1epcid   = a2epcid.flatten()
    
   
                #-- Read Output variables --
                if varName  in ['t2m']:
                    baseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.%s'%(varName)
                    srcDir  = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    srcPath = srcDir + '/%s.%s.npy'%(varName, obtnum)
                    Dat     = np.load(srcPath)[:,cx-w:cx+w+1].flatten()
 
                else:
                    print 'check varName',varName 
                    sys.exit()

                #-- stack to epc-classifiled database ---
                lepcidset = sort(list(set(a1epcid)))
                lepcidset = ma.masked_less(lepcidset, epcid_min)
                lepcidset = ma.masked_greater_equal(lepcidset, epcid_max)
                lepcidset = lepcidset.compressed()

                #print 'epc_min,epc_max,setlen=',epcid_min,epcid_max,len(lepcidset)
                for epcid in lepcidset:
                    if epcid ==-9999: continue
                    a1bool = ma.masked_equal(a1epcid, epcid).mask
                    DatTmp  = Dat[a1bool]
                    dstack[epcid].extend(DatTmp.tolist())
                    #-- save large arrays ---
                    if (len(dstack[epcid]))>maxrec:
                        dnum[epcid] = dnum[epcid] +1 
        
                        outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                        outPath = outDir + '/%s.%05d.%d.npy'%(varName,epcid, dnum[epcid])
                        util.mk_dir(outDir)
                        np.save(outPath, array(dstack[epcid]).astype(dattype[var]))
        
                        del dstack[epcid]
                        dstack[epcid] = deque([])
                        gc.collect()
                         
                        #print outPath
        
            #--- save ---
            for epcid in range(epcid_min,epcid_max):
                dnum[epcid] = dnum[epcid] + 1 
                if len(dstack[epcid])>0:
                    outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                    outPath = outDir + '/%s.%05d.%d.npy'%(varName,epcid,dnum[epcid])
                    util.mk_dir(outDir)
    
                    np.save(outPath, array(dstack[epcid]).astype(dattype[var]))
    
                    del dstack[epcid]
                    gc.collect()
    
    
            #--- Joint segments ---
            for epcid in range(epcid_min,epcid_max):
                outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                ssearch = outDir + '/%s.%05d.*.npy'%(varName,epcid)
                lsrcPath= sort(glob.glob(ssearch))
    
                if len(lsrcPath)==0: continue
    
                aout = deque([])
                for srcPath in lsrcPath:
                    atmp = np.load(srcPath)
                    aout.extend(atmp)
    
                outPath =  outDir + '/%s.%05d.npy'%(varName,epcid)
    
                np.save(outPath, aout)
                print outPath
    
                #-- delete temporary files --
                for srcPath in lsrcPath:
                    os.remove(srcPath)
    


