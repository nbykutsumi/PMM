from numpy import *
import calendar
import myfunc.util as util
import numpy as np
from collections import deque
import sys, os, glob
import gc
import h5py

prod      = '1C'
verGMI    = '05'
subverGMI = 'A'
fullverGMI= '%s%s'%(verGMI, subverGMI)
verDPR    = '06'
subverDPR = 'A'
fullverDPR= '%s%s'%(verDPR, subverDPR)

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15624]]  # 25*25*25 = 15625



cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
#cw  = 27    # extract this width around center
w   = int(cw/2)
worg= 221  # GMI total angle bins


#lvar = [['Ku','NS/SLV/zFactorCorrected'],['Ku','NS/SLV/precipRate']]
#lvar = [['Ku','NS/SLV/zFactorCorrected']]
#lvar = [['Ku','NS/SLV/precipRate']]
lvar = [['Ku','NS/SLV/precipRateESurface']]
#lvar = [['Ku','NS/PRE/elevation']]



'''
int8 : -128 ~ +127
int16: -32768 ~ +32767
int32: -2147483648 ~ +2147483647
'''
dattype={
 'NS/SLV/zFactorCorrected': 'float32'
,'NS/SLV/precipRate':       'float32'
,'NS/SLV/precipRateESurface':'float32'
,'NS/PRE/elevation':        'float32'
}

dnvect ={
 'NS/SLV/zFactorCorrected': 50
,'NS/SLV/precipRate':       50
,'NS/SLV/precipRateESurface':1
,'NS/PRE/elevation':        1

}

#---------------------
maxmem = 3*1000*1000*1000 # N GB
dmaxrec= {}
for (radar, var) in lvar:
    nvect = dnvect[var]
    dmaxrec[var] = int(maxmem / (4*nvect*10000))

#-----------------------

listDir  = '/work/hk01/utsumi/PMM/EPCDB/list'
extractidDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, 'S1', cx-w, cx+w, 'epcid-s1')
extractKuIdxDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/S1.ABp083-137.Ku.V06A.IDX'%(fullverGMI)

#gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
#gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'

outbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)

#-- Functions --
def csv2list(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split(',')
        lout.append(line)
    return lout

def average_2ranges_3d(a3in,miss=None,dtype=float32, fill=True):
    '''
    a3in: (ny, nx, nz) --> a2out: (ny, nx, nz/2)
    nz should be an even number
    '''
    ny,nx,nz = a3in.shape
    a4out = empty([2,ny,nx,nz/2], dtype)
    a4out[0] = a3in[:,:,::2]
    a4out[1] = a3in[:,:,1::2]
    if fill==True:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).filled(miss).astype(dtype)
    else:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).astype(dtype)
    return a3out

def average_4ranges_3d(a3in,miss=None,dtype=float32, fill=True):
    '''
    a3in: (ny, nx, nz) --> a2out: (ny, nx, nz/4)
    nz should be an even number
    '''
    ny,nx,nz = a3in.shape
    a4out = empty([4,ny,nx,nz/4], dtype)
    a4out[0] = a3in[:,:,::4]
    a4out[1] = a3in[:,:,1::4]
    a4out[2] = a3in[:,:,2::4]
    a4out[3] = a3in[:,:,3::4]
    if fill==True:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).filled(miss).astype(dtype)
    else:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).astype(dtype)
    return a3out

def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr,nzdpr= a3in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask
        
        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,miss).mean(axis=0)
    a2datTmp[a1dprmask,:] = miss
    return a2datTmp


def ave_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr = a2in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a2datTmp    = empty([9,len(a1y)], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask
        
        a1datTmp= a2in[a1yTmp.filled(0),a1xTmp.filled(0)]

        a2datTmp[itmp,:] = a1datTmp


    a1datTmp = ma.masked_equal(a2datTmp,miss).mean(axis=0)
    a1datTmp[a1dprmask] = miss
    return a1datTmp



#---------------
for Year,Mon in lYM:
    for (radar, var) in lvar:
        for epcid_range in lepcid_range:
            epcid_min, epcid_max = epcid_range
            maxrec   = dmaxrec[var]
            grp      = '/'.join(var.split('/')[:-1])
            varName  = var.split('/')[-1]
            listPath = listDir + '/list.%s.V%s.%04d%02d.csv'%(prod,verGMI,Year,Mon)
            lobt = csv2list(listPath)
        
            dstack = {}
            dnum   = {}
            for epcid in range(epcid_min,epcid_max): 
                dstack[epcid] = deque([])
                dnum  [epcid] = 0
    
            #-- test --
            #lobt = lobt[:4]
            #lobt = lobt[3:3+1]
            #----------    
            for (obtnum, Year,Mon,Day,time0, time1) in lobt:   
                Year,Mon,Day = map(int, [Year,Mon,Day])
                print obtnum,Year,Mon,Day
                #-- Read EPC-id --
                extractDir = extractidDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                epcidPath = extractDir + '/epcid-s1.%s.npy'%(obtnum)
                a2epcid   = np.load(epcidPath)
                a1epcid   = a2epcid.flatten()

                #-- Read DPR.IDX --
                extractDir = extractKuIdxDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                xPath     = extractDir + '/Xpy.1.%s.npy'%(obtnum)
                yPath     = extractDir + '/Ypy.1.%s.npy'%(obtnum)
                a2dprx    = np.load(xPath)[:,cx-w-83:cx+w+1-83]
                a2dpry    = np.load(yPath)[:,cx-w-83:cx+w+1-83]

                a1dprx    = a2dprx.flatten()
                a1dpry    = a2dpry.flatten()

                #-- Search DPR granules --
                dprbaseDir   = '/work/hk01/PMM/NASA/GPM.%s/2A/V%s'%(radar, verDPR)
                if  (radar=='Ku'):
                    srcDir = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch= srcDir + '/*.%s.V%s.HDF5'%(obtnum,fullverDPR)
                    srcPath= glob.glob(ssearch)[0]
                else:
                    print 'check radar',radar
                    sys.exit()

                #-- Read DPR variables --
                if   varName =='zFactorCorrected':
                    h    = h5py.File(srcPath)
                    Dat0 = h[var][:]
                    h.close() 

                    #-- Convert dBZ --> Z ---
                    '''
                    #XdBZ = 10*log_10(Z/Z0)
                    #Z    = Z0* 10^(X*0.1)
                    #No need to multiply Z0.
                    #(Because they are devided by Z0 when converted to dBZ at the end)
                    '''
                    DatLN = np.power(10, 0.1*ma.masked_equal(Dat0,-9999.9)).filled(-9999.9)

                    #-- Average vertical ranges (over Linearlized Z)--
                    '''
                    # original: vertical 176bins (with 125m resol)
                    # 0-10km  (-80:last bin): 125m ->250m (total 40bin)
                    # 10-15km (-120:-80 bin): 125m ->500m (total 10bin)
                    '''
                    Dat_up  = average_4ranges_3d(DatLN[:,:,-120:-80],miss=-9999.9,dtype=float32, fill=False)
                    Dat_low = average_2ranges_3d(DatLN[:,:,-80:],miss=-9999.9, dtype=float32, fill=False)
                    DatLN   = concatenate([Dat_up, Dat_low],axis=2)

                    #-- Average 9 grids (over Linearlized Z)--
                    a2datTmp = ave_9grids_3d(DatLN, a1dpry, a1dprx, miss=-9999.9)

                    #-- Convert to dBZ --
                    Dat = 10*np.log10(ma.masked_equal(a2datTmp, -9999.9))
                    Dat = ma.masked_invalid(Dat).filled(-9999.9).astype(float32)


                elif   varName =='precipRate':
                    h    = h5py.File(srcPath)
                    Dat0 = h[var][:]
                    h.close() 

                    #-- Average vertical ranges (over Linearlized Z)--
                    Dat_up  = average_4ranges_3d(Dat0[:,:,-120:-80],miss=-9999.9,dtype=float32, fill=False)
                    Dat_low = average_2ranges_3d(Dat0[:,:,-80:],miss=-9999.9, dtype=float32, fill=False)
                    DatCnc  = concatenate([Dat_up, Dat_low],axis=2)

                    #-- Average 9 grids (over Linearlized Z)--

                    Dat     = ave_9grids_3d(DatCnc, a1dpry, a1dprx, miss=-9999.9)

                elif varName in ['precipRateESurface']:
                    h   = h5py.File(srcPath)
                    Dat0 = h[var][:]
                    h.close()
                    if   dattype[var]=='int32':
                        miss=-9999
                    elif dattype[var]=='float32':
                        miss=-9999.9
   
                    Dat     = ave_9grids_2d(Dat0, a1dpry, a1dprx, miss=miss)
                else:
                    print 'check varName',varName
                    sys.exit()
                
                #-- stack to epc-classifiled database ---
                lepcidset = sort(list(set(a1epcid)))
                lepcidset = ma.masked_less(lepcidset, epcid_min)
                lepcidset = ma.masked_greater_equal(lepcidset, epcid_max)
                lepcidset = lepcidset.compressed()

                print 'epc_min,epc_max,setlen=',epcid_min,epcid_max,len(lepcidset)
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
                         
                        print outPath
        
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
    


