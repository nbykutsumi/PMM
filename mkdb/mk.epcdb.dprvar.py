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

iYM = [2017,12]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15624]]  # 25*25*25 = 15625
#lepcid_range = [[0,5000],[5000,10000],[10000,15000],[15000,20000],[20000,25000]]  # 29*29*29 = 24389
lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15000],[15000,17500],[17500,20000],[20000,22500],[22500,25000]]  # 29*29*29 = 24389

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
#cw  = 27    # extract this width around center
w   = int(cw/2)
worg= 221  # GMI total angle bins


#lvar = [['Ku','NS/SLV/zFactorCorrected'],['Ku','NS/SLV/precipRate']]
#lvar = [['Ku','NS/SLV/zFactorCorrected']]
#lvar = [['Ku','NS/PRE/zFactorMeasured'],['Ka','MS/PRE/zFactorMeasured'],['Ku','NS/PRE/elevation']]
#lvar = [['Ka','MS/PRE/zFactorMeasured']]
#lvar = [['Ku','NS/PRE/zFactorMeasured']]
#lvar = [['DPRGMI','NS/precipTotWaterCont']]
#lvar = [['DPRGMI','NS/surfPrecipTotRate']]
lvar = [['DPRGMI','MS/surfPrecipTotRate']]
#lvar = [['DPRGMI','NS/surfPrecipTotRate'],['DPRGMI','NS/precipTotWaterCont']]
#lvar = [['Ku','NS/SLV/precipRate']]
#lvar = [['Ku','NS/SLV/precipRateESurface']]
#lvar = [['Ku','NS/SLV/precipRateNearSurface']]
#lvar = [['Ka','MS/SLV/precipRateNearSurface']]
#lvar = [['Ku','NS/PRE/elevation'],['Ku','NS/CSF/typePrecip'],['Ku','NS/PRE/heightStormTop']]
#lvar = [['Ku','NS/PRE/heightStormTop']]



'''
int8 : -128 ~ +127
int16: -32768 ~ +32767
int32: -2147483648 ~ +2147483647
'''
dattype={
 'NS/SLV/zFactorCorrected': 'int16'  # scaled by 100, (125m->250m)
,'MS/SLV/zFactorCorrected': 'int16'  # scaled by 100, (125m->250m)
,'NS/PRE/zFactorMeasured' : 'int16'  # scaled by 100, (125m->250m)
,'MS/PRE/zFactorMeasured' : 'int16'  # scaled by 100, (125m->250m)
,'NS/SLV/precipRate':       'int16'  # scaled by 100, (125m->250m)
,'NS/SLV/precipRateESurface':'float32'
,'NS/PRE/elevation':        'int16'  # float32 --> int16
,'NS/CSF/typePrecip':       'int16'  # convert to int16 (-32768, 32767)
,'NS/PRE/heightStormTop':   'int16'  # save as int16 (-32768, 32767)

,'NS/SLV/precipRateNearSurface': 'float32'
,'MS/SLV/precipRateNearSurface': 'float32'

,'NS/precipTotWaterCont':   'float32'  # original:250m
,'MS/precipTotWaterCont':   'float32'  # original:250m
,'MS/surfPrecipTotRate' :   'float32' # DPRGMI
,'NS/surfPrecipTotRate' :   'float32' # DPRGMI

}

dnvect ={
 'NS/SLV/zFactorCorrected': 60
,'MS/SLV/zFactorCorrected': 60
,'NS/PRE/zFactorMeasured': 60
,'MS/PRE/zFactorMeasured': 60
,'NS/SLV/precipRate':       60
,'NS/SLV/precipRateESurface':1
,'NS/SLV/precipRateNearSurface':1
,'MS/SLV/precipRateNearSurface':1
,'NS/PRE/elevation':        1
,'NS/CSF/typePrecip':       1
,'NS/PRE/heightStormTop':   1

,'NS/precipTotWaterCont':  60
,'MS/precipTotWaterCont':  60
,'MS/surfPrecipTotRate' :  1
,'NS/surfPrecipTotRate' :  1
}

#---------------------
maxmem = 2*1000*1000*1000 # N GB
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

def sum_9grids_2d(a2in, a1y, a1x, miss):
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

    a1datTmp = ma.masked_equal(a2datTmp,miss).sum(axis=0)
    a1datTmp[a1dprmask] = miss
    return a1datTmp



#---------------
for Year,Mon in lYM:
    for (radar, var) in lvar: 
        for epcid_range in lepcid_range:

            #lepcskip = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15000],[15000,17500],[17500,20000],[20000,22500],[22500,25000]]  # 29*29*29 = 24389
            #lepcskip = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15000],[15000,17500]]  # 29*29*29 = 24389
            #if ((Mon ==2) & (epcid_range in lepcskip)): continue # test

            epcid_min, epcid_max = epcid_range
            maxrec   = dmaxrec[var]
            grp      = '/'.join(var.split('/')[:-1])
            scanName = var.split('/')[0]
            varName  = var.split('/')[-1]
            varNameOut='%s_%s_%s'%(radar,scanName,varName)
            nvect    = dnvect[var]
            listPath = listDir + '/list.shuffle.%s.V%s.%04d%02d.csv'%(prod,verGMI,Year,Mon)
            lobt = csv2list(listPath)
        
            dstack = {}
            dnum   = {}
            for epcid in range(epcid_min,epcid_max): 
                dstack[epcid] = deque([])
                dnum  [epcid] = 0
    
            #-- test --
            #lobt = lobt[:3]
            #lobt = lobt[3:3+1]
            #----------    
            for (obtnum, Year,Mon,Day,time0, time1) in lobt:   
                Year,Mon,Day = map(int, [Year,Mon,Day])
                print var, obtnum,Year,Mon,Day
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

                if scanName =='MS':
                    a2dprx = ma.masked_less(a2dprx,0)-12
                    a2dprx = ma.masked_outside(a2dprx, 0, 24).filled(-9999)

                a1dprx    = a2dprx.flatten()
                a1dpry    = a2dpry.flatten()

                #-- Search DPR granules --
                if radar in ['Ku','Ka']:
                    prdLev = '2A'
                elif radar in ['DPRGMI']:
                    prdLev = '2B'

                dprbaseDir   = '/work/hk01/PMM/NASA/GPM.%s/%s/V%s'%(radar, prdLev, verDPR)
                if  (radar in ['Ku','Ka','DPRGMI']):
                    srcDir = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch= srcDir + '/*.%s.V%s.HDF5'%(obtnum,fullverDPR)
                    srcPath= glob.glob(ssearch)[0]
                else:
                    print 'check radar',radar
                    sys.exit()

                #-- Read DPR variables --
                if   varName in ['zFactorCorrected','zFactorMeasured']:
                    with h5py.File(srcPath) as h:
                        Dat0 = h[var][:,:,-nvect*2:]

                    #-- Convert dBZ --> Z ---
                    '''
                    #XdBZ = 10*log_10(Z/Z0)
                    #Z    = Z0* 10^(X*0.1)
                    #No need to multiply Z0.
                    #(Because they are devided by Z0 when converted to dBZ at the end)
                    '''
                    DatLN = np.power(10, 0.1*ma.masked_less_equal(Dat0,-9999.9)).filled(-9999.9)

                    #-- Average vertical ranges (over Linearlized Z)--
                    '''
                    # original: vertical 176bins (with 125m resol)
                    # 125m ->250m (total 40bin)
                    '''
                    DatLN = average_2ranges_3d(DatLN, miss=-9999.9, dtype=float32, fill=False)

                    #-- Average 9 grids (over Linearlized Z)--
                    a2datTmp = ave_9grids_3d(DatLN, a1dpry, a1dprx, miss=-9999.9)

                    #-- Convert to dBZ --
                    Dat = 10*np.log10(ma.masked_equal(a2datTmp, -9999.9))

                    #-- Scale by 100, with int16 ---
                    Dat = (100*ma.masked_invalid(Dat)).filled(-9999.9).astype(dattype[var])


                elif   varName =='precipRate':
                    with h5py.File(srcPath) as h:
                        Dat0 = h[var][:,:,-nvect*2:]

                    #-- Average vertical ranges (over Linearlized Z)--
                    DatCnc = average_2ranges_3d(Dat0, miss=-9999.9, dtype=float32, fill=False)

                    #-- Average 9 grids (over Linearlized Z)--

                    Dat     = ave_9grids_3d(DatCnc, a1dpry, a1dprx, miss=-9999.9)

                elif varName in ['precipTotWaterCont']:
                    with h5py.File(srcPath) as h:
                        Dat0 = h[var][:,:,-nvect:]

                    #-- Average 9 grids --

                    Dat     = ave_9grids_3d(Dat0, a1dpry, a1dprx, miss=-9999.9)

                elif varName in ['precipRateESurface','surfPrecipTotRate','precipRateNearSurface']:
                    with h5py.File(srcPath) as h:
                        Dat0 = h[var][:]

                    Dat     = ave_9grids_2d(Dat0, a1dpry, a1dprx, miss=-9999.9)


                elif varName in ['elevation']:
                    with h5py.File(srcPath) as h:
                        Dat0 = (h[var][:]).astype(dattype[var])

                    a1dpryTmp = ma.masked_less(a1dpry,0).filled(0)
                    a1dprxTmp = ma.masked_less(a1dprx,0).filled(0)
                    Dat = Dat0[a1dpryTmp, a1dprxTmp]
                    if dattype[var] in ['int32','int16']:
                        miss = -9999
                    elif dattype[var] in ['float32']:
                        miss = -9999.9
                    Dat = ma.masked_where(a1dpry==-9999, Dat).filled(miss)
   
                elif varName in ['typePrecip']:
                    with h5py.File(srcPath) as h:
                        Dat0 = (h[var][:]/10000000).astype(dattype[var])
                    DatStrat = ma.masked_not_equal(Dat0,1).filled(0)
                    DatConv  = ma.masked_not_equal(Dat0,2).filled(0)/2
                    DatOther = ma.masked_not_equal(Dat0,3).filled(0)/3

                    DatStrat = sum_9grids_2d(DatStrat, a1dpry, a1dprx, miss=0).filled(0)
                    DatConv  = sum_9grids_2d(DatConv,  a1dpry, a1dprx, miss=0).filled(0)
                    DatOther = sum_9grids_2d(DatOther, a1dpry, a1dprx, miss=0).filled(0)
                    Dat = (DatStrat + DatConv*10 + DatOther*100)


                elif varName in ['heightStormTop']:
                    with h5py.File(srcPath) as h:
                        Dat0 = ma.masked_greater(h[var][:], 32767).filled(32767).astype(dattype[var])

                    a1dpryTmp = ma.masked_less(a1dpry,0).filled(0)
                    a1dprxTmp = ma.masked_less(a1dprx,0).filled(0)
                    Dat = Dat0[a1dpryTmp, a1dprxTmp]
                    if dattype[var] in ['int32','int16']:
                        miss = -9999
                    elif dattype[var] in ['float32']:
                        miss = -9999.9
                    Dat = ma.masked_where(a1dpry==-9999, Dat).filled(miss)
                    print Dat.min(),Dat.max()

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
        
                        outDir  = outbaseDir + '/%s/%04d%02d'%(varNameOut,Year,Mon)
                        outPath = outDir + '/%s.%05d.%d.npy'%(varNameOut,epcid, dnum[epcid])
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
                    outDir  = outbaseDir + '/%s/%04d%02d'%(varNameOut,Year,Mon)
                    outPath = outDir + '/%s.%05d.%d.npy'%(varNameOut,epcid,dnum[epcid])
                    util.mk_dir(outDir)
    
                    np.save(outPath, array(dstack[epcid]).astype(dattype[var]))
    
                    del dstack[epcid]
                    gc.collect()
    
    
            #--- Joint segments ---
            for epcid in range(epcid_min,epcid_max):
                outDir  = outbaseDir + '/%s/%04d%02d'%(varNameOut,Year,Mon)
                ssearch = outDir + '/%s.%05d.*.npy'%(varNameOut,epcid)
                lsrcPath= sort(glob.glob(ssearch))
    
                if len(lsrcPath)==0: continue
    
                aout = deque([])
                for srcPath in lsrcPath:
                    atmp = np.load(srcPath)
                    aout.extend(atmp)
    
                outPath =  outDir + '/%s.%05d.npy'%(varNameOut,epcid)
    
                np.save(outPath, aout)
                print outPath
    
                #-- delete temporary files --
                for srcPath in lsrcPath:
                    os.remove(srcPath)
    


