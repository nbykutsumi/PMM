from numpy import *
import h5py
import calendar
import myfunc.util as util
import random, glob
import numpy as np
import sys


prod     = '1C'
verGMI   = '05'
subverGMI= 'A'
fullverGMI='%s%s'%(verGMI, subverGMI)

iYM = [2017,1]
eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
nrand= 10

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)
nsample = 5

verDPR    = '06'
subverDPR = 'A'
fullverDPR= '%s%s'%(verDPR, subverDPR)

#prodDPR = '2A'
prodDPR = '2B'
#radar= 'Ku'
#radar= 'Ka'
radar= 'DPRGMI'

#gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
#gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
dprbaseDir   = '/work/hk01/PMM/NASA/GPM.%s/%s/V%s'%(radar, prodDPR, verDPR)

dbbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)

def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''

    if ma.is_masked(a3in):
        a3in = a3in.filled(miss)  # 2019/12/02
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

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
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



for Year,Mon in lYM:
    #-- Read granule list --
    listPath = '/work/hk01/utsumi/PMM/EPCDB/list/list.1C.V05.%04d%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines=f.readlines(); f.close()
    dymd = {}
    for line in lines:
        line = map(int, line.strip().split(','))
        oid, yyyy,mm,dd,itime,etime = line
        dymd[oid] = [yyyy,mm,dd]


    #-- make epcid list --
    ssearch = dbbaseDir + '/pYXpmw/%04d%02d/pYXpmw.*.npy'%(Year,Mon)
    lyxPath = glob.glob(ssearch)
    lyxPath = sort(lyxPath)
    #lipath  = random.sample(range(len(lyxPath)), k=5)
    lipath  = range(4)
    lipath  = sort(lipath)

    for ipath in lipath:
        yxPath = lyxPath[ipath]
        #print yxPath
        epcid  = yxPath.split('/')[-1].split('.')[1]
        #if epcid != '02994':
        #    continue
        print 'epcid=',epcid

        
        gnumPath = dbbaseDir + '/gNum/%04d%02d/gNum.%s.npy'%(Year,Mon, epcid)
        #mdhmsPath= dbbaseDir + '/mdhms/%04d%02d/mdhms.%s.npy'%(Year,Mon, epcid)
        #yearPath = dbbaseDir + '/Year/%04d%02d/Year.%s.npy'%(Year,Mon, epcid)
        #esurfPath= dbbaseDir + '/precipRateESurface/%04d%02d/precipRateESurface.%s.npy'%(Year,Mon, epcid)

        #surfPath= dbbaseDir + '/DPRGMI_NS_surfPrecipTotRate/%04d%02d/DPRGMI_NS_surfPrecipTotRate.%s.npy'%(Year,Mon, epcid)

        watcontPath= dbbaseDir + '/DPRGMI_NS_precipTotWaterCont/%04d%02d/DPRGMI_NS_precipTotWaterCont.%s.npy'%(Year,Mon, epcid)
        zmkaPath = dbbaseDir + '/Ka_MS_zFactorMeasured/%04d%02d/Ka_MS_zFactorMeasured.%s.npy'%(Year,Mon, epcid)
        #zmkuPath = dbbaseDir + '/Ku_NS_zFactorMeasured/%04d%02d/Ku_NS_zFactorMeasured.%s.npy'%(Year,Mon, epcid)

        #-- Read files --
        a1gnum   = np.load(gnumPath)
        a2yx     = np.load(yxPath)
        #a2zmka   = np.load(zmkaPath)
        #a2zmku   = np.load(zmkuPath)
        a2watcont= np.load(watcontPath)

        #a2mdhms  = np.load(mdhmsPath)
        #a2year   = np.load(yearPath)
        #a1surf   = np.load(surfPath)

        #-- Loop in a single EPC subset database
        if len(a1gnum)>nsample:
            ksample = nsample
        else:
            ksample = len(a1gnum)
        lient  = random.sample(range(len(a1gnum)), k=ksample)
        lient  = sort(lient)
      
        for ient in lient:
            gnum = a1gnum[ient]
            gmiy,gmix  = a2yx[ient]
            #zmka1  = a2zmka[ient]
            #zmku1  = a2zmku[ient]
            watcont1  = a2watcont[ient]
            #Mon1,Day1,Hour1,Mnt1,Sec1 = a2mdhms[ient]
            #Year1= a2year[ient]
            #surf1 = a1surf[ient]
            Year1,Mon1,Day1 = dymd[gnum]

            #-- location matching file (GMI-->DPR) --
            dpridxDir= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year1,Mon1, Day1)
            dprxPath= dpridxDir + '/Xpy.1.%06d.npy'%(gnum)
            dpryPath= dpridxDir + '/Ypy.1.%06d.npy'%(gnum)
    
            a2dprx  = np.load(dprxPath)
            a2dpry  = np.load(dpryPath)

            if radar=='Ka':
                a2dprx = (ma.masked_less(a2dprx,0) - 12)
                a2dprx = ma.masked_outside(a2dprx, 0,24).filled(-9999)
            #----------------------------------------


            dpry = a2dpry[gmiy, gmix-83]
            dprx = a2dprx[gmiy, gmix-83]
            

            #-- Read DPR granule HDF file --
            dprDir = dprbaseDir + '/%04d/%02d/%02d'%(Year1,Mon1,Day1)
            ssearch= dprDir + '/*.%06d.V???.HDF5'%(gnum)
            dprPath= glob.glob(ssearch)[0]

            with h5py.File(dprPath) as h:
                #a3zmka2 = h['/MS/PRE/zFactorMeasured'][:,:,-60*2:]
                #a3zmku2 = h['/NS/PRE/zFactorMeasured'][:,:,-60*2:]
                a3watcont2 = h['/NS/precipTotWaterCont'][:,:,-60:]
    
                #a2hour2 = h['/NS/ScanTime/Hour'][:]
                #a2mnt2  = h['/NS/ScanTime/Minute'][:]
                #a2sec2  = h['/NS/ScanTime/Second'][:]
                #a2esurf2= h['/NS/SLV/precipRateESurface'][:]
                #a2surf2= h['/NS/surfPrecipTotRate'][:]
    
                #Hour2   = a2hour2[dpry]
                #Mnt2    = a2mnt2[dpry]
                #Sec2    = a2sec2[dpry]

            #-- average 9-grids --
            #zmka2 = a3zmka2[dpry-1:dpry+1+1,dprx-1:dprx+1+1]
            #zmku2 = a3zmku2[dpry-1:dpry+1+1,dprx-1:dprx+1+1]
            watcont2 = ave_9grids_3d(a3watcont2, array([dpry]), array([dprx]), miss=-9999.9).astype(float32)[0]
            #esurf2  = ave_9grids_2d(a2esurf2, array([dpry]), array([dprx]), miss=-9999.9)[0]
            #surf2  = ave_9grids_2d(a2surf2, array([dpry]), array([dprx]), miss=-9999.9)[0]
 

            ##****** Only Zm ****************
            ##-- Convert dBZ --> Z ---
            #'''
            ##XdBZ = 10*log_10(Z/Z0)
            ##Z    = Z0* 10^(X*0.1)
            ##No need to multiply Z0.
            ##(Because they are devided by Z0 when converted to dBZ at the end)
            #'''
            #Dat0  = zmka2
            ##Dat0  = zmku2
            #DatLN = np.power(10, 0.1*ma.masked_less_equal(Dat0,-9999.9)).filled(-9999.9)

            ##-- Average vertical ranges (over Linearlized Z)--
            #'''
            ## original: vertical 176bins (with 125m resol)
            ## 0-10km  (-80:last bin): 125m ->250m (total 40bin)
            ## 10-15km (-120:-80 bin): 125m ->500m (total 10bin)
            #'''
            #DatLN = average_2ranges_3d(DatLN,miss=-9999.9, dtype=float32, fill=False)

            ##-- Average 9 grids (over Linearlized Z)--
            #a2datTmp = ave_9grids_3d(DatLN, array([1]), array([1]), miss=-9999.9)

            ##-- Convert to dBZ --
            #Dat = 10*np.log10(ma.masked_equal(a2datTmp, -9999.9))

            ##-- Scale by 100, with int16 ---
            #Dat = (100*ma.masked_invalid(Dat)).filled(-9999.9).astype(int16)
            ##**********************************


            #print '-'*10
            #print zmka1[-1:]
            ##print zmku1[-1:]
            #print ''
            #print Dat[0,-1:]
            ##print zmka1.shape
            #print zmka2.shape
            #print zmku1.shape
            #print zmku2.shape
            #print Dat.shape
            print watcont1
            print watcont2
            print watcont1-watcont2
            #print 'Hour  ',Hour1,Hour2
            #print 'Mnt   ',Mnt1, Mnt2
            #print 'Sec   ',Sec1, Sec2
            #print 'esurf ',esurf, esurf2
            #print 'surf ',surf1, surf2
            print '-'*10
            
        
