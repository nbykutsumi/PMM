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

iYM = [2017,2]
eYM = [2017,2]
lYM = util.ret_lYM(iYM,eYM)
nrand= 10

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)
nsample = 5

verDPR    = '06'
subverDPR = 'A'
fullverDPR= '%s%s'%(verDPR, subverDPR)

radar= 'Ku'

#gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
#gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
dprbaseDir   = '/work/hk01/PMM/NASA/GPM.%s/2A/V%s'%(radar, verDPR)

dbbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)

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


for Year,Mon in lYM:
    #-- make epcid list --
    ssearch = dbbaseDir + '/pYXpmw/%04d%02d/pYXpmw.*.npy'%(Year,Mon)
    lyxPath = glob.glob(ssearch)
    lipath  = random.sample(range(len(lyxPath)), k=5)
    lipath  = sort(lipath)

    for ipath in lipath:
        yxPath = lyxPath[ipath]
        #print yxPath
        epcid  = yxPath.split('/')[-1].split('.')[1]
        print 'epcid=',epcid

        gnumPath = dbbaseDir + '/gNum/%04d%02d/gNum.%s.npy'%(Year,Mon, epcid)
        mdhmsPath= dbbaseDir + '/mdhms/%04d%02d/mdhms.%s.npy'%(Year,Mon, epcid)
        yearPath = dbbaseDir + '/Year/%04d%02d/Year.%s.npy'%(Year,Mon, epcid)
        esurfPath= dbbaseDir + '/precipRateESurface/%04d%02d/precipRateESurface.%s.npy'%(Year,Mon, epcid)

        #-- Read files --
        a1gnum   = np.load(gnumPath)
        a2yx     = np.load(yxPath)
        a2mdhms  = np.load(mdhmsPath)
        a2year   = np.load(yearPath)
        a2esurf  = np.load(esurfPath)


        #-- Loop in a single EPC subset database
        if len(a1gnum)>nsample:
            ksample = nsample
        else:
            ksample = len(a1gnum)
        lient  = random.sample(range(len(a1gnum)), k=ksample)
        lient  = sort(lient)
      
        for ient in lient:
            gnum = a1gnum[ient]
            Mon1,Day1,Hour1,Mnt1,Sec1 = a2mdhms[ient]
            Year1= a2year[ient]
            gmiy,gmix  = a2yx[ient]
            esurf= a2esurf[ient]

            #-- location matching file (GMI-->DPR) --
            dpridxDir= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year1,Mon1, Day1)
            dprxPath= dpridxDir + '/Xpy.1.%06d.npy'%(gnum)
            dpryPath= dpridxDir + '/Ypy.1.%06d.npy'%(gnum)
    
            a2dprx  = np.load(dprxPath)
            a2dpry  = np.load(dpryPath)
            #----------------------------------------


            dpry = a2dpry[gmiy, gmix-83]
            dprx = a2dprx[gmiy, gmix-83]
            

            #-- Read GMI granule HDF file --
            dprDir = dprbaseDir + '/%04d/%02d/%02d'%(Year1,Mon1,Day1)
            ssearch= dprDir + '/*.%06d.V???.HDF5'%(gnum)
            dprPath= glob.glob(ssearch)[0]

            h = h5py.File(dprPath)
            a2hour2 = h['/NS/ScanTime/Hour'][:]
            a2mnt2  = h['/NS/ScanTime/Minute'][:]
            a2sec2  = h['/NS/ScanTime/Second'][:]
            a2esurf2= h['/NS/SLV/precipRateESurface'][:]

            Hour2   = a2hour2[dpry]
            Mnt2    = a2mnt2[dpry]
            Sec2    = a2sec2[dpry]

            #-- average 9-grids --
            esurf2  = ave_9grids_2d(a2esurf2, array([dpry]), array([dprx]), miss=-9999.9)[0]
 

            print '-'*10
            print 'Hour  ',Hour1,Hour2
            print 'Mnt   ',Mnt1, Mnt2
            print 'Sec   ',Sec1, Sec2
            print 'esurf ',esurf, esurf2
            print '-'*10
            
        
