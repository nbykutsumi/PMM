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
NDB = 10000

iYM = [2017,1]
eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
nrand= 10

cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)
nsample = 5
gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'


dbbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)
varNameFull = ''
varName     = varNameFull.split('/')[-1]

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
    lipath  = random.sample(range(len(lyxPath)), k=5)
    lipath  = sort(lipath)

    for ipath in lipath:
        yxPath = lyxPath[ipath]
        #print yxPath
        epcid  = yxPath.split('/')[-1].split('.')[1]
        print epcid

        gnumPath = dbbaseDir + '/gNum/%04d%02d/gNum.%s.npy'%(Year,Mon, epcid)
        #mdhmsPath= dbbaseDir + '/mdhms/%04d%02d/mdhms.%s.npy'%(Year,Mon, epcid)
        #yearPath = dbbaseDir + '/Year/%04d%02d/Year.%s.npy'%(Year,Mon, epcid)
        tcPath  = dbbaseDir + '/%s/%04d%02d/%s.%s.npy'%('Tc', Year,Mon, 'Tc', epcid)
        #varPath  = dbbaseDir + '/%s/%04d%02d/%s.%s.npy'%(varName, Year,Mon, epcid, varName)
        #-- Read files --
        a1gnum   = np.load(gnumPath)
        a2yx     = np.load(yxPath)
        #a2mdhms  = np.load(mdhmsPath)
        #a2year   = np.load(yearPath)
        a2tc     = np.load(tcPath)
        #avar     = np.load(varPath)
        #-- Loop in a single EPC subset database
        if len(a1gnum)>nsample:
            ksample = nsample
        else:
            ksample = len(a1gnum)
        lient  = random.sample(range(len(a1gnum)), k=ksample)
        lient  = sort(lient)
      
        #print 'a2yx' 
        #print a2yx
        #sys.exit() 
        for ient in lient:
            gnum = a1gnum[ient]
            y,x  = a2yx[ient]
            tc1  = a2tc[ient]
            #Mon1,Day1,Hour1,Mnt1,Sec1 = a2mdhms[ient]
            #Year1= a2year[ient]
            #vdat = avar[ient]

            Year1,Mon1,Day1 = dymd[gnum] 
            #-- Read GMI matched data --
            matchbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
            tcDir1 = matchbaseDir + '/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d'%(Year1,Mon1,Day1)
            tcDir2 = matchbaseDir + '/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d'%(Year1,Mon1,Day1)
            tcPath1= tcDir1 + '/Tc.%06d.npy'%(gnum)
            tcPath2= tcDir2 + '/TcS2.1.%06d.npy'%(gnum)
            a3tc2  = np.concatenate([ np.load(tcPath1), np.load(tcPath2)], axis=2)
            tc2    = a3tc2[y,x-103,:]
            ##-- Read GMI granule HDF file --
            #gmiDir = gmibaseDir + '/%04d/%02d/%02d'%(Year1,Mon1,Day1)
            #ssearch= gmiDir + '/*.%06d.V???.HDF5'%(gnum)
            #gmiPath= glob.glob(ssearch)[0]

            #with h5py.File(gmiPath) as h:
            #   a2hour2 = h['/S1/ScanTime/Hour'][:]
            #   a2mnt2  = h['/S1/ScanTime/Minute'][:]
            #   a2sec2  = h['/S1/ScanTime/Second'][:]
            #   avar2   = h[varNameFull][:]

            #Hour2   = a2hour2[y]
            #Mnt2    = a2mnt2[y]
            #Sec2    = a2sec2[y]


            print '-'*10
            print tc1
            print tc2
            print tc1==tc2 
