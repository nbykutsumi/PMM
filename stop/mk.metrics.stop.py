import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import random

myhost = socket.gethostname()
if myhost == 'shui':
    listDir    = '/work/hk01/utsumi/PMM/US/obtlist'
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    epcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'

    figDir   = '/home/utsumi/temp/stop'

elif myhost == 'well':
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()

lsurftype = ['ocean','vegetation','coast','snow']
#lsurftype = ['vegetation']

dsurflabel={ 'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

calcflag= True
#calcflag= False
#region = 'US'
region = 'GLB'
#lseason = ['JJADJF']
#lseason = ['JJA']
lseason = [8]
DB_MAXREC = 10000
DB_MINREC = 1000

lrettype = ['ml','NScmb']
#lrettype = ['ml']
#lrettype = ['NScmb']
prmin = 0.2
#prmin = 1
#prmin = 0.01
#prmin = 0.0
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
stopexpr = 'best01'
#stopact  = 'LTQZ'
stopact  = 'HTQZ'
nsample  = 100 # number of sampled orbits per month (for region=='GLB')
def read_orbitlist_us(Year,Mon):
    listPath = listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath); lines=f.readlines(); f.close()
    lout=[]
    for line in lines:
        line = map(int,line.strip().split(','))
        lout.append(line)
    return lout

def read_orbitlist_glob(Year,Mon):
    iDay = 1
    eDay = calendar.monthrange(Year,Mon)[1]
    lout = []
    for Day in range(iDay,eDay+1):
        retDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
        ssearch= retDir + '/nsurf%s.??????.y-9999--9999.nrec%d.npy'%('NScmb', DB_MAXREC)
        lsrcPath = np.sort(glob.glob(ssearch))
        for srcPath in lsrcPath:
            oid = int(srcPath.split('.')[-4])

            lout.append([Year,Mon,Day,oid,-9999,-9999])
    return lout

def ret_lym(season):
    if season=='JJADJF':
        lYM = util.ret_lYM([2014,6],[2014,8]) + util.ret_lYM([2014,12],[2015,2])
    if season=='ALL':
        lYM = util.ret_lYM([2014,6],[2015,5])
    elif season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
    elif season=='AD':
        lYM = [[2014,8],[2014,12]]


    elif type(season)==int:
        if season <6:
            lYM = [[2015,season]]
        else:
            lYM = [[2014,season]]
    elif season =='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    elif season =='DJ':
        lYM = util.ret_lYM([2014,12],[2015,1])

    else:
        print 'check season',season
        sys.exit()
    return lYM

def ave_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
    #-- Average 9 grids --
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

    #------------
    a1datTmp = ma.masked_equal(a2datTmp,miss).mean(axis=0)
    a1datTmp[a1dprmask] = miss


    return a1datTmp
#*************************************************


for season in lseason:
    lYM = ret_lym(season)
    for rettype in lrettype:
        if rettype =='ml':
            exprTmp = '%s-%s'%(stopexpr, stopact)
        elif rettype in ['NS','MS','NScmb','MScmb']:
            exprTmp = expr

        for surftype in lsurftype:
    
            for Year,Mon in lYM:
                if calcflag==False: continue

                a1ret = array([])
                a1obs = array([])


                print surftype,Year,Mon
                if region=='US':
                    lorbit = read_orbitlist_us(Year,Mon)
                elif region=='GLB':
                    lorbit = read_orbitlist_glob(Year,Mon)
                    random.seed(0)
                    lorbit = random.sample(lorbit, nsample)
                else:
                    print 'check region',region
                    sys.exit()

                lorbit.sort()
                for orbit in lorbit:
                    Day,oid,iy,ey = orbit[2:]
                    print Year,Mon,Day ,oid

                    #-- Read surface type ----
                    ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.20140701-S030429-E043702.001921.V05A.HDF5'
                    ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.*.%06d.V05A.HDF5'%(Year,Mon,Day,oid)
                    lgprofPath = glob.glob(ssearch)
                    if len(lgprofPath)==0:
                        print 'No file for',Year,Mon,Day,oid
                        print ssearch
                        continue
                    gprofPath= lgprofPath[0]

                    with h5py.File(gprofPath,'r') as h:
                        a2surftype = h['S1/surfaceTypeIndex'][:]

                    #-- Read EPC precip (for Rain/No Rain) ----
                    epcDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
                    a2prcp = np.load(epcDir + '/nsurf%s.%06d.y-9999--9999.nrec%d.npy'%('NScmb', oid, DB_MAXREC))

                    #-- Read storm top ----
                    if rettype =='ml':
                        srcbaseDir = tankbaseDir + '/utsumi/PMM/stop/orbit/%s-%s-ssn%s'%(stopexpr, stopact, 0)
                        srcDir     = srcbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                        srcPath = srcDir + '/stop.%06d.npy'%(oid)
                        a2ret   = np.load(srcPath)

                    elif rettype in ['NS','MS','NScmb','MScmb']:
                        epcDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(exprTmp,Year,Mon,Day)
                        a2ret  = np.load(epcDir + '/top-heightStormTopNS.%06d.y-9999--9999.nrec%d.npy'%(oid, DB_MAXREC))

                    #-- Read DPR-Combined -----------------------------------------------
                    dprbaseDir = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06'
                    dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch = dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)

                    try:
                        dprPath = glob.glob(ssearch)[0]
                    except:
                        print 'No DPR file for oid=',oid
                        continue
                
                    with h5py.File(dprPath, 'r') as h:
                        a2obsOrg = h['NS/PRE/heightStormTop'][:]

                    #-- index file for colocation ---
                    xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
                    xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
                    yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)
                
                    if not os.path.exists(xPath):
                        continue
                    a2x = np.load(xPath)
                    a2y = np.load(yPath)

                    ny,nx= a2x.shape
                    a1x = a2x.flatten()
                    a1y = a2y.flatten()

                    a2obsOrg = ma.masked_less(a2obsOrg,0).filled(-9999.)
                    a2obs = ave_9grids_2d(a2obsOrg, a1y, a1x, miss=-9999.).filled(-9999.).reshape(ny,nx)


                    #-- Extract regional data ---------------------
                    if iy >0:
                        a2surftype = a2surftype[iy:ey+1]
                        a2prcp= a2prcp[iy:ey+1]
                        a2ret = a2ret[iy:ey+1]
                        a2obs = a2obs[iy:ey+1]  # radar

                    
                    #-- Extract x=83-137 of PMW  ------------------
                    a2surftype = a2surftype[:,83:137+1]
                    a2prcp= a2prcp[:,83:137+1]
                    a2ret = a2ret[:,83:137+1] 
 
                    if a2ret.shape[0]==0: continue 

                    #-- Screen by surface types -------------------
                    if surftype=='ocean':
                        a2masksurf = ma.masked_not_equal(a2surftype,1).mask
                    elif surftype=='vegetation':
                        a2masksurf = ma.masked_outside(a2surftype,3,7).mask
                    elif surftype=='snow':
                        a2masksurf = ma.masked_outside(a2surftype,8,11).mask
                    elif surftype=='coast':
                        a2masksurf = ma.masked_not_equal(a2surftype,13).mask
                    else:
                        print '\n'+'check surftype',surftype
                        sys.exit() 
                    #-- Screeen missing pixels both datasets--
                    a2mask1 = ma.masked_less_equal(a2ret, 0).mask
                    a2mask2 = ma.masked_less_equal(a2obs,0).mask
                    a2mask  = a2mask1 + a2mask2
                    a2mask  = a2mask + a2masksurf

                    #-- Screeen No precip pixels --
                    a2maskP = ma.masked_less(a2prcp,prmin).mask
                    a2mask  = a2mask + a2maskP
                    #------------------------------

                    a1retTmp = ma.masked_where(a2mask, a2ret).compressed()
                    a1obsTmp = ma.masked_where(a2mask, a2obs).compressed()
    
                    a1ret = np.concatenate([a1ret, a1retTmp])
                    a1obs = np.concatenate([a1obs, a1obsTmp])

                 

                #*******************************
                # Save
                #-------------------------------
                stamp =  '%s.%s.%s.%04d.%02d'%(exprTmp, region, surftype, Year, Mon)

                pickleDir  = '/home/utsumi/temp/stop/pickle'
                util.mk_dir(pickleDir)

                retPath = pickleDir + '/stop.%s.ret.%s.npy'%(rettype,stamp)
                obsPath = pickleDir + '/stop.%s.obs.%s.npy'%(rettype,stamp)
                if calcflag == True:
                    np.save(retPath, a1ret) 
                    np.save(obsPath, a1obs)                    
            

                print surftype,a1ret.shape 
            #*******************************
            # Load
            #-------------------------------
            a1ret = np.array([])
            a1obs = np.array([])
            for Year,Mon in lYM:
                stamp =  '%s.%s.%s.%04d.%02d'%(exprTmp, region, surftype, Year, Mon)
                pickleDir  = '/home/utsumi/temp/stop/pickle'

                retPathTmp = pickleDir + '/stop.%s.ret.%s.npy'%(rettype,stamp)
                obsPathTmp = pickleDir + '/stop.%s.obs.%s.npy'%(rettype,stamp)

                a1ret = np.concatenate([a1ret, np.load(retPathTmp).astype(float32)]) 
                a1obs = np.concatenate([a1obs, np.load(obsPathTmp).astype(float32)]) 

            #*******************************
            # Metrics
            #-------------------------------
            cc = np.ma.corrcoef(a1ret,a1obs, allow_masked=True)[0][1] 
            nbias= (a1ret - a1obs).mean()/a1obs.mean()
            rmse = np.sqrt(((a1ret-a1obs)**2).mean())

            stampOut = '%s.%s.%s.%s'%(exprTmp, region, surftype, season)
            outPath = figDir + '/metrics.stop.%s.%s.csv'%(rettype, stampOut)
            lout = [['cc','nbias','rmse']]
            lout.append([cc,nbias,rmse])
            sout = util.list2csv(lout)
            f=open(outPath,'w'); f.write(sout); f.close()
            print ''
            print outPath
            print 'cc=',cc
            print 'bias=',nbias
            print 'rmse=',rmse

