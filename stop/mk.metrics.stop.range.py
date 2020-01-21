#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
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
from sklearn import linear_model
from scipy import stats

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

#lsurftype = ['ocean','vege','coast','snow']
#lsurftype = ['vege']
lsurftype = ['ocean','seaice','vege','snow','swater','coast','siedge']

#calcflag= True
calcflag= False
#region = 'US'
region = 'GLB'
iYM = [2014,6]
eYM = [2015,5]
#iYM = [2017,1]
#eYM = [2017,12]
lYMorg = util.ret_lYM(iYM,eYM)
#lseason = ['ALL']
lseason = ['JJA']
rnrtype = 'epc'
#rnrtype = 'dpr'

lrettype = ['ml','NScmb']
#lrettype = ['ml-2017','NScmb']
#lrettype = ['ml']
#lrettype = ['NScmb']
prmin = 0.2
#prmin = 1
#prmin = 0.01
#prmin = 0.0

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
stopexpr = 'best01'
#stopexpr = 'best01cr'
#stopact  = 'LTQZ'
stopact  = 'HTQZ'
nsample  = 100 # number of sampled orbits per month (for region=='GLB')


lminmax = [[0,99],[0,4],[4,8],[8,15]]


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
        srcDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch= srcDir + '/1C.GPM.GMI.*.??????.V05A.HDF5'
        lsrcPath = np.sort(glob.glob(ssearch))
        for srcPath in lsrcPath:
            oid = int(srcPath.split('.')[-3])

            lout.append([Year,Mon,Day,oid,-9999,-9999])
    return lout

def read_orbitlist_testdata(Year,Mon):
    lDTime  = np.load('/home/utsumi/bin/PMM/stop/ldtime-%04d-test.npy'%(Year), allow_pickle=True)
    lDTime  = np.sort(lDTime)
    lout = []
    for DTime in lDTime:
        if DTime.month != Mon: continue
        Day = DTime.day
        srcDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch= srcDir + '/1C.GPM.GMI.*.??????.V05A.HDF5'
        lsrcPath = np.sort(glob.glob(ssearch))
        for srcPath in lsrcPath:
            oid = int(srcPath.split('.')[-3])

            lout.append([Year,Mon,Day,oid,-9999,-9999])

    return lout



def ret_lmon(season):
    if season=='JJADJF':
        lmon = [6,7,8,12,1,2]
    if season=='ALL':
        lmon = range(1,12+1)
    elif season=='JJA':
        lmon = [6,7,8]
    elif season=='SON':
        lmon = [9,10,11]
    elif season=='DJF':
        lmon = [12,1,2]
    elif season=='MAM':
        lmon = [3,4,5]
    elif season=='AD':
        lmon = [8,12]
    elif season =='JJ':
        lmon = [6,7]
    elif season =='DJ':
        lmon = [12,1]
    elif type(season)==int:
        lmon = [int(season)]
    else:
        print 'check season',season
        sys.exit()
    return lmon


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
    lmon = ret_lmon(season)
    lYM = [[Y,M] for [Y,M] in lYMorg if M in lmon]
    for rettype in lrettype:
        if rettype in ['ml', 'ml-2017']:
            exprTmp = '%s-%s'%(stopexpr, stopact)

        elif rettype in ['NS','MS','NScmb','MScmb']:
            exprTmp = expr

        for surftype in lsurftype:
    
            for Year,Mon in lYM:
                if Mon not in lmon: continue
                if calcflag==False: continue

                a1ret = array([])
                a1obs = array([])

                if region=='US':
                    lorbit = read_orbitlist_us(Year,Mon)

                elif (region=='GLB')and(rettype =='ml-%04d'%(Year)):
                    lorbit = read_orbitlist_testdata(Year,Mon)

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

                    #-- Read Rain/No Rain data ----
                    if rnrtype=='epc':
                        epcDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
                        a2rnr = np.load(epcDir + '/nsurf%s.%06d.y-9999--9999.nrec%d.npy'%('NScmb', oid, DB_MAXREC))

                    elif rnrtype=='dpr':
                        rnrDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.DPRGMI.V06A.9ave.surfPrecipTotRate/%04d/%02d/%02d'%(Year,Mon,Day)
                        rnrPath = rnrDir + '/surfPrecipTotRate.%06d.npy'%(oid)
                        a2rnr = np.load(rnrPath)

                    #-- Read storm top ----
                    if rettype in ['ml','ml-2017']:
                        srcbaseDir = tankbaseDir + '/utsumi/PMM/stop/orbit/%s-%s-ssn%s'%(stopexpr, stopact, 0)
                        srcDir     = srcbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                        srcPath = srcDir + '/stop.%06d.npy'%(oid)
                        a2ret   = np.load(srcPath)

                    elif rettype in ['NS','MS','NScmb','MScmb']:
                        epcDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(exprTmp,Year,Mon,Day)
                        a2ret  = np.load(epcDir + '/top-heightStormTopNS.%06d.y-9999--9999.nrec%d.npy'%(oid, DB_MAXREC))

                    #-- Read DPR-Ku -----------------------------------------------
                    dprbaseDir = workbaseDir + '/hk02/PMM/NASA/GPM.Ku/2A/V06'
                    dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch = dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)

                    try:
                        dprPath = glob.glob(ssearch)[0]
                    except:
                        print 'No DPR file for oid=',oid
                        print ssearch
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
                        a2rnr = a2rnr[iy:ey+1]
                        a2ret = a2ret[iy:ey+1]
                        a2obs = a2obs[iy:ey+1]  # radar

                    
                    #-- Extract x=83-137 of PMW  ------------------
                    a2surftype = a2surftype[:,83:137+1]
                    a2ret = a2ret[:,83:137+1] 
                    if rnrtype=='epc':
                        a2rnr = a2rnr[:,83:137+1]

 
                    if a2ret.shape[0]==0: continue 

                    #-- Screen by surface types -------------------
                    if surftype=='ocean':
                        a2masksurf = ma.masked_not_equal(a2surftype,1).mask
                    elif surftype=='seaice':
                        a2masksurf = ma.masked_not_equal(a2surftype,2).mask
                    elif surftype=='vege':
                        a2masksurf = ma.masked_outside(a2surftype,3,7).mask
                    elif surftype=='snow':
                        a2masksurf = ma.masked_outside(a2surftype,8,11).mask
                    elif surftype=='swater':
                        a2masksurf = ma.masked_not_equal(a2surftype,12).mask
                    elif surftype=='coast':
                        a2masksurf = ma.masked_not_equal(a2surftype,13).mask
                    elif surftype=='siedge':
                        a2masksurf = ma.masked_not_equal(a2surftype,14).mask
                    else:
                        print '\n'+'check surftype',surftype
                        sys.exit() 
                    #-- Screeen missing pixels both datasets--
                    a2mask1 = ma.masked_less_equal(a2ret, 0).mask
                    a2mask2 = ma.masked_less_equal(a2obs,0).mask
                    a2mask  = a2mask1 + a2mask2
                    a2mask  = a2mask + a2masksurf

                    #-- Screeen No precip pixels --
                    a2maskP = ma.masked_less(a2rnr,prmin).mask

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
            # Linear regression
            #-------------------------------
            slope, intercept, r_value, p_value, std_err = stats.linregress(a1obs, a1ret)
            #*******************************
            # Metrics
            #-------------------------------
            lout = [['','cc','nbias','rmse','slope','intercept']]
            for minmax in lminmax:
                vmin,vmax = minmax
                a1mask = ma.masked_outside(a1obs, vmin*1000., vmax*1000.).mask
                a1retTmp = ma.masked_where(a1mask, a1ret).compressed()
                a1obsTmp = ma.masked_where(a1mask, a1obs).compressed()

                cc = np.ma.corrcoef(a1retTmp,a1obsTmp, allow_masked=True)[0][1] 
                nbias= (a1retTmp - a1obsTmp).mean()/a1obsTmp.mean()
                rmse = np.sqrt(((a1retTmp-a1obsTmp)**2).mean())

                ltmp = ['%d~%d'%(vmin,vmax),cc,nbias,rmse,slope,intercept]
                lout.append(ltmp)
 
            stampOut = '%s.%s.%s.%s'%(exprTmp, region, surftype, season)
            outPath = figDir + '/metrics.stop.%s.%s.csv'%(rettype, stampOut)
            sout = util.list2csv(lout)
            f=open(outPath,'w'); f.write(sout); f.close()
            print ''
            print outPath
            print 'cc=',cc
            print 'bias=',nbias
            print 'rmse=',rmse
            print 'slope=',slope
            print 'intercept=',intercept
