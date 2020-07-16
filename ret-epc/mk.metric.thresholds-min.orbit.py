import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from mpl_toolkits.basemap import Basemap
#import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import random

myhost = socket.gethostname()
if myhost == 'shui':
    listDir    = '/work/hk02/utsumi/PMM/US/obtlist'
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    srcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk02/PMM/NASA/GPM.GMI/2A/V05'

    mrmsDir  = '/work/hk02/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/2A/V05'
    mrmsDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()

lsurftype = ['ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

calcflag= True
#calcflag= False
lobstype = ['cmb','mrms']
dregion = {'cmb':'GLB', 'mrms':'US'}
#lseason = ['JJADJF']
#lseason = ['JJA']
lseason = ['ALL']
#lseason = [6]
DB_MAXREC = 10000
DB_MINREC = 1000

#lrettype = ['NS','MS','NScmb','MScmb']
#lrettype = ['NS','MS','NScmb','MScmb','GPROF']
lrettype = ['NScmb','GPROF']
#lrettype = ['GPROF'] 
#lrettype = ['NScmb']
lcmbmin = [0, 0.1, 1, 5, 10]
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-wgt-ret-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-wgt-cor-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
nsample = 100  # number of sampled orbits per month (for region=='GLB')

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
        retDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
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


for obstype in lobstype:
    region = dregion[obstype]
    for season in lseason:
        lYM = ret_lym(season)
        for rettype in lrettype:
            if rettype=='GPROF':
                exprTmp = 'gprof'
            else:
                exprTmp = expr
    
            for surftype in lsurftype:
                a1ret = array([])
                a1obs = array([])
    
                for Year,Mon in lYM:
                    if calcflag==False: continue
    
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
    
                    lorbit = sorted(lorbit, key=lambda x:x[3])
                    #lorbit = lorbit[:10]   # test
                    for orbit in lorbit:
                        Day,oid,iy,ey = orbit[2:]
                        print Year,Mon,Day ,oid
                        #-- Read surface type ----
                        ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.20140701-S030429-E043702.001921.V05A.HDF5'
                        ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.*.%06d.V05A.HDF5'%(Year,Mon,Day,oid)
                        lgprofPath = glob.glob(ssearch)
                        if len(lgprofPath)==0:
                            print 'No file for',Year,Mon,Day,oid
                            continue
                        gprofPath= lgprofPath[0]
    
                        with h5py.File(gprofPath,'r') as h:
                            a2surftype = h['S1/surfaceTypeIndex'][:]
    
                        #-- Read retrieved precipitation ----
                        if rettype =='GPROF':
                            with h5py.File(gprofPath,'r') as h:
                                a2ret = h['S1/surfacePrecipitation'][:]
    
                        else:
                            retDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(exprTmp,Year,Mon,Day)
                            a2ret  = np.load(retDir + '/nsurf%s.%06d.y-9999--9999.nrec%d.npy'%(rettype, oid, DB_MAXREC))
                            a2ret = ma.masked_less(a2ret, 0.01).filled(0.0) 
    
    
                        #-- Read MRSM over the orbit ------
                        if obstype=='mrms':
                            ssearch = mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.????????.%06d.?????-?????.npy'%(oid)
                            lmrmsPath= glob.glob(ssearch)
                            if len(lmrmsPath)==0: continue
                            mrmsPath= lmrmsPath[0]
                            a2obs  = np.load(mrmsPath)
                            a2obs  = ma.masked_less(a2obs,0).filled(0) 
    
                        #-- Read DPR-Combined -----------------------------------------------
                        elif obstype=='cmb':
                            dprbaseDir = workbaseDir + '/hk02/PMM/NASA/GPM.DPRGMI/2B/V06'
                            dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                            ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
                            try:
                                dprPath = glob.glob(ssearch)[0]
                            except:
                                print 'No DPR file for oid=',oid
                                continue
                    
                            with h5py.File(dprPath, 'r') as h:
                                a2obsOrg = h['NS/surfPrecipTotRate'][:]
    
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
                            a2obs = ave_9grids_2d(a2obsOrg, a1y, a1x, miss=-9999.).filled(-9999.).reshape(ny,nx)
    
    
                        else:
                            print 'check reftyp',obstype
                            sys.exit()
    
                        #-- Extract regional data ---------------------
                        if iy >0:
                            a2surftype = a2surftype[iy:ey+1]
                            a2ret = a2ret[iy:ey+1]
    
                            if obstype !='mrms':
                                a2obs = a2obs[iy:ey+1]
    
                        if obstype=='cmb':
                            a2surftype = a2surftype[:,83:137+1]
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
                        
                        #-- Screen no-precip cases for both datasets--
                        a2mask1 = ma.masked_less_equal(a2ret, 0).mask
                        a2mask2 = ma.masked_less_equal(a2obs,0).mask
                        a2mask  = a2mask1 * a2mask2
    
                        #-- Screen invalid cases ---
                        a2maskinvalid1 = ma.masked_invalid(a2ret).mask 
                        a2maskinvalid2 = ma.masked_invalid(a2obs).mask 
                        a2maskinvalid3 = ma.masked_less(a2ret,0).mask 
                        a2maskinvalid4 = ma.masked_less(a2obs,0).mask 
    
                        a2maskinvalid = a2maskinvalid1 + a2maskinvalid2 + a2maskinvalid3 + a2maskinvalid4  
    
                        #-- mask ------
                        a2mask  = a2mask + a2masksurf + a2maskinvalid
    
                        
    
    
                        a1retTmp = ma.masked_where(a2mask, a2ret).compressed()
                        a1obsTmp = ma.masked_where(a2mask, a2obs).compressed()
        
                        a1ret = np.concatenate([a1ret, a1retTmp])
                        a1obs = np.concatenate([a1obs, a1obsTmp])
    
                #-- Histograms log-scale -------
                #-- shift 0.01mm/h for log-scale --
    
                if calcflag is True:
                    ddat = {}
                    for cmbmin in lcmbmin:
                        a1flag = ma.masked_greater(a1obs, cmbmin).mask
                        a1obstmp = a1obs[a1flag]
                        a1rettmp = a1ret[a1flag]

                        num  = len(a1obstmp)
                        if num==0:
                            cc  = -9999
                            rmse= -9999
                            nme = -9999
                        else:
                            cc   = np.corrcoef(a1obstmp, a1rettmp)[0,1]
                            rmse = np.sqrt(((a1rettmp - a1obstmp)**2).mean())
        
                            nme  = (a1rettmp - a1obstmp).mean() / a1obstmp.mean()
    
                        ddat[cmbmin,'num']  = num
                        ddat[cmbmin,'cc']   = cc
                        ddat[cmbmin,'rmse'] = rmse
                        ddat[cmbmin,'nme']  = nme
    
                #*******************************
                # Save
                #-------------------------------
                #stamp =  '%s.%s.%s.%s.%04d.%02d'%(exprTmp, region, obstype, surftype, Year, Mon)
                stamp =  '%s.%s.%s.%s.%s'%(exprTmp, region, obstype, surftype, season)
    
                pickleDir  = '/home/utsumi/temp/ret/pickle/metrics'
                util.mk_dir(pickleDir)
    
                picklepath = pickleDir + '/dict.%s.bfile'%(stamp)
                if calcflag == True:
                    with open(picklepath, 'wb') as f:
                        pickle.dump(ddat, f)
           
                        print picklepath 

#**************************************
# Make table
#--------------------------------------
for season in lseason:
    lout = []
    llabel= ['',   '',    '']+['rmse']*len(lcmbmin)+['nme']*len(lcmbmin) + ['cc']*len(lcmbmin)
    lout.append(llabel)

    llabel = ['Ref','Surf','PMW']+lcmbmin*3
    lout.append(llabel)

    print lout
    for obstype in lobstype:
        region = dregion[obstype]
        for surftype in lsurftype:
            for rettype in lrettype:
                if rettype=='GPROF':
                    exprTmp = 'gprof'
                else:
                    exprTmp = expr
        
                #*******************************
                # Load
                #-------------------------------
                pickleDir  = '/home/utsumi/temp/ret/pickle/metrics'
                stamp =  '%s.%s.%s.%s.%s'%(exprTmp, region, obstype, surftype, season)
                picklepath = pickleDir + '/dict.%s.bfile'%(stamp)
    
                with open(picklepath, 'r') as f:
                    ddat = pickle.load(f)


                line = [obstype,surftype,rettype]
                for var in ['rmse','nme','cc']:
                    for cmbmin in lcmbmin:
                        dat = ddat[cmbmin,var]
                        line.append(dat)
                  
 
                lout.append(line)

    #---------------------------
    sout = util.list2csv(lout)
    csvpath = figDir + '/metrics.prec.%s.csv'%(season)
    f=open(csvpath,'w'); f.write(sout); f.close()
    print csvpath
