from numpy import *
import myfunc.util as util
import os, sys
import myfunc.IO.GPM as GPM
from datetime import datetime, timedelta
import numpy as np
import socket
import glob

#prj      = "TRMM.PR"
prdName = "L2A25"
prdVer  = "07"
#var     = "e_SurfRain"
var     = "rain"

#iYM = [1998,4]
#eYM = [2005,3]

iYM = [2005,4]
#eYM = [2014,10]
eYM = [2014,8]
lYM = util.ret_lYM(iYM,eYM)

lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]

gpm = GPM.L2A25(version=prdVer)
gpm2a23 = GPM.L2A23(version=prdVer)

compressed = False

check_obs_periods= True
#check_obs_periods= False

if check_obs_periods==False:
    ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']


hostname = socket.gethostname()
if hostname == 'mizu':
    #listDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
    listDir  = '/home/utsumi/mnt/wellshare/data/GPMGV/sitelist'
    obaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)

elif hostname=='well':
    listDir  = '/media/disk2/share/data/GPMGV/sitelist'
    obaseDir = "/media/disk2/share/GPMGV/%s"%(prdName)

#-- read sitelist_summary ---
listPath = listDir + "/sitelist_summary.csv"
f = open(listPath, "r"); lines=f.readlines(); f.close()

dnwCode  = {}
dBBox    = {}
dyyyymm  = {}
for line in lines[1:]:
    line = line.strip().split(",")
    region = line[0]
    nwName = line[1]
    nwCode = line[2]
    domain = line[3]
    lllat  = float(line[4])
    urlat  = float(line[5])
    lllon  = float(line[6])
    urlon  = float(line[7])
    sYearDat  = int(line[8])
    sMonDat   = int(line[9])
    eYearDat  = int(line[10])
    eMonDat   = int(line[11])
    lyyyymm   = line[12:]

    key         = domain
    dBBox[key]  = [[lllat,lllon],[urlat,urlon]]
    dyyyymm[key]= lyyyymm

#----------------------------
for YM in lYM:
    print YM
    Year,Mon = YM
    lsrcPath = gpm.list_granule(Year,Mon)
    lsrcPath = sorted(lsrcPath)
    for srcPath in lsrcPath:
        print srcPath
        if compressed==False:
            if srcPath[-3:]==".gz": continue
        else:
            if srcPath[-3:]!=".gz": continue

        fileName = os.path.basename(srcPath)
        gNum     = fileName[12:12+5]

        #-- path 2A23 ---
        srcDir2a23  = gpm2a23.rootDir + '/TRMM.PR/L2A23/%s/%04d/%02d'%(prdVer, Year, Mon)
        #srcPath2a23 = srcDir2a23 + '/T1PR2001051619953_2A23F0007.01'
        #srcPath2a23 = srcDir2a23 + '/T1PR%04d%02d%02d%s_2A23F0007.01'%(Year,Mon,Day,gNum)
        searchPath2a23 = srcDir2a23 + '/T1PR%04d%02d??%s_2A23?????.*'%(Year,Mon,gNum)
        #print searchPath2a23 
        srcPath2a23 = glob.glob(searchPath2a23)[0]
        #print srcPath2a23
        #----------------

        #if int(gNum) < 95262: continue

        Lat  = gpm.load_var_granule(srcPath, "Latitude")
        Lon  = gpm.load_var_granule(srcPath, "Longitude")
        dtime= gpm.load_dtime_granule(srcPath)
        Dat  = gpm.load_var_granule(srcPath, var)
        rangeBinNum= gpm.load_var_granule(srcPath, "rangeBinNum")
        eSurf= gpm.load_var_granule(srcPath, "e_SurfRain")
        method= gpm.load_var_granule(srcPath, "method")
        ##nSurf= gpm.load_var_granule(srcPath, "nearSurfRain")
        rainType= gpm.load_var_granule(srcPath, "rainType") 

        stormH  = gpm2a23.load_var_granule(srcPath2a23, "stormH") 

        #print rangeBinNum.shape
        #print rangeBinNum[:,:,3]

        #- check every network
        if check_obs_periods == True:
            lkey = dBBox.keys()
        elif check_obs_periods == False:
            lkey = ldomain
        else:
            print 'check_obs_periods=',check_obs_periods
            sys.exit()

        for key in lkey:
            #if key != ("IOWA","SFB"): continue

            BBox    = dBBox[key]
            lyyyymm = dyyyymm[key]       
            yyyymm  = "%04d%02d"%(YM[0],YM[1])

            if check_obs_periods==True:
                if yyyymm not in lyyyymm: continue

            # filter
            a1mask= GPM.functions.ret_extract_a1mask(Lat=Lat, Lon=Lon,  BBox=BBox)
            a1mask= ~a1mask   # masked =True --> masked=False
            if a1mask.sum()==0: continue
 
            LatTmp= Lat[a1mask]
            LonTmp= Lon[a1mask]
            dtimeTmp = dtime[a1mask]
            DatTmp= Dat[a1mask]
            rangeBinNumTmp = rangeBinNum[a1mask]
            eSurfTmp = eSurf[a1mask]
            methodTmp= method[a1mask]
            #nSurfTmp = nSurf[a1mask]
            rainTypeTmp = rainType[a1mask]
            stormHTmp   = stormH[a1mask]

            if len(dtimeTmp) <2: continue
            sDTime = dtimeTmp[0]
            eDTime = dtimeTmp[-1]

            stime  = "%04d%02d%02d%02d%02d%02d"%(sDTime.year, sDTime.month, sDTime.day, sDTime.hour, sDTime.minute, sDTime.second)    
            etime  = "%04d%02d%02d%02d%02d%02d"%(eDTime.year, eDTime.month, eDTime.day, eDTime.hour, eDTime.minute, eDTime.second)    

 
            # save
            print key
            outDir   = obaseDir + "/%s"%(key)
            util.mk_dir(outDir)
            outDir   = outDir + "/%04d%02d"%(Year,Mon)
            util.mk_dir(outDir)

            latPath  = outDir + "/lat.%s-%s.%s.npy"%(stime,etime,gNum)
            lonPath  = outDir + "/lon.%s-%s.%s.npy"%(stime,etime,gNum)
            dtimePath= outDir + "/dtime.%s-%s.%s.npy"%(stime,etime,gNum)
            prcpPath = outDir + "/prcp.%s-%s.%s.npy"%(stime,etime,gNum)
            rbinPath = outDir + "/rangeBinNum.%s-%s.%s.npy"%(stime,etime,gNum)
            eSurfPath = outDir + "/eSurf.%s-%s.%s.npy"%(stime,etime,gNum)
            methodPath = outDir + "/method.%s-%s.%s.npy"%(stime,etime,gNum)
            #nSurfPath = outDir + "/nSurf.%s-%s.%s.npy"%(stime,etime,gNum)
            rainTypePath = outDir + "/rainType.%s-%s.%s.npy"%(stime,etime,gNum)
            stormHPath   = outDir + "/stormH.%s-%s.%s.npy"%(stime,etime,gNum)

            np.save(latPath,  LatTmp.astype(float32))
            np.save(lonPath,  LonTmp.astype(float32))
            np.save(dtimePath,dtimeTmp)
            np.save(prcpPath, DatTmp.astype(int16))
            np.save(rbinPath, rangeBinNumTmp.astype(int16))
            np.save(eSurfPath,eSurfTmp.astype(float32))
            np.save(methodPath,methodTmp.astype(int16))
            ##np.save(nSurfPath,nSurfTmp.astype(float32))
            np.save(rainTypePath,rainTypeTmp.astype(int16))
            np.save(stormHPath,stormHTmp.astype(int16))

            print '-- out ---'
            print rainTypePath
