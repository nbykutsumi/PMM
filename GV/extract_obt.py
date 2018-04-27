from numpy import *
import myfunc.util as util
import os, sys
import myfunc.IO.GPM as GPM
from datetime import datetime, timedelta
import numpy as np

#prj      = "TRMM.PR"
prdName = "L2A25"
prdVer  = "07"
#var     = "e_SurfRain"
var     = "rain"

iYM = [2014,1]
eYM = [2014,9]
lYM = util.ret_lYM(iYM,eYM)

gpm = GPM.L2A25(version=prdVer)

compressed = False
#-- read sitelist_summary ---
listDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
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
    Year,Mon = YM
    lsrcPath = gpm.list_granule(Year,Mon)
    lsrcPath = sorted(lsrcPath)
    for srcPath in lsrcPath:
        if compressed==False:
            if srcPath[-3:]==".gz": continue
        else:
            if srcPath[-3:]!=".gz": continue

        fileName = os.path.basename(srcPath)
        gNum     = fileName[12:12+5]

        Dat  = gpm.load_var_granule(srcPath, var)
        Lat  = gpm.load_var_granule(srcPath, "Latitude")
        Lon  = gpm.load_var_granule(srcPath, "Longitude")
        dtime= gpm.load_dtime_granule(srcPath)
        rangeBinNum= gpm.load_var_granule(srcPath, "rangeBinNum")

 

        #print rangeBinNum.shape
        #print rangeBinNum[:,:,3]

        #- check every network
        for key in dBBox.keys():
            #if key != ("IOWA","SFB"): continue

            BBox    = dBBox[key]
            lyyyymm = dyyyymm[key]       
            yyyymm  = "%04d%02d"%(YM[0],YM[1])

            if yyyymm not in lyyyymm: continue

            # filter
            a1mask= GPM.functions.ret_extract_a1mask(Lat=Lat, Lon=Lon,  BBox=BBox)
            a1mask= ~a1mask   # masked =True --> masked=False
            if a1mask.sum()==0: continue
 
            DatTmp= Dat[a1mask]
            LatTmp= Lat[a1mask]
            LonTmp= Lon[a1mask]
            dtimeTmp = dtime[a1mask]
            rangeBinNumTmp = rangeBinNum[a1mask]

            if len(dtimeTmp) <2: continue
            sDTime = dtimeTmp[0]
            eDTime = dtimeTmp[1]

            stime  = "%04d%02d%02d%02d%02d%02d"%(sDTime.year, sDTime.month, sDTime.day, sDTime.hour, sDTime.minute, sDTime.second)    
            etime  = "%04d%02d%02d%02d%02d%02d"%(eDTime.year, eDTime.month, eDTime.day, eDTime.hour, eDTime.minute, eDTime.second)    

 
            # save
            print key, DatTmp.shape, DatTmp.max()
            #obaseDir = "/work/a01/utsumi/GPMGV/%s"%(prdName)
            obaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
            outDir   = obaseDir + "/%s"%(key)
            util.mk_dir(outDir)
            outDir   = outDir + "/%04d%02d"%(Year,Mon)
            util.mk_dir(outDir)

            latPath  = outDir + "/lat.%s-%s.%s.npy"%(stime,etime,gNum)
            lonPath  = outDir + "/lon.%s-%s.%s.npy"%(stime,etime,gNum)
            dtimePath= outDir + "/dtime.%s-%s.%s.npy"%(stime,etime,gNum)
            prcpPath = outDir + "/prcp.%s-%s.%s.npy"%(stime,etime,gNum)
            rbinPath = outDir + "/rangeBinNum.%s-%s.%s.npy"%(stime,etime,gNum)
            np.save(latPath,  LatTmp.astype(float32))
            np.save(lonPath,  LonTmp.astype(float32))
            np.save(prcpPath, DatTmp.astype(int16))
            np.save(rbinPath, rangeBinNumTmp.astype(int16))

            print prcpPath