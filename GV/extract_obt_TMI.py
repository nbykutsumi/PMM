from numpy import *
import myfunc.util as util
import os, sys
import myfunc.IO.GPM as GPM
from datetime import datetime, timedelta
import numpy as np

sensor  = 'TRMM.TMI'
prdName = '2A-CLIM'
prdVer  = 'V05'
minorVer= 'A'

iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM,eYM)

lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]

gpm = GPM.L2A_GPROF_HDF5(sensor=sensor, prdName=prdName, version=prdVer, minorversion=minorVer)

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
    print YM
    Year,Mon = YM
    lsrcPath = gpm.list_granule(Year,Mon)
    lsrcPath = sorted(lsrcPath)
    for ipath, srcPath in enumerate(lsrcPath):
    #for ipath, srcPath in enumerate(lsrcPath[10:]):
        print srcPath

        fileName = os.path.basename(srcPath)
        gNum     = fileName.split('.')[-3]
        print gNum

        #if int(gNum) < 95262: continue

        Lat  = gpm.load_var_granule(srcPath, 'S1/Latitude')
        Lon  = gpm.load_var_granule(srcPath, 'S1/Longitude')
        dtime= gpm.load_dtime_granule(srcPath)

        #-- load profile data for reconstruction ---
        if ipath ==0:
            clusterProf= gpm.load_var_granule(srcPath, 'GprofDHeadr/clusterProfiles')
            hgtTopLayer= gpm.load_var_granule(srcPath, 'GprofDHeadr/hgtTopLayer')
            species    = gpm.load_var_granule(srcPath, 'GprofDHeadr/speciesDescription')
            species    = [''.join( map(chr, line) ) for line in species] 

            for isp, sp in enumerate(species):
                print isp, sp
                '''
                0 Rain Water Content
                1 Cloud Water Content
                2 Ice Water Content
                3 Snow Water Content
                4 Grauple/Hail Content
                Hydrometeor Unit:[g/m3] see filespec.TRMM.V7.2A12
                '''

            # save
            obaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
            clusterProfPath= obaseDir + "/clusterProf.npy"
            speciesPath    = obaseDir + "/species.npy"
            hgtTopLayerPath= obaseDir + "/hgtTopLayer.npy" 

            np.save(clusterProfPath, clusterProf)
            np.save(speciesPath, species)
            np.save(hgtTopLayerPath, hgtTopLayer)



        #---------------------------------------------

        qFlag    = gpm.load_var_granule(srcPath, 'S1/qualityFlag')
        eSurf    = gpm.load_var_granule(srcPath, 'S1/surfacePrecipitation')
        mlPrecip = gpm.load_var_granule(srcPath, 'S1/mostLikelyPrecipitation')
        tIndex   = gpm.load_var_granule(srcPath, 'S1/profileTemp2mIndex') 
        profNum  = gpm.load_var_granule(srcPath, 'S1/profileNumber') 
        profScale= gpm.load_var_granule(srcPath, 'S1/profileScale') 


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
 
            LatTmp= Lat[a1mask]
            LonTmp= Lon[a1mask]
            dtimeTmp = dtime[a1mask]

            qFlagTmp    = qFlag[a1mask]
            eSurfTmp    = eSurf[a1mask]
            mlPrecipTmp = mlPrecip[a1mask]

            # for reconstruction
            tIndexTmp    = tIndex[a1mask]
            profNumTmp   = profNum[a1mask]
            profScaleTmp = profScale[a1mask]

            '''
            ny,nx = tIndexTmp.shape
            nz    = 28
            a4DatTmp = empty([5,ny,nx,nz]).astype(float32) 
            for ispecies in [0,1,2,3,4]:
                a1profNum   = profNumTmp[:,:,ispecies].flatten() -1
                a1profScale = profScaleTmp[:,:,ispecies].flatten() 
                a2prof = a1profScale.reshape(-1,1) * clusterProf[a1profNum,:,a1tIndex,ispecies]
                print a2prof.shape
                a4DatTmp[ispecies] = a2prof.reshape(ny,nx,-1)

            DatTmp = a4DatTmp.sum(axis=0)
            print DatTmp.shape
            #-- quality flag -- 

            if eSurfTmp.max() > 5:
                print '-----------------'
                imax = np.argmax(eSurfTmp)
                y    = int(imax/208)
                x    = imax - y*208 
                print eSurfTmp[y,x], eSurfTmp.max()
                print DatTmp[y,x,:]
                #print ''
                #print 'tIndex',tIndexTmp[y,x]
                #print 'profNum',profNumTmp[y,x]
                #print 'profScale',profScaleTmp[y,x]
                T = tIndexTmp[y,x] -1

                a1prof = zeros(28)
                for S in [0,1,2,3,4]: 
                    scale = profScaleTmp[y,x,S]
                    P     = profNumTmp[y,x,S] -1
                    a1tmp = scale*clusterProf[P,:,T,S] 
                    print P,T,S,scale
                    print S, a1tmp
                    a1prof = a1prof + a1tmp

                print 'all',a1prof  # [g/m3]
                a1prof_rate = a1prof * 4.17*60*60/1000.  # g/m3 --> mm/h, falling velocity 4.17m/h is assumed.
                print 'mm/h',a1prof_rate
                print '-----------------'
            '''

            #------------------

            if len(dtimeTmp) <2: continue
            sDTime = dtimeTmp[0]
            eDTime = dtimeTmp[-1]

            stime  = "%04d%02d%02d%02d%02d%02d"%(sDTime.year, sDTime.month, sDTime.day, sDTime.hour, sDTime.minute, sDTime.second)    
            etime  = "%04d%02d%02d%02d%02d%02d"%(eDTime.year, eDTime.month, eDTime.day, eDTime.hour, eDTime.minute, eDTime.second)    

 
            # save
            #obaseDir = "/work/a01/utsumi/GPMGV/%s"%(prdName)
            obaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
            outDir   = obaseDir + "/%s"%(key)
            util.mk_dir(outDir)
            outDir   = outDir + "/%04d%02d"%(Year,Mon)
            util.mk_dir(outDir)


            latPath  = outDir + "/lat.%s-%s.%s.npy"%(stime,etime,gNum)
            lonPath  = outDir + "/lon.%s-%s.%s.npy"%(stime,etime,gNum)
            dtimePath= outDir + "/dtime.%s-%s.%s.npy"%(stime,etime,gNum)

            qFlagPath = outDir + "/qFlag.%s-%s.%s.npy"%(stime,etime,gNum)
            eSurfPath = outDir + "/eSurf.%s-%s.%s.npy"%(stime,etime,gNum)
            mlPrecipPath = outDir + "/mlPrecip.%s-%s.%s.npy"%(stime,etime,gNum)
            tIndexPath   = outDir + "/tIndex.%s-%s.%s.npy"%(stime,etime,gNum)
            profNumPath  = outDir + "/profNum.%s-%s.%s.npy"%(stime,etime,gNum)
            profScalePath= outDir + "/profScale.%s-%s.%s.npy"%(stime,etime,gNum)


            np.save(latPath,  LatTmp.astype(float32))
            np.save(lonPath,  LonTmp.astype(float32))
            np.save(dtimePath, dtimeTmp)
            np.save(qFlagPath, qFlagTmp)
            np.save(eSurfPath, eSurfTmp)
            np.save(mlPrecipPath,  mlPrecipTmp)
            np.save(tIndexPath,    tIndexTmp)
            np.save(profNumPath,   profNumTmp)
            np.save(profScalePath, profScaleTmp)


            print eSurfPath
