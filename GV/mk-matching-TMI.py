import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
from datetime import datetime, timedelta
import numpy as np
import myfunc.util as util
import GPMGV
import os, sys
import glob
import pickle
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
from gv_fsub import *
from bisect import bisect_left


iYM = [2005,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM
prdName = '2A-CLIM'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()


ldomain = gv.domains
#ldomain = ['FLORIDA-STJ']
#ldomain = ['FLORIDA-SFL-N']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']


miss= -9999.
offset_bef = 15    # minutes
offset_aft = 45    # minutes
ldMnt = range(-offset_bef, offset_aft+1)

nlev  = 19 # at highest 10 km
#thdist  = 2.5 # km
#thdist  = 7.5 # km
thdist  = 15 # km
ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/%s'%(prdName)
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)


def load_gauge(domain, gName, Year,Mon):

    if len(domain.split('-'))==2:
        region, nwName= domain.split('-')
    else:
        region, nwName, tmp = domain.split('-')

    nwCode   = gv.dnwCode[domain]
    gaugebaseDir= '/work/a01/utsumi/data/GPMGV/2A56'
    gaugeDir = gaugebaseDir + '/%s/%s/%04d'%(region,nwName,Year)

    #-- asc type --
    ssearch  = gaugeDir + '/2A56_%s_%s_%04d%02d_?.asc'%(nwCode, gName, Year, Mon)
    lgaugePath= glob.glob(ssearch)

    #-- 2a56 type --
    if len(lgaugePath)==0:
        #ssearch  = gaugeDir + '/%s-%s-%04d%02d.2a56'%(nwName.split('_')[0], gName, Year, Mon)
        ssearch  = gaugeDir + '/%s-%s-%04d%02d.2a56'%(nwCode, gName, Year, Mon)
        lgaugePath= glob.glob(ssearch)

        #if len(lgaugePath) == 0:
        #    ssearch  = gaugeDir + '/%s-%s-%04d%02d.2a56'%(nwName, gName, Year, Mon)
        #    lgaugePath= glob.glob(ssearch)

    if len(lgaugePath) !=0:
        gaugePath = lgaugePath[0]
        

    else:
        print 'no file',domain,gName,YM
        print 'search=',ssearch
        print 'in',gaugeDir
        sys.exit()

    # load
    sfx = gaugePath.split('.')[-1]
    if sfx =='asc':
        head, aDTime, aPrcp = gv.load_obsfile_asc(gaugePath)
    elif sfx =='2a56':
        head, aDTime, aPrcp = gv.load_obsfile_2a56(gaugePath)
    else:
        print 'Check gauge file',gaugePath
        sys.exit()

    return aPrcp




def load_gtopo(lat,lon):
    orogDir  = "/data1/hjkim/GTOPO30"

    ullat = int( (lat - (-60))/50. )*50. + 50 -60.
    ullon = int( (lon - (-180))/40.)*40. -180.

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<ullon)&(ullon<0):
        WE = "W"
    else:
        WE = "E"

    #orogPath = orogDir + "/E060N40.DEM"
    orogPath = orogDir + "/%s%03d%s%02d.DEM"%(WE, abs(ullon), SN, abs(ullat))
    dmwPath  = orogDir + "/%s%03d%s%02d.DMW"%(WE, abs(ullon), SN, abs(ullat))
    print dmwPath
    """
    BYTEORDER      M
    LAYOUT       BIL
    NROWS         6000
    NCOLS         4800
    NBANDS        1
    NBITS         16
    BANDROWBYTES         9600
    TOTALROWBYTES        9600
    BANDGAPBYTES         0
    NODATA        -9999
    ULXMAP        60.00416666666667
    ULYMAP        39.99583333333333
    XDIM          0.00833333333333
    YDIM          0.00833333333333
    """

    ny, nx = 6000, 4800
    a = flipud(fromfile(orogPath, "int16").byteswap().reshape(ny,nx))

    # load DMW file
    f=open(dmwPath, "r"); lines=f.readlines(); f.close()
    lonmin = float(lines[4])
    latmax = float(lines[5])
    lonmax = lonmin + 0.00833333333333*(nx-1)
    latmin = latmax - 0.00833333333333*(ny-1)

    dlat = 0.00833333333333
    dlon = 0.00833333333333
    para = {}
    para["ny"] = ny
    para["nx"] = nx
    para["lllat"]=latmin
    para["lllon"]=lonmin
    para["urlat"]=latmax
    para["urlon"]=lonmax
    para["dlat" ]=dlat
    para["dlon" ]=dlon

    para["Lat"]  = arange(latmin,latmax+0.5*dlat, dlat)
    para["Lon"]  = arange(lonmin,lonmax+0.5*dlon, dlon)

    return a,para

def ret_elevation(lat,lon):
    a2elev, para = load_gtopo(lat, lon)
    lllat= para['lllat']
    lllon= para['lllon']
    dlat = para['dlat']
    dlon = para['dlon']
    y    = int((lat - lllat)/dlat)
    x    = int((lon - lllon)/dlon)
    return a2elev[y,x]


#*******************************************************

for YM in lYM:
    Year, Mon = YM
    for domain in ldomain:
        if (domain,Year,Mon) not in dgName.keys():
            continue

        #-- satellite overpasses
        sateDir = satebaseDir + '/%s/%04d%02d'%(domain, Year, Mon)


        if not os.path.exists(sateDir):
            print 'No directory', sateDir
            continue
        ssearch = sateDir + '/profNum.*.npy'
        lsatePath = glob.glob(ssearch)
        lsatePath = sorted(lsatePath)

        #-- load gauge data
        lgName  = dgName[domain, Year, Mon]
        daPrcp  = {}
        dLat    = {}
        dLon    = {}
        for gName in lgName:
            dLat[gName], dLon[gName] = gv.dlatlon[domain,gName]

            daPrcp[gName] = load_gauge(domain, gName, Year, Mon)


        #-- load elevation data
        delev   = {}
        for gName in lgName:
            lat = dLat[gName]
            lon = dLon[gName]
            delev[gName] = ret_elevation(lat,lon)  

 
        #-- set ground bin
        hgtTopLayerPath= satebaseDir + '/hgtTopLayer.npy'
        a1hgtTop = np.load(hgtTopLayerPath)*1000. #[m] 

        print a1hgtTop
        dgroundBin={}
        for gName in lgName:
            elev = delev[gName]
            dgroundBin[gName] = bisect_left(a1hgtTop, elev) 


        #-- empty container
        aqFlag     = []
        aprofNum   = []
        atIndex    = []
        aprofScale = []
        aeSurf     = []  # surface precip by GPROF
        amlPrcp    = []  # most likely surface precip by GPROF
        asatelat = []
        asatelon = []
        aglat    = []
        aglon    = []
        agName   = []
        adtime   = [] 
        agroundBin= []

        agv    = []


        #-- empty container for pool
        pqFlag     = []
        pprofNum   = []
        ptIndex    = []
        pprofScale = []
        peSurf     = []  # surface precip by GPROF
        pmlPrcp    = []  # most likely surface precip by GPROF
        psatelat = []
        psatelon = []
        pglat    = []
        pglon    = []
        pgName   = []
        pdtime   = [] 
        pgroundBin= []
        pgelev    = []

        pgv    = []
        pngv   = []  # only for pool


        #for satePath in lsatePath:
        for satePath in lsatePath:
            print satePath
            fileName = os.path.basename(satePath)
            gNum    = fileName.split('.')[2]
            ietime  = fileName.split('.')[1]

            satelatPath = sateDir + '/lat.%s.%s.npy'%(ietime, gNum)
            satelonPath = sateDir + '/lon.%s.%s.npy'%(ietime, gNum)
            satetimePath= sateDir + '/dtime.%s.%s.npy'%(ietime, gNum)
            qFlagPath   = sateDir + '/qFlag.%s.%s.npy'%(ietime, gNum)
            eSurfPath   = sateDir + '/eSurf.%s.%s.npy'%(ietime, gNum)
            mlPrcpPath  = sateDir + '/mlPrecip.%s.%s.npy'%(ietime, gNum)
            tIndexPath  = sateDir + '/tIndex.%s.%s.npy'%(ietime, gNum)
            profNumPath = sateDir + '/profNum.%s.%s.npy'%(ietime, gNum)
            profScalePath = sateDir + '/profScale.%s.%s.npy'%(ietime, gNum)

            # load sateprcp
            a2satelat  = np.load(satelatPath)
            a2satelon  = np.load(satelonPath)
            a1satetime = np.load(satetimePath)

            a2qFlag    = np.load(qFlagPath)
            a2eSurf    = np.load(eSurfPath)
            a2mlPrcp   = np.load(mlPrcpPath)
            a2tIndex   = np.load(tIndexPath)
            a3profNum  = np.load(profNumPath)
            a3profScale= np.load(profScalePath)

            #--------------------
            dtimeA = dtime_round_Mnt(a1satetime.min())
            dtimeB = dtime_round_Mnt(a1satetime.max()) 
            ldtime0 = util.ret_lDTime(dtimeA, dtimeB, ddtime_1min)

            #-- time loop in each granule
            for dtime0 in ldtime0:
                if dtime0.month != Mon:
                    continue
                # mask satellite data
                mskTmp = ma.masked_outside(a1satetime, dtime0, dtime0+ddtime_1min - timedelta(microseconds=1)).mask
                a1timeidx  = ma.masked_where(mskTmp, arange(len(a1satetime)).astype(int32)).compressed()


                a2satelatTmp  = a2satelat[a1timeidx,:]
                a2satelonTmp  = a2satelon[a1timeidx,:]
                a2qFlagTmp    = a2qFlag[a1timeidx,:]
                a2eSurfTmp    = a2eSurf[a1timeidx,:]
                a2mlPrcpTmp   = a2mlPrcp[a1timeidx,:]
                a2tIndexTmp   = a2tIndex[a1timeidx,:]
                a3profNumTmp  = a3profNum[a1timeidx,:]

                #------------------
                idtime_offset =  dtime0 - offset_bef * ddtime_1min + ddtime_1min
                edtime_offset =  dtime0 + offset_aft * ddtime_1min + ddtime_1min

                imin    = int( (idtime_offset - datetime(Year,Mon,1,0)).total_seconds()/60.)
                emin    = int( (edtime_offset - datetime(Year,Mon,1,0)).total_seconds()/60.)

                if idtime_offset.month != Mon: continue
                if edtime_offset.month != Mon: continue



                # Dictionary for pooling 
                lxy     = []
                dngv    = {}
                dla1gv  = {} 
                dlgName = {} 
                dlglat  = {} 
                dlglon  = {} 
                dlgelev = {} 


                # gauge loop--- 

                lgName  = dgName[domain, Year, Mon]
                for gName in lgName:

                    # find gauge-matching y and x index over a2satelat & a2satelon
                    glat, glon = dLat[gName], dLon[gName]
                    a1xsate, a1ysate = gv_fsub.gauge_match_pyxy(a2satelonTmp.T, a2satelatTmp.T, glon, glat, thdist)  
                    
                    a1xsate = ma.masked_less(a1xsate,0).compressed()
                    a1ysate = ma.masked_less(a1ysate,0).compressed()


                    if len(a1xsate)==0:
                        continue

                    #-------------------------------------
                    # elevation and groundbin
                    #-------------------------------------
                    a1elev = (ones(len(a1xsate))*delev[gName]).astype(float32)
                    a1groundbin=(ones(len(a1xsate))*dgroundBin[gName]).astype(float32)

                    #-------------------------------------
                    # gauge data
                    #-------------------------------------
                    a2gvprcp  = empty([len(a1xsate), len(ldMnt)]).astype(float32)


                    a1gvTmp   = daPrcp[gName][imin:emin+1]

                    for itmp in range(a2gvprcp.shape[0]):

                        a2gvprcp[itmp] = a1gvTmp

                    #-------------------------------------
                    # extract from satellite data
                    #-------------------------------------
                    a1qFlag  = a2qFlagTmp[a1ysate, a1xsate] 
                    a1eSurf  = a2eSurfTmp[a1ysate, a1xsate] 
                    a1mlPrcp = a2mlPrcpTmp[a1ysate, a1xsate] 
                    a1tIndex = a2tIndexTmp[a1ysate, a1xsate] 
                    a2profNum= a3profNum[a1ysate,a1xsate,:]
                    a2profScale= a3profScale[a1ysate,a1xsate,:]

                    # lat & lon
                    a1satelat= a2satelatTmp[a1ysate, a1xsate]
                    a1satelon= a2satelonTmp[a1ysate, a1xsate]                     
                    a1glat   = ones(len(a1satelat), float32)*glat
                    a1glon   = ones(len(a1satelon), float32)*glon

                    # gName
                    a1gName  = array([gName] * len(a1satelat))

                    # DTime
                    a1dtime  = array([dtime0]* len(a1satelat))

                    # append to container
                    agv.append(a2gvprcp)

                    aqFlag.append(a1qFlag)
                    aprofNum.append(a2profNum)
                    aprofScale.append(a2profScale)
                    atIndex.append(a1tIndex)

                    aeSurf.append(a1eSurf)
                    amlPrcp.append(a1mlPrcp)

                    asatelat.append(a1satelat)
                    asatelon.append(a1satelon)
                    aglat.append(a1glat)
                    aglon.append(a1glon)
                    agName.append(a1gName)
                    adtime.append(a1dtime)

                    # groundBin 
                    agroundBin.append(a1groundbin)

                    #**********************
                    # pool gauge obs in the same
                    # satellite footprint (xy)
                    #**********************
                    for xy in zip(a1xsate, a1ysate):
                        if xy not in lxy:
                            lxy.append(xy)
                            dngv[xy]   = 1
                            dla1gv[xy] = [a1gvTmp]
                            dlgName[xy]= [gName]
                            dlglat[xy] = [glat]
                            dlglon[xy] = [glon]
                            dlgelev[xy]= [delev[gName]]

                        else:
                            dngv[xy] = dngv[xy]+1
                            dla1gv[xy].append(a1gvTmp)
                            dlgName[xy].append(gName)
                            dlglat[xy].append(glat)
                            dlglon[xy].append(glon)
                            dlgelev[xy].append(delev[gName])
                #-- shrink pooled data --

                if len(lxy)==0: continue

                for xy in lxy:
                    if dngv[xy]==1:
                        a1gvTmp = array(dla1gv[xy][0])
                        elev =dlgelev[xy][0]

                    else:
                        a2gvTmp = array(dla1gv[xy])

                        a2gvAbsTmp= abs(a2gvTmp)
                        a1gvTmp = a2gvAbsTmp.mean(axis=0)
                         
                        a1nnega = ma.masked_less(a2gvTmp,0).mask.sum(axis=0)
                        a1gvTmp = (ma.masked_where(a1nnega==0, a1gvTmp)*(-1)).data

                        elev = mean(dlgelev[xy])

                    # append pooled pauge info to be pickled
                    pgv.append(a1gvTmp)
                    groundBin = bisect_left(a1hgtTop, elev)
                    pgroundBin.append(groundBin)
                    pngv.append(dngv[xy])

                    pglat .append(dlglat[xy])
                    pglon .append(dlglon[xy])
                    pgName.append(dlgName[xy])
                    pgelev.append(dlgelev[xy])



 

                # extract satellite for pool
                a1xsate, a1ysate = zip(*lxy)

                #-------------------------------------
                # extract from satellite data for pool
                #-------------------------------------
                p1qFlag  = a2qFlagTmp[a1ysate, a1xsate] 
                p1eSurf  = a2eSurfTmp[a1ysate, a1xsate] 
                p1mlPrcp = a2mlPrcpTmp[a1ysate, a1xsate] 
                p1tIndex = a2tIndexTmp[a1ysate, a1xsate] 
                p2profNum= a3profNum[a1ysate,a1xsate,:]
                p2profScale= a3profScale[a1ysate,a1xsate,:]

                # lat & lon for pool
                p1satelat= a2satelatTmp[a1ysate, a1xsate]
                p1satelon= a2satelonTmp[a1ysate, a1xsate]                     

                # DTime for pool
                p1dtime  = array([dtime0]* len(p1satelat))

                # append to container for pool

                pqFlag.append(p1qFlag)
                pprofNum.append(p2profNum)
                pprofScale.append(p2profScale)
                ptIndex.append(p1tIndex)

                peSurf.append(p1eSurf)
                pmlPrcp.append(p1mlPrcp)

                psatelat.append(p1satelat)
                psatelon.append(p1satelon)
                pdtime.append(p1dtime)



        # make output array
        agv        = concatenate(agv, axis=0).astype(float32)
        
        aqFlag     = concatenate(aqFlag,axis=0)
        aprofNum   = concatenate(aprofNum,axis=0)
        aprofScale = concatenate(aprofScale,axis=0)
        atIndex    = concatenate(atIndex, axis=0)

        aeSurf     = concatenate(aeSurf)
        amlPrcpf   = concatenate(amlPrcp)
        agroundBin = concatenate(agroundBin)

        asatelat = concatenate(asatelat).astype(float32)
        asatelon = concatenate(asatelon).astype(float32)
        aglat    = concatenate(aglat)
        aglon    = concatenate(aglon)

        agName   = concatenate(agName)
        adtime   = concatenate(adtime)


        # make output array for pool
        pgv        = array(pgv).astype(float32)
 
        pqFlag     = concatenate(pqFlag,axis=0)
        pprofNum   = concatenate(pprofNum,axis=0)
        pprofScale = concatenate(pprofScale,axis=0)
        ptIndex    = concatenate(ptIndex, axis=0)

        peSurf     = concatenate(peSurf)
        pmlPrcpf   = concatenate(pmlPrcp)
        #pgroundBin = concatenate(pgroundBin)  # no need to concatenate

        psatelat = concatenate(psatelat).astype(float32)
        psatelon = concatenate(psatelon).astype(float32)

        pdtime   = concatenate(pdtime)
        #sys.exit()       

        # save satellite obs to file
        obaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.%s'%(prdName)
        outDir   = obaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)
        util.mk_dir(outDir)

        ogvPath        = outDir + '/gvprcp.npy'
        oqFlagPath     = outDir + '/qFlag.npy'
        oprofNumPath   = outDir + '/profNum.npy'
        oprofScalePath = outDir + '/profScale.npy'
        otIndexPath    = outDir + '/tIndex.npy'

        oeSurfPath    = outDir + '/eSurf.npy'
        omlPrcpPath   = outDir + '/mlPrecip.npy'
        ogroundBinPath= outDir + '/groundBin.npy'
        osatelatPath  = outDir + '/sateLat.npy'
        osatelonPath  = outDir + '/sateLon.npy'
        oglatPath     = outDir + '/gLat.npy'
        oglonPath     = outDir + '/gLon.npy'

        ogNamePath    = outDir + '/gName.npy'
        odtimePath    = outDir + '/dtime.npy'


        #---------------------

        np.save(ogvPath,    agv)
        np.save(oqFlagPath,    aqFlag)
        np.save(oprofNumPath,    aprofNum)
        np.save(oprofScalePath,    aprofScale)
        np.save(otIndexPath,    atIndex)
        np.save(oeSurfPath, aeSurf)
        np.save(omlPrcpPath, amlPrcp)
        np.save(ogroundBinPath, agroundBin)


        np.save(osatelatPath, asatelat)
        np.save(osatelonPath, asatelon)
        np.save(oglatPath,    aglat)   
        np.save(oglonPath,    aglon) 

        np.save(ogNamePath,   agName) 
        np.save(odtimePath,   adtime)

        print oprofNumPath


        # save satellite obs to file for pool
        pgvPath        = outDir + '/p_gvprcp.npy'
        pqFlagPath     = outDir + '/p_qFlag.npy'
        pprofNumPath   = outDir + '/p_profNum.npy'
        pprofScalePath = outDir + '/p_profScale.npy'
        ptIndexPath    = outDir + '/p_tIndex.npy'

        peSurfPath    = outDir + '/p_eSurf.npy'
        pmlPrcpPath   = outDir + '/p_mlPrecip.npy'
        pgroundBinPath= outDir + '/p_groundBin.npy'
        psatelatPath  = outDir + '/p_sateLat.npy'
        psatelonPath  = outDir + '/p_sateLon.npy'
        pdtimePath    = outDir + '/p_dtime.npy'
        pngvPath      = outDir + '/p_ngv.npy'

        pglatPath     = outDir + '/p_gLat.pickle'
        pglonPath     = outDir + '/p_gLon.pickle'
        pgNamePath    = outDir + '/p_gName.pickle'

        #---------------------

        np.save(pgvPath,    pgv)
        np.save(pqFlagPath,    pqFlag)
        np.save(pprofNumPath,    pprofNum)
        np.save(pprofScalePath,    pprofScale)
        np.save(ptIndexPath,    ptIndex)
        np.save(peSurfPath, peSurf)
        np.save(pmlPrcpPath, pmlPrcp)
        np.save(pgroundBinPath, pgroundBin)


        np.save(psatelatPath, psatelat)
        np.save(psatelonPath, psatelon)
        np.save(pdtimePath,   pdtime)
        np.save(pngvPath,   pngv)

        # pickling
        with open(pglatPath, 'wb') as f:
            pickle.dump(pglat, f) 

        with open(pglonPath, 'wb') as f:
            pickle.dump(pglon, f) 

        with open(pgNamePath, 'wb') as f:
            pickle.dump(pgName, f) 



