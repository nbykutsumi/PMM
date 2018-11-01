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
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
from gv_fsub import *
import pickle

iYM = [1998,4]
eYM = [2004,10]
#iYM = [2005,4]
#eYM = [2014,9]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()


ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
#ldomain = gv.domains
#ldomain = ['FLORIDA-SFL-N']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

res = 0.02
miss= -9999.
offset_bef = 15    # minutes
offset_aft = 45    # minutes
ldMnt = range(-offset_bef, offset_aft+1)

nlev  = 40

#thdist  = 2.5 # km
#thdist  = 5.0 # km
thdist  = 10.0 # km

ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

def ret_dBBoxes(res):
    #-- read sitelist_summary ---
    listDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
    listPath = listDir + "/sitelist_summary.csv"
    f = open(listPath, "r"); lines=f.readlines(); f.close()
    dnynx       = {}
    dBBox       = {}
    dBBoxBnd    = {}
    for line in lines[1:]:
        line = line.strip().split(",")
        domain = line[3]
        lllat  = float(line[4])
        urlat  = float(line[5])
        lllon  = float(line[6])
        urlon  = float(line[7])
   
        #-- pixel centers -- 
        lat0 = float( Decimal(lllat).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
        lon0 = float( Decimal(lllon).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
        lat1 = float( Decimal(urlat).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )
        lon1 = float( Decimal(urlon).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )

        a1lat_map = arange(lat0,lat1+res*0.1, res)
        a1lon_map = arange(lon0,lon1+res*0.1, res)

        lllat = a1lat_map.min()
        lllon = a1lon_map.min()
        urlat = a1lat_map.max()
        urlon = a1lon_map.max()

        ny = len(a1lat_map)
        nx = len(a1lon_map)
        #-- pixcel edge --
        lllatBnd = lllat - 0.5*res
        lllonBnd = lllon - 0.5*res
        urlatBnd = urlat + 0.5*res
        urlonBnd = urlon + 0.5*res

        dBBox[domain]    = [[lllat, lllon],[urlat,urlon]]
        dBBoxBnd[domain] = [[lllatBnd,lllonBnd],[urlatBnd,urlonBnd]]
        dnynx[domain]    = [ny,nx] 
    
    return dnynx, dBBox, dBBoxBnd



def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)


def load_a2mindist(domain,Year,Mon):
    gvmapDir = gvmapbaseDir + '/%s/%04d%02d'%(domain,Year,Mon)
    mindistPath = gvmapDir + '/map.mindist.%s.%04d%02d.npy'%(domain,Year,Mon)
    return np.load(mindistPath)




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


def ret_gNameLatLon(domain):
    iYM = [2005,4]
    eYM = [2014,9]
    
    lYM = util.ret_lYM(iYM, eYM)
    lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
    lYM = lYM[::-1]
    #-- load gauge data
    daPrcp  = {}
    dLat    = {}
    dLon    = {}
    lgName  = []
    for YM in lYM:
        Year,Mon = YM
        try:
            lgNameTmp  = dgName[domain, Year, Mon]
        except KeyError:
            continue
        for gName in lgNameTmp:
            if gName not in lgName:
                lgName.append(gName)
                dLat[gName], dLon[gName] = gv.dlatlon[domain,gName]

    return lgName, dLat, dLon



#*******************************************************
dnynx, dBBox, dBBoxBnd = ret_dBBoxes(res)


for domain in ldomain:
    BBox    = dBBox[domain]
    BBoxBnd = dBBoxBnd[domain]
    [[lllatBnd, lllonBnd], [urlatBnd,urlonBnd]] = dBBoxBnd[domain]

    lgName, dLat, dLon = ret_gNameLatLon(domain)

    ##-- load elevation data
    #delev   = {}
    #for gName in lgName:
    #    lat = dLat[gName]
    #    lon = dLon[gName]
    #    delev[gName] = ret_elevation(lat,lon)

    for YM in lYM:
        Year, Mon = YM

        #-- satellite overpass
        sateDir = satebaseDir + '/%s/%04d%02d'%(domain, Year, Mon)

        if not os.path.exists(sateDir):
            print 'No directory', sateDir
            continue
        ssearch = sateDir + '/prcp.*.npy'
        lsatePath = glob.glob(ssearch)
        lsatePath = sorted(lsatePath)


        #-- empty container
        aprof  = []
        aeSurf = []
        anSurf = []
        anSurfBin= []
        asatelat = []
        asatelon = []
        arainType= []

        aglat    = []
        aglon    = []
        agName   = []
        adtime   = [] 

        #for satePath in lsatePath:
        for satePath in lsatePath:
            print satePath
            fileName = os.path.basename(satePath)
            gNum    = fileName.split('.')[2]
            ietime  = fileName.split('.')[1]

            satelatPath = sateDir + '/lat.%s.%s.npy'%(ietime, gNum)
            satelonPath = sateDir + '/lon.%s.%s.npy'%(ietime, gNum)
            rangebinPath= sateDir + '/rangeBinNum.%s.%s.npy'%(ietime, gNum)
            satetimePath= sateDir + '/dtime.%s.%s.npy'%(ietime, gNum)

            eSurfPath   = sateDir + '/eSurf.%s.%s.npy'%(ietime, gNum)
            rainTypePath= sateDir + '/rainType.%s.%s.npy'%(ietime, gNum)

            # load sateprcp
            a3sateprcp = np.load(satePath)
            
            a2satelat  = np.load(satelatPath)
            a2satelon  = np.load(satelonPath)
            a1satetime = np.load(satetimePath)
            a3rangebin = np.load(rangebinPath)

            a2eSurf    = np.load(eSurfPath)
            a2rainType = np.load(rainTypePath)

            # gauge loop--- 
            for gName in lgName:

                # find gauge-matching y and x index over a2satelat & a2satelon
                glat, glon = dLat[gName], dLon[gName]
                a1xsate, a1ysate = gv_fsub.gauge_match_pyxy(a2satelon.T, a2satelat.T, glon, glat, thdist)  
                
                a1xsate = ma.masked_less(a1xsate,0).compressed()
                a1ysate = ma.masked_less(a1ysate,0).compressed()


                if len(a1xsate)==0:
                    continue

                #-------------------------------------
                # elevation and groundbin
                #-------------------------------------
                #a1elev = (ones(len(a1xsate))*delev[gName]).astype(float32)

                #-------------------------------------
                # extract from satellite data
                #-------------------------------------
                a2profile= a3sateprcp[a1ysate, a1xsate,:] 
                a1eSurf  = a2eSurf[a1ysate, a1xsate] 
                a1rainType= a2rainType[a1ysate, a1xsate] 

                # bins
                a1groundbin= a3rangebin[a1ysate, a1xsate,2]
                a1nsurfbin = a3rangebin[a1ysate, a1xsate,6]

                # near surface rain
                a1ytmp   = range(a2profile.shape[0])
                a1nSurf  = a2profile[a1ytmp,a1nsurfbin]

                # lat & lon
                a1satelat= a2satelat[a1ysate, a1xsate]
                a1satelon= a2satelon[a1ysate, a1xsate]                     
                a1glat   = ones(len(a1satelat), float32)*glat
                a1glon   = ones(len(a1satelon), float32)*glon

                # gName
                a1gName  = array([gName] * len(a1satelat))

                # DTime
                a1dtime  = a1satetime[a1ysate]

                # profile of nlev range bins above ground
                a2profout= empty([a2profile.shape[0], nlev]).astype(int16)
                for iy in range(a2profile.shape[0]):
                    groundbin     = a1groundbin[iy]
                    if (nlev<=groundbin)&(groundbin <=80):
                        ''' set 20 levels above ground. Note that range bin with ground level is not included. '''
                        a1tmp  = a2profile[iy, groundbin-nlev:groundbin]
                        a2profout[iy] = a1tmp[::-1] # flip order. New order: Low level to High level.

                    else:
                        a2profout[iy] = -8888

                # append to container
                aprof.append(a2profout)
                aeSurf.append(a1eSurf)
                anSurf.append(a1nSurf)
                arainType.append(a1rainType)

                asatelat.append(a1satelat)
                asatelon.append(a1satelon)
                aglat.append(a1glat)
                aglon.append(a1glon)
                adtime.append(a1dtime)

                # corrected nSurfBin 
                anSurfBin.append(a1groundbin-a1nsurfbin-1)
        # make output array
        if len(aprof)==0: continue
 
        aprof   = concatenate(aprof, axis=0).astype(int16)
        aeSurf  = concatenate(aeSurf).astype(float32)
        anSurf  = concatenate(anSurf).astype(int16)
        arainType= concatenate(arainType).astype(int16)
        anSurfBin  = concatenate(anSurfBin).astype(int16)

        asatelat = concatenate(asatelat).astype(float32)
        asatelon = concatenate(asatelon).astype(float32)
        aglat    = concatenate(aglat)
        aglon    = concatenate(aglon)

        adtime   = concatenate(adtime)

        # save satellite obs to file
        obaseDir = '/home/utsumi/mnt/wellshare/GPMGV/GLOC.L2A25'
        outDir   = obaseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)

        util.mk_dir(outDir)
        oprofPath  = outDir  + '/prof.npy'
        oeSurfPath = outDir  + '/eSurf.npy'
        onSurfPath = outDir  + '/nSurf.npy'
        orainTypePath = outDir + '/rainType.npy'
        onSurfBinPath = outDir + '/nSurfBin.npy'
        osatelatPath  = outDir + '/sateLat.npy'
        osatelonPath  = outDir + '/sateLon.npy'
        oglatPath     = outDir + '/gLat.npy'
        oglonPath     = outDir + '/gLon.npy'

        ogNamePath    = outDir + '/gName.npy'
        odtimePath    = outDir + '/dtime.npy'


        #---------------------

        np.save(oprofPath, aprof)
        np.save(oeSurfPath, aeSurf)
        np.save(onSurfPath, anSurf)
        np.save(orainTypePath, arainType)
        np.save(onSurfBinPath, anSurfBin)

        np.save(osatelatPath, asatelat)
        np.save(osatelonPath, asatelon)
        np.save(oglatPath,    aglat)   
        np.save(oglonPath,    aglon) 

        np.save(ogNamePath,   agName) 
        np.save(odtimePath,   adtime)

        print aprof.shape
        print oprofPath


