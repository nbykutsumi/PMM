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


iYM = [2013,4]
eYM = [2014,9]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()


ldomain = gv.domains
#ldomain = ['FLORIDA-STJ']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

res = 0.02
miss= -9999.
offset_bef = 15    # minutes
offset_aft = 45    # minutes
ldMnt = range(-offset_bef, offset_aft+1)

nlev  = 40

thdist  = 2.5 # km
lgvtype = ['All','Gd']

ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'
def ret_aidx(Lat, Lon, lllatBnd, lllonBnd, nymap, nxmap, a2radmask=None):
    Xmap = ( (Lon - lllonBnd)/res ).astype(int32)
    Ymap = ( (Lat - lllatBnd)/res ).astype(int32)

    msk1 = ma.masked_outside(Xmap, 0, nxmap-1).mask
    msk2 = ma.masked_outside(Ymap, 0, nymap-1).mask
    msk  = ~(msk1 + msk2)

    nysate, nxsate = Lat.shape
    Xsate, Ysate = meshgrid(range(nxsate), range(nysate))

    a1xmap = Xmap[msk] # automatically flatten
    a1ymap = Ymap[msk]
    a1xsate= Xsate[msk]
    a1ysate= Ysate[msk]

    #- mask with radmask
    if type(a2radmask) ==np.ndarray:
        a1idx_sc = arange(len(a1xmap))
        a1idx_sc = ma.masked_where(a2radmask[a1ymap,a1xmap], a1idx_sc).compressed()
        a1xmap    = a1xmap[a1idx_sc]
        a1ymap    = a1ymap[a1idx_sc]
        a1xsate   = a1xsate[a1idx_sc]
        a1ysate   = a1ysate[a1idx_sc]
 
    return a1xmap, a1ymap, a1xsate, a1ysate    

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



#*******************************************************
dnynx, dBBox, dBBoxBnd = ret_dBBoxes(res)

for domain in ldomain:
    nymap,nxmap = dnynx[domain]
    BBox    = dBBox[domain]
    BBoxBnd = dBBoxBnd[domain]
    [[lllatBnd, lllonBnd], [urlatBnd,urlonBnd]] = dBBoxBnd[domain]


    for YM in lYM:
        Year, Mon = YM

        if (domain,Year,Mon) not in dgName.keys():
            continue

        #-- load gauge data
        lgName  = dgName[domain, Year, Mon]
        daPrcp  = {}
        dLat    = {}
        dLon    = {}
        for gName in lgName:
            dLat[gName], dLon[gName] = gv.dlatlon[domain,gName]

            daPrcp[gName] = load_gauge(domain, gName, Year, Mon)
        

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
        aglat    = []
        aglon    = []
        agName   = []
        adtime   = [] 


        agv    = []
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



            # load sateprcp
            a3sateprcp = np.load(satePath)
            
            a2satelat  = np.load(satelatPath)
            a2satelon  = np.load(satelonPath)
            a1satetime = np.load(satetimePath)
            a3rangebin = np.load(rangebinPath)

            a2eSurf    = np.load(eSurfPath)

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


                a3sateprcpTmp = a3sateprcp[a1timeidx,:,:]
                a2satelatTmp  = a2satelat[a1timeidx,:]
                a2satelonTmp  = a2satelon[a1timeidx,:]
                a3rangebinTmp = a3rangebin[a1timeidx,:,:]

                a2eSurfTmp    = a2eSurf[a1timeidx,:]

                ## test -----------------
                #nttmp,nxtmp,nhtmp = a3sateprcpTmp.shape 
                #a1xtmp = range(nxtmp)
                #a1htmp = range(nhtmp)
                #H,X  = meshgrid(a1htmp, a1xtmp)

                #H    = H[:,::-1]    # flip

                #a3sateprcpTmp[:] = H



                idtime_offset =  dtime0 - offset_bef * ddtime_1min + ddtime_1min
                edtime_offset =  dtime0 + offset_aft * ddtime_1min + ddtime_1min

                imin    = int( (idtime_offset - datetime(Year,Mon,1,0)).total_seconds()/60.)
                emin    = int( (edtime_offset - datetime(Year,Mon,1,0)).total_seconds()/60.)

                if idtime_offset.month != Mon: continue
                if edtime_offset.month != Mon: continue

 
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
                    # gauge data
                    #-------------------------------------
                    a2gvprcp  = empty([len(a1xsate), len(ldMnt)]).astype(float32)


                    a1gvTmp   = daPrcp[gName][imin:emin+1]

                    for itmp in range(a2gvprcp.shape[0]):

                        a2gvprcp[itmp] = a1gvTmp

                    #-------------------------------------
                    # extract from satellite data
                    #-------------------------------------
                    a2profile= a3sateprcpTmp[a1ysate, a1xsate,:] 
                    a1eSurf  = a2eSurfTmp[a1ysate, a1xsate] 
                    # bins
                    a1groundbin= a3rangebin[a1ysate, a1xsate,2]
                    a1nsurfbin = a3rangebin[a1ysate, a1xsate,6]

                    # near surface rain
                    a1ytmp   = range(a2profile.shape[0])
                    a1nSurf  = a2profile[a1ytmp,a1nsurfbin]

                    # lat & lon
                    a1satelat= a2satelatTmp[a1ysate, a1xsate]
                    a1satelon= a2satelonTmp[a1ysate, a1xsate]                     
                    a1glat   = ones(len(a1satelat), float32)*glat
                    a1glon   = ones(len(a1satelon), float32)*glon

                    # gName
                    a1gName  = array([gName] * len(a1satelat))

                    # DTime
                    a1dtime  = array([dtime0]* len(a1satelat))

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
                    agv.append(a2gvprcp)
                    aprof.append(a2profout)
                    aeSurf.append(a1eSurf)
                    anSurf.append(a1nSurf)

                    asatelat.append(a1satelat)
                    asatelon.append(a1satelon)
                    aglat.append(a1glat)
                    aglon.append(a1glon)
                    agName.append(a1gName)
                    adtime.append(a1dtime)

                    # corrected nSurfBin 
                    anSurfBin.append(a1groundbin-a1nsurfbin-1)



        # make output array
        agv     = concatenate(agv, axis=0).astype(float32)
        aprof   = concatenate(aprof, axis=0).astype(int16)
        aeSurf  = concatenate(aeSurf).astype(float32)
        anSurf  = concatenate(anSurf).astype(int16)
        anSurfBin  = concatenate(anSurfBin).astype(int16)

        asatelat = concatenate(asatelat).astype(float32)
        asatelon = concatenate(asatelon).astype(float32)
        aglat    = concatenate(aglat)
        aglon    = concatenate(aglon)

        agName   = concatenate(agName)
        adtime   = concatenate(adtime)
        

        # save satellite obs to file
        obaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        outDir   = obaseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        util.mk_dir(outDir)
        oprofPath  = outDir  + '/prof.npy'
        oeSurfPath = outDir  + '/eSurf.npy'
        onSurfPath = outDir  + '/nSurf.npy'
        ogvPath    = outDir  + '/gvprcp.npy'
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
        np.save(ogvPath,    agv)
        np.save(onSurfBinPath, anSurfBin)


        np.save(osatelatPath, asatelat)
        np.save(osatelonPath, asatelon)
        np.save(oglatPath,    aglat)   
        np.save(oglonPath,    aglon) 

        np.save(ogNamePath,   agName) 
        np.save(odtimePath,   adtime)

        print aprof.shape
        print oprofPath

