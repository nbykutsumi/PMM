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


iYM = [2014,1]
eYM = [2014,1]
lYM = util.ret_lYM(iYM, eYM)


#ldomain = gv.domains
ldomain = ['FLORIDA-STJ']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

res = 0.02
miss= -9999.
offset_bef = 5    # minutes
offset_aft = 30   # minutes

mindist = 10 # km
lgvtype = ['All','Gd']

ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()

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

#*******************************************************
dnynx, dBBox, dBBoxBnd = ret_dBBoxes(res)

for domain in ldomain:
    nymap,nxmap = dnynx[domain]
    BBox    = dBBox[domain]
    BBoxBnd = dBBoxBnd[domain]
    [[lllatBnd, lllonBnd], [urlatBnd,urlonBnd]] = dBBoxBnd[domain]

    a2gvzero = zeros([nymap,nxmap]).astype(float32)

    for YM in lYM:
        Year, Mon = YM

        if (domain,Year,Mon) not in dgName.keys():
            continue
        #-- load minimum distance file for the month
        a2mindist = load_a2mindist(domain, Year, Mon)
        a2radmask = ma.masked_greater(a2mindist, mindist).mask

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

        agvAll = []
        agvGd  = []
        for satePath in lsatePath:
            print satePath
            fileName = os.path.basename(satePath)
            gNum    = fileName.split('.')[2]
            ietime  = fileName.split('.')[1]
            #itime   = ietime.split('-')[0]
            #etime   = ietime.split('-')[1]

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



                # index arrays
                '''independent from dMnt in the following loop'''

                a1xmap, a1ymap, a1xsate, a1ysate = ret_aidx(a2satelatTmp, a2satelonTmp, lllatBnd, lllonBnd, nymap, nxmap, a2radmask=a2radmask) # Use monthly a2radmask

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

                ## lat & lon for test
                #a1satelatTmp = a2satelatTmp[a1ysate, a1xsate]
                #a1satelonTmp = a2satelonTmp[a1ysate, a1xsate]

                #for itmp, esurf in enumerate(a1eSurf):
                #    if esurf <5: continue
                #    print a1satelatTmp[itmp], a1satelonTmp[itmp], esurf

                # profile of 20 range bins above ground
                a2profout= empty([a2profile.shape[0], 20]).astype(int16)
                for iy in range(a2profile.shape[0]):
                    groundbin     = a1groundbin[iy]
                    if (20<=groundbin)&(groundbin <=80):
                        ''' set 20 levels above ground. Note that range bin with ground level is not included. '''
                        a1tmp  = a2profile[iy, groundbin-20:groundbin]
                        a2profout[iy] = a1tmp[::-1] # flip order. New order: Low level to High level.

                    else:
                        a2profout[iy] = -8888

                # append to container
                aeSurf.append(a1eSurf)
                anSurf.append(a1nSurf)
                aprof.append(a2profout)

                # corrected nSurfBin 
                anSurfBin.append(a1groundbin-a1nsurfbin-1)

                #--------------------------------
                # load and extract gauge data
                #--------------------------------
                ldMnt     = arange(-offset_bef, offset_aft +1).astype(int32)
                for gvtype in lgvtype:
                    a2gvprcp  = empty([len(a1eSurf), len(ldMnt)]).astype(float32)
                    for iMnt, dMnt in enumerate(ldMnt):
    
                        dtimeGV = dtime0 + dMnt*ddtime_1min + ddtime_1min
                        YearGV = dtimeGV.year
                        MonGV  = dtimeGV.month
                        DayGV  = dtimeGV.day
                        HourGV = dtimeGV.hour
                        MntGV  = dtimeGV.minute

                        # load gvmap
                        gvmapDir = gvmapbaseDir + '/%s/%04d%02d'%(domain,YearGV,MonGV)
                        gvmapPath = gvmapDir + '/map.%s.%s.%04d%02d%02d.%02d%02d.npy'%(gvtype, domain,YearGV,MonGV,DayGV,HourGV,MntGV)

 
                        a2gvmap = np.load(gvmapPath)
                        if a2gvmap.shape==(1,):
                            a2gvmap = a2gvzero.copy() 
                        
                        # extract from gvmap
                        a1gvprcp = a2gvmap[a1ymap, a1xmap]
   
 
                        a2gvprcp[:,iMnt] = a1gvprcp 

                    # append gvprcp
                    if   gvtype=='All':
                        agvAll.append(a2gvprcp)

                    elif gvtype=='Gd':
                        agvGd.append(a2gvprcp)

                    else:
                        print 'check gvtype',gvtype
                        sys.exit()


        # make output array
        aprof   = concatenate(aprof, axis=0).astype(int16)
        aeSurf  = concatenate(aeSurf).astype(float32)
        anSurf  = concatenate(anSurf).astype(int16)
        anSurfBin  = concatenate(anSurfBin).astype(int16)

        agvAll  = concatenate(agvAll, axis=0).astype(float32)
        agvGd   = concatenate(agvGd,  axis=0).astype(float32)

        # save satellite obs to file
        obaseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        outDir   = obaseDir + '/%s/%04d%02d'%(domain, Year,Mon)
        util.mk_dir(outDir)
        oprofPath  = outDir  + '/prof.npy'
        oeSurfPath = outDir  + '/eSurf.npy'
        onSurfPath = outDir  + '/nSurf.npy'
        ogvAllPath= outDir  + '/gvAll.npy'
        ogvGdPath = outDir  + '/gvGd.npy'
        onSurfBinPath = outDir  + '/nSurfBin.npy'

        #---------------------

        np.save(oprofPath, aprof)
        np.save(oeSurfPath, aeSurf)
        np.save(onSurfPath, anSurf)
        np.save(ogvAllPath, agvAll)
        np.save(ogvGdPath,  agvGd)
        np.save(onSurfBinPath, anSurfBin)
        print oprofPath

        ##''' 
        #aidxtmp= range(aprof.shape[0]) 
        #print aprof[aidxtmp,anSurfBin+1]
        #figDir = '/work/a01/utsumi/GPMGV/fig'
        #figPath = figDir + '/tmp.prof.%04d.%02d.png'%(Year,Mon)
        #plt.imshow(ma.masked_less(aprof.T, 0),origin='lower')
        #plt.colorbar()

        #plt.savefig(figPath)
        #plt.clf()
        #print figPath
        ##''' 
