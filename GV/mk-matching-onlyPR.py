import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
from datetime import datetime, timedelta
import numpy as np
import myfunc.util as util
import myfunc.IO.GPM as GPM
import GPMGV
import os, sys, socket
import glob
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
import pickle

prdName = 'L2A25'
iYM = [1998,4]
#eYM = [2003,10]
#iYM = [2004,4]
eYM = [2004,10]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()


#ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
ldomain = ['N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA']
#ldomain = ['FLORIDA-KSC']
#ldomain = ['N.Carolina-IPHEx_Duke']
#ldomain = gv.domains

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

res = 0.02
miss= -9999.

nlev  = 40

swathbins= 49
cbins = 11
icbin = (swathbins-1)/2 - (cbins-1)/2

ddtime_1min = timedelta(seconds=60)

hostname = socket.gethostname()
if hostname == 'mizu':
    listDir  = '/home/utsumi/mnt/wellshare/data/GPMGV/sitelist'
    satebaseDir = "/home/utsumi/mnt/wellshare/GPMGV/%s"%(prdName)
    obaseDir    = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25'

elif hostname=='well':
    listDir  = '/media/disk2/share/data/GPMGV/sitelist'
    satebaseDir = "/media/disk2/share/GPMGV/%s"%(prdName)
    obaseDir    = '/media/disk2/share/GPMGV/DOM.L2A25'

#satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'

def ret_dBBoxes():
    #-- read sitelist_summary ---
    listPath = listDir + "/sitelist_summary.csv"
    f = open(listPath, "r"); lines=f.readlines(); f.close()
    dBBox       = {}
    for line in lines[1:]:
        line = line.strip().split(",")
        domain = line[3]
        lllat  = float(line[4])
        urlat  = float(line[5])
        lllon  = float(line[6])
        urlon  = float(line[7])
   
        dBBox[domain]    = [[lllat, lllon],[urlat,urlon]]
    return dBBox


def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)

#*******************************************************
dBBox = ret_dBBoxes()


for domain in ldomain:
    BBox    = dBBox[domain]
    [[lllat,lllon],[urlat,urlon]] = BBox

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
        aprof    = []
        aeSurf   = []
        anSurf   = []
        anSurfBin= []
        asatelat = []
        asatelon = []
        arainType= []
        amethod  = []
        astormH  = []
        afreezH  = []

        adtime   = [] 

        #for satePath in lsatePath:
        for satePath in lsatePath:
            #print satePath
            fileName = os.path.basename(satePath)
            gNum    = fileName.split('.')[2]
            ietime  = fileName.split('.')[1]

            satelatPath = sateDir + '/lat.%s.%s.npy'%(ietime, gNum)
            satelonPath = sateDir + '/lon.%s.%s.npy'%(ietime, gNum)
            rangebinPath= sateDir + '/rangeBinNum.%s.%s.npy'%(ietime, gNum)
            satetimePath= sateDir + '/dtime.%s.%s.npy'%(ietime, gNum)

            eSurfPath   = sateDir + '/eSurf.%s.%s.npy'%(ietime, gNum)
            rainTypePath= sateDir + '/rainType.%s.%s.npy'%(ietime, gNum)
            stormHPath = sateDir  + '/stormH.%s.%s.npy'%(ietime, gNum)
            freezHPath = sateDir  + '/freezH.%s.%s.npy'%(ietime, gNum)
            methodPath = sateDir  + '/method.%s.%s.npy'%(ietime, gNum)

            # load sateprcp
            a3sateprcp = np.load(satePath)
            
            a2satelat  = np.load(satelatPath)
            a2satelon  = np.load(satelonPath)
            a1satetime = np.load(satetimePath)
            a3rangebin = np.load(rangebinPath)

            a2eSurf    = np.load(eSurfPath)
            a2rainType = np.load(rainTypePath)
            a2stormH   = np.load(stormHPath)
            a2freezH   = np.load(freezHPath)
            a2method   = np.load(methodPath)

            #-- constrain angle bin  -----
            a3sateprcp = a3sateprcp[:,icbin:-icbin,:]
            
            a2satelat  = a2satelat[:,icbin:-icbin]
            a2satelon  = a2satelon[:,icbin:-icbin]
            a3rangebin = a3rangebin[:,icbin:-icbin,:]
            a2eSurf    = a2eSurf[:,icbin:-icbin]
            a2rainType = a2rainType[:,icbin:-icbin]
            a2stormH   = a2stormH[:,icbin:-icbin]
            a2freezH   = a2freezH[:,icbin:-icbin]
            a2method   = a2method[:,icbin:-icbin]

            #-- duplicate dtime -----
            a2satetime = np.tile(a1satetime, (cbins,1)).T

            #------------------------

            ##-- extract domain ---
            #nysate,nxsate = a2satelat.shape
            #a2xsate,a2ysate = meshgrid(range(nxsate),range(nysate))

            #a2mskLat = ma.masked_outside(a2satelat, lllat, urlat).mask
            #a2mskLon = ma.masked_outside(a2satelon, lllon, urlon).mask
            #a2msk    = a2mskLat + a2mskLon
            #a1ysate  = ma.masked_where(a2msk, a2ysate).compressed()
            #a1xsate  = ma.masked_where(a2msk, a2xsate).compressed()


            #a2profile= a3sateprcp[a1ysate, a1xsate,:] 
            #a1eSurf  = a2eSurf[a1ysate, a1xsate] 
            #a1rainType= a2rainType[a1ysate, a1xsate] 

            ## bins
            #a1groundbin= a3rangebin[a1ysate, a1xsate,2]
            #a1nsurfbin = a3rangebin[a1ysate, a1xsate,6]

            ## near surface rain
            #a1ytmp   = range(a2profile.shape[0])
            #a1nSurf  = a2profile[a1ytmp,a1nsurfbin]

            ## lat & lon
            #a1satelat= a2satelat[a1ysate, a1xsate]
            #a1satelon= a2satelon[a1ysate, a1xsate]                     


            #   
            ## method, stormH
            #a1method = a2method[a1ysate, a1xsate]
            #a1stormH = a2stormH[a1ysate, a1xsate]
            #
            ## DTime
            #a1dtime  = a1satetime[a1ysate]



            nztemp   = a3sateprcp.shape[2]
            a2profile= a3sateprcp.reshape(-1,nztemp)
            a1eSurf  = a2eSurf.flatten()
            a1rainType= a2rainType.flatten()

            # bins
            a1groundbin= a3rangebin[:, :,2].flatten()
            a1nsurfbin = a3rangebin[:, :,6].flatten()

            # near surface rain
            a1ytmp   = range(a2profile.shape[0])
            a1nSurf  = a2profile[a1ytmp,a1nsurfbin]

            # lat & lon
            a1satelat= a2satelat.flatten()
            a1satelon= a2satelon.flatten()
               
            # method, stormH, freezH
            a1method = a2method.flatten()
            a1stormH = a2stormH.flatten()
            a1freezH = a2freezH.flatten()
            
            # DTime
            a1dtime  = a2satetime.flatten()


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
            amethod.append(a1method)
            astormH.append(a1stormH)
            afreezH.append(a1freezH)

            asatelat.append(a1satelat)
            asatelon.append(a1satelon)
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

        amethod  = concatenate(amethod).astype(int16)
        astormH  = concatenate(astormH).astype(int16)
        afreezH  = concatenate(afreezH).astype(int16)

        asatelat = concatenate(asatelat).astype(float32)
        asatelon = concatenate(asatelon).astype(float32)

        adtime   = concatenate(adtime)

        # save satellite obs to file
        outDir   = obaseDir + '/cbin.%d/%s/%04d%02d'%(cbins, domain, Year,Mon)

        util.mk_dir(outDir)
        oprofPath  = outDir  + '/prof.npy'
        oeSurfPath = outDir  + '/eSurf.npy'
        onSurfPath = outDir  + '/nSurf.npy'
        orainTypePath = outDir + '/rainType.npy'
        onSurfBinPath = outDir + '/nSurfBin.npy'
        omethodPath   = outDir + '/method.npy'
        ostormHPath   = outDir + '/stormH.npy'
        ofreezHPath   = outDir + '/freezH.npy'


        osatelatPath  = outDir + '/sateLat.npy'
        osatelonPath  = outDir + '/sateLon.npy'

        odtimePath    = outDir + '/dtime.npy'

        #---------------------
        np.save(oprofPath, aprof)
        np.save(oeSurfPath, aeSurf)
        np.save(onSurfPath, anSurf)
        np.save(orainTypePath, arainType)
        np.save(onSurfBinPath, anSurfBin)
        np.save(omethodPath, amethod)
        np.save(ostormHPath, astormH)
        np.save(ofreezHPath, afreezH)

        np.save(osatelatPath, asatelat)
        np.save(osatelonPath, asatelon)

        np.save(odtimePath,   adtime)

        print oprofPath


