import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from numpy import *
import numpy as np
import myfunc.util as util
import glob
from datetime import datetime, timedelta
import calendar
import h5py
import socket, os, sys
#------------------------------------------
myhost = socket.gethostname()
if myhost == 'shui':
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    figDir   = '/home/utsumi/temp/geo'

elif myhost == 'well':
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()



sateName= 'GPM'
#BBox   = [[20,-130],[55,-55]] # MRMS
BBox   = [[-55,90],[55,200]] # MRMS
[[lllat,lllon],[urlat,urlon]] = BBox
iYM    = [2017,1]
eYM    = [2017,12]
#iYM    = [2014,10]
#eYM    = [2014,10]
lYM    = util.ret_lYM(iYM,eYM)
#iDTime = datetime(2014,5,1)
#eDTime = datetime(2015,5,31)
#iDTime = datetime(2016,4,18)
#eDTime = datetime(2016,4,18)
#eDTime = datetime(2015,5,31)
iabin = 12  # extracted angle bin (start)
eabin = 36  # extracted angle bin (end)

for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    #eDay  =1 # test
    iDTime = datetime(Year,Mon,1,0)
    eDTime = datetime(Year,Mon,eDay,0)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    baseDir = workbaseDir + '/hk02/PMM/NASA/GPM.Ku/2A/V06'

    loid = []
    ly   = []
    lout = []

    #lDTime = lDTime[:2]  # test
    for DTime in lDTime:
        Day = DTime.day
        srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = srcDir + '/2A.GPM.Ku.*.V06A.HDF5'
        lsrcPath= sort(glob.glob(ssearch))
        if len(lsrcPath)==0:
            print 'No Ku file', DTime
            print ssearch
            continue
        for srcPath in lsrcPath:
            oid = int(srcPath.split('.')[-3])
            print Year,Mon,Day,oid
 
            with h5py.File(srcPath, 'r') as h:
                a2lat = h['NS/Latitude'][:,iabin:eabin+1]   # 25 angle bins at the center (Same as Ka) 
                a2lon = h['NS/Longitude'][:,iabin:eabin+1]  # Ku: 49 range bins. Ka: 25 range bins

            a2latTmp = ma.masked_outside(a2lat, lllat, urlat)
            a2lonTmp = ma.masked_outside(a2lon, lllon, urlon)
            a1latmask = a2latTmp.mask.any(axis=1)
            a1lonmask = a2lonTmp.mask.any(axis=1)
   
            a1maskTmp = a1latmask + a1lonmask
            a1y = range(a2lat.shape[0])
            a1y = ma.masked_where(a1maskTmp, a1y).compressed()
            if len(a1y)==0:continue
   
            iy = a1y.min()
            ey = a1y.max()

            with h5py.File(srcPath, 'r') as h:
                a1year = h['NS/ScanTime/Year'][:]
                a1mon  = h['NS/ScanTime/Month'][:]
                a1day  = h['NS/ScanTime/DayOfMonth'][:]
                a1hour = h['NS/ScanTime/Hour'][:]
                a1min  = h['NS/ScanTime/Minute'][:]
                a1sec  = h['NS/ScanTime/Second'][:]

            iyyyy = a1year[iy]
            eyyyy = a1year[ey]
            imm   = a1mon[iy]
            emm   = a1mon[ey]
            idd   = a1day[iy]
            edd   = a1day[ey]
            ihh   = a1hour[iy]
            ehh   = a1hour[ey]
            imin  = a1min[iy]
            emin  = a1min[ey]
            isec  = a1sec[iy]
            esec  = a1sec[ey]

 
            loid.append(oid)
            ly.append([iy,ey])
            ltmp = [Year,Mon,Day,oid,iy,ey,iyyyy,imm,idd,ihh,imin,isec,eyyyy,emm,edd,ehh,emin,esec]
            lout.append(ltmp)
            print Year,Mon,Day,oid, 'Found overpass'

    #-- Save file --
    sout   = util.list2csv(lout)
    outDir = tankbaseDir + '/utsumi/PMM/himawari/obtlist'
    util.mk_dir(outDir)
    outPath= outDir + '/overpass.%s.ABp%03d-%03d.%04d.%02d.csv'%(sateName,iabin, eabin, Year,Mon)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath

    ##'''
    ##-- test draw ----
    #fig = plt.figure(figsize=(8,8))
    #ax  = fig.add_axes([0.2,0.2,0.6,0.6])
    #M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    #
    #for ltmp in lout:

    #    Year,Mon,Day,oid,iy,ey = ltmp[:6]
    #    srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    #    ssearch = srcDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)
    #    srcPath = sort(glob.glob(ssearch))[0]
 
    #    with h5py.File(srcPath, 'r') as h:
    #        a2lat = h['NS/Latitude'][:]             
    #        a2lon = h['NS/Longitude'][:]
    #
    #        a2lat = a2lat[iy:ey+1,:]
    #        a2lon = a2lon[iy:ey+1,:]
    #
    #        M.scatter(a2lon, a2lat, c=a2lat, cmap='jet', s=2)
    #        M.drawcoastlines()
    #
    #figPath = '/home/utsumi/temp/geo/temp.himawri.png'
    #plt.savefig(figPath)
    #print figPath
    #plt.clf()
    ##'''
