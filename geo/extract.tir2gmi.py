# %%
import numpy as np
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, subprocess
import os,sys, socket
import calendar
import h5py
from collections import deque
import myfunc.util as util
from bisect import bisect_left

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

lch = [['tir',1],['tir',2],['tir',3],['tir',4],['tir',5],['tir',6],['tir',7],['tir',8],['tir',9],['tir',10],['sir',1],['sir',2]]
#lch = [['tir',1]]
miss = -9999.
ldy   = np.arange(-10,10+1)
ldx   = np.arange(-10,10+1)

ny,nx = len(ldy),len(ldx)
cy,cx = int(ny/2), int(nx/2)
ldydx = [[dy,dx] for dy in ldy for dx in ldx]
iabin = 103  # extracted angle bin (start)
eabin = 117  # extracted angle bin (end)
 
iDTime = datetime(2017,1,1)
eDTime = datetime(2017,12,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

#Year,Mon,Day = 2017,7,1
#oid = 18971
#iy,ey = 5633, 7358
BBox   = [[-55,90],[55,200]]
[[lllat,lllon],[urlat,urlon]] = BBox

lon0 = 85   # Himawari data boundary
lat0 = -60  # Himawari data boundary
dlon = 0.02
dlat = 0.02

#************************
# Read orbit list
#------------------------
sateName = 'GPM.GMI'
iYM = iDTime.timetuple()[:2]
eYM = eDTime.timetuple()[:2]
lYM = util.ret_lYM(iYM,eYM)
lorbit = []
for Year,Mon in lYM:
    listDir = tankbaseDir + '/utsumi/PMM/himawari/obtlist'
    listPath= listDir + '/overpass.%s.ABp%03d-%03d.%04d.%02d.csv'%(sateName,iabin, eabin, Year,Mon)

    f=open(listPath,'r'); lines=f.readlines(); f.close()
    for line in lines:
        line = line.split(',')
        Year,Mon,Day,oid,iy,ey = map(int, line[:6])
        lorbit.append([Year,Mon,Day,oid,iy,ey])

#*******************************
# Start Loop
#-------------------------------
for orbit in lorbit:
    Year,Mon,Day,oid,iy,ey = orbit
    print 'orbit-info=',orbit
    if oid <= 16156: continue  # test
    gmiDir = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = gmiDir + '/1C.GPM.GMI.*.%06d.V05A.HDF5'%(oid)

    lgmiPath= glob.glob(ssearch)
    if len(lgmiPath)==0:
        print 'No file'
        print ssearch
        continue

    gmiPath = lgmiPath[0]


    outDir = tankbaseDir + '/utsumi/PMM/himawari/obt.gmi.ABp%03d-%03d/%04d/%02d/%02d/%06d'%(iabin, eabin, Year,Mon,Day,oid)
    util.mk_dir(outDir)
 
    #*******************************
    # Read GMI lat, lon
    #-------------------------------
   
    print ''
    print 'Read GMI'
    with h5py.File(gmiPath, 'r') as h:
        a2lat = h['S1/Latitude'][:,iabin:eabin+1]   #
        a2lon = h['S1/Longitude'][:,iabin:eabin+1]  #
        a1year = h['S1/ScanTime/Year'][:]
        a1mon  = h['S1/ScanTime/Month'][:]
        a1day  = h['S1/ScanTime/DayOfMonth'][:]
        a1hour = h['S1/ScanTime/Hour'][:]
        a1min  = h['S1/ScanTime/Minute'][:]
    
    a1min10 = (a1min/10.).astype('int32')*10
    
    print 'Done GMI'
    print 'make lDTime'
    a1dtimeOrg = np.array([datetime(yyyy,mm,dd,hh,mn)
                for (yyyy,mm,dd,hh,mn) in
                zip(a1year, a1mon, a1day, a1hour, a1min10)
                ])
    
    #--- X and Y ----
    nyOrg, nxOrg = a2lat.shape
    a2x,a2y = np.meshgrid(np.arange(iabin,eabin+1).astype('int32'), np.arange(nyOrg).astype('int32'))
    
    a2lat  = a2lat[iy:ey+1]
    a2lon  = a2lon[iy:ey+1]
    a2x    = a2x[iy:ey+1]
    a2y    = a2y[iy:ey+1]
    
    a1dtime= a1dtimeOrg[iy:ey+1]
    
    setdtime = np.sort(list(set(a1dtime)))
    print setdtime
    ltmpDir = []
    for itime,dtime0 in enumerate(setdtime):
        yyyy,mm,dd,hh,mn = dtime0.timetuple()[:5]
        dtime1 = dtime0 + timedelta(minutes=10)
        y0 = bisect_left(a1dtime, dtime0)
        y1 = bisect_left(a1dtime, dtime1) -1
        
        a2latTmp  = a2lat[y0:y1+1]
        a2lonTmp  = a2lon[y0:y1+1]
        a2xTmp    = a2x[y0:y1+1]
        a2yTmp    = a2y[y0:y1+1]
    
        a1lat  = a2latTmp.flatten()
        a1lon  = a2lonTmp.flatten()
        a1xgmi = a2xTmp.flatten()   
        a1ygmi = a2yTmp.flatten()   
        a1dtimeTmp = a1dtimeOrg[a1ygmi]

        #--- Save GMI variables ----------
        np.save(outDir + '/lat.%02d.npy'%(itime), a2lat)
        np.save(outDir + '/lon.%02d.npy'%(itime), a2lon)
        np.save(outDir + '/dtime.%02d.npy'%(itime), a1dtimeTmp)


        #--- Corresponding GEO x and y ----
        a1xgeo = ((a1lon - lon0)/dlon).astype('int32')
        a1ygeo = ((a1lat - lat0)/dlat).astype('int32')
    
        for (ch,chnum) in lch:
    
            #-- count2tbb ---------
            prog = '/home/utsumi/bin/PMM/geo/count2tbb.py'
            scmd = 'python %s %s %s %s %s %s %s %s'%(prog,yyyy,mm,dd,hh,mn, ch, chnum)
            print scmd
            subprocess.call(scmd.split(' '))
    
            #-- Extract geo data at DPR pixel --
            timestamp = '%04d%02d%02d%02d%02d'%(yyyy,mm,dd,hh,mn)
            tmpDir = tankbaseDir + '/utsumi/PMM/himawari/temp/temp.%s.%s.%02d'%(timestamp,ch,chnum)
            ltmpDir.append(tmpDir)
    
            geoPath = tmpDir + '/grid20.dat'
            if not os.path.exists(geoPath):
                a2geo = np.ones([6000,6000],'float32')*(-9999.)
            else:
                a2geo  = np.flipud(np.fromfile(geoPath, 'float32').reshape(6000,6000))
    
            nl = len(a1ygeo)
            a3dat = np.empty([nl,ny,nx], float32)
            for [dy,dx] in ldydx:
                a1geo = a2geo[a1ygeo+dy, a1xgeo+dx]
                a3dat[:,cy+dy,cx+dx] = a1geo
    

            #-- Reshape output ---
            nygmi, nxgmi = a2latTmp.shape
            a4dat = a3dat.reshape(nygmi,nxgmi,ny,nx)

            #k=20

            #a2xtmp2,a2ytmp2 = np.meshgrid(np.arange(nxgmi),np.arange(nygmi))
            #a1xtmp2 = a2xtmp2.flatten()
            #a1ytmp2 = a2ytmp2.flatten()
            #ytmp,xtmp = a1ytmp2[k], a1xtmp2[k]
            #print ytmp,xtmp
            #print ''
            #print a3dat[k]
            #print a3dat[k].shape
            #print ''
            #print a4dat[ytmp,xtmp]
            #print a4dat[ytmp,xtmp].shape
            #sys.exit()

            #-- Save -------------
            outPath= outDir + '/%s.%02d.%02d.npy'%(ch,chnum,itime)
            np.save(outPath, a4dat.astype('float32'))
            print ''
            print 'Extracted'
            print outPath
    
    #-- Connect DPR files ------------
    for var in ['lat','lon','dtime']:
        for itime in range(len(setdtime)):
            srcPath = outDir + '/%s.%02d.npy'%(var,itime)
            if itime ==0:
                a1out = np.load(srcPath, allow_pickle=True)
            else:
                a1out = np.concatenate([a1out, np.load(srcPath, allow_pickle=True)], axis=0)  # set allow_pickle for loading dtime array
    
            os.remove(srcPath)
    
        outPath = outDir + '/%s.npy'%(var)
        np.save(outPath, a1out)
        print ''
        print outPath
        print a1out.shape

    #-- Connect GEO files ------------
    for [ch,chnum] in lch:
        a3out = None
        for itime in range(len(setdtime)):
            srcPath= outDir + '/%s.%02d.%02d.npy'%(ch,chnum,itime)
    
            if not os.path.exists(srcPath):
                continue
     
            if a3out is None:
                a3out = np.load(srcPath)
            else:
                a3out = np.concatenate([a3out, np.load(srcPath)],axis=0)
    
            os.remove(srcPath)
    
        outPath= outDir + '/%s.%02d.npy'%(ch,chnum)
        np.save(outPath, a3out)
        print ''
        print 'Connected'
        print outPath
        print a3out.shape

    #-- Save iy and ey ----
    outPath= outDir + '/iyey.npy'
    aout = np.array([iy,ey]).astype('int32')
    np.save(outPath, aout)


    ##-- Remove temporay files ----- 
    #for tmpDir in ltmpDir:
    #    ch    = tmpDir.split('.')[-2]
    #    chnum = int(tmpDir.split('.')[-1])
    #    print ''
    #    print 'Remove files'
    #    ltmpPath = glob.glob(tmpDir + '/*.geoss*')
    #    for tmpPath in ltmpPath: os.remove(tmpPath)
    
    #    tmpPath = tmpDir + '/grid20.dat'
    #    if os.path.exists(tmpPath):
    #        os.remove(tmpPath) 
    
    #    tmpPath = tmpDir + '/%s.%02d'%(ch, chnum)
    #    if os.path.exists(tmpPath):
    #        os.remove(tmpPath) 
    
    
    #    os.rmdir(tmpDir) 
    



# %%
