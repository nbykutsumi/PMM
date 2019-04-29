#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

from numpy import *
import numpy as np
import myfunc.util as util
import glob
from datetime import datetime, timedelta
import calendar
import h5py
import socket, os, sys
#------------------------------------------
#-- Read GTOPO ----
demPath = '/work/hk01/utsumi/gtopo/0.1deg/dem.lat.-90.090.lon.-180.0180.1800x3600.npy'
a2dem   = np.load(demPath)
latdem0 = -90
londem0 = -180
dlatdem = 0.1
dlondem = 0.1
nydem = int(180/dlatdem)
nxdem = int(360/dlondem)
sateName= 'GPM'
BBox   = [[20,-130],[55,-55]] # MRMS
[[lllat,lllon],[urlat,urlon]] = BBox
iYM    = [2014,12]
eYM    = [2015,5]
#iYM    = [2014,10]
#eYM    = [2014,10]
lYM    = util.ret_lYM(iYM,eYM)
#iDTime = datetime(2014,5,1)
#eDTime = datetime(2015,5,31)
#iDTime = datetime(2016,4,18)
#eDTime = datetime(2016,4,18)
#eDTime = datetime(2015,5,31)


for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1,0)
    eDTime = datetime(Year,Mon,eDay,0)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    baseDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'

    loid = []
    ly   = []
    lout = []

    #lDTime = lDTime[:2]  # test
    for DTime in lDTime:
        Day = DTime.day
        srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        lsrcPath= sort(glob.glob(srcDir + '/1C.GPM.GMI.XCAL2016-C.*.HDF5'))
    
        for srcPath in lsrcPath:
            oid = int(srcPath.split('.')[-3])
            print Year,Mon,Day,oid
 
            #-- Known missing orbits ---
            if oid in range(3683,3717+1): continue # 2014/10/22-24. No GMI data
            #---------------------------
            #if oid != 4447: continue  # test
 
            with h5py.File(srcPath, 'r') as h:
                a2lat = h['S1/Latitude'][:]             
                a2lon = h['S1/Longitude'][:]

            a2latTmp = ma.masked_outside(a2lat, lllat, urlat)
            a2lonTmp = ma.masked_outside(a2lon, lllon, urlon)
            a1latmask = a2latTmp.mask.any(axis=1)
            a1lonmask = a2lonTmp.mask.any(axis=1)
   
            a1y = range(a2lat.shape[0])
    
            a1maskTmp = a1latmask + a1lonmask
            a1yTmp = ma.masked_where(a1maskTmp, a1y).compressed()
            if len(a1yTmp)==0:continue
    
            #-- Mask for Ocean -------
            ny,nx  = a2lat.shape
            X,Y    = np.meshgrid(range(nx),range(ny))
    
            a1saty = ((a2lat - latdem0)/dlatdem).astype(int16).flatten()
            a1satx = ((a2lon - londem0)/dlondem).astype(int16).flatten()
            a1satx = ma.masked_equal(a1satx,nxdem).filled(0)  # for the case lon==180.0
            # For missing lat&lon info
            a2latmissmask = ma.masked_equal(a2lat, -9999.9).mask
            a1latmissmask = a2latmissmask.flatten()
            if np.any(a1latmissmask):
                a1saty = ma.masked_where(a1latmissmask, a1saty).filled(0)
                a1satx = ma.masked_where(a1latmissmask, a1satx).filled(0)
            ##-- test ---
            #idxtmp = np.argmax(a1saty)
            #print idxtmp,a2lat.flatten()[idxtmp], a1saty[idxtmp]
            ##----------- 
            a2elev = a2dem[a1saty, a1satx].reshape(ny,nx)

            if ma.isarray(a1latmissmask):
                a2elev = ma.masked_where(a2latmissmask, a2elev).filled(-9999) 

 
            a2oceanmask = ma.masked_less_equal(a2elev, 0)
            a1oceanmask = a2oceanmask.mask.all(axis=1)
    
            #-- END Mask for Ocean -------
            a1mask = a1latmask + a1lonmask + a1oceanmask
            a1y = ma.masked_where(a1mask, a1y).compressed()
   
            if len(a1y)==0: continue
 
            iy = a1y.min()
            ey = a1y.max()
    
            loid.append(oid)
            ly.append([iy,ey])
            ltmp = [Year,Mon,Day,oid,iy,ey]
            lout.append(ltmp)
            print Year,Mon,Day,oid, 'Found overpass'

    #-- Save file --
    sout   = util.list2csv(lout)
    outDir = '/work/hk01/utsumi/PMM/US/obtlist'
    outPath= outDir + '/overpass.%s.%04d.%02d.csv'%(sateName,Year,Mon)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath

    '''
    #-- test draw ----
    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_axes([0.2,0.2,0.6,0.6])
    M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    
    for ltmp in lout:

        Year,Mon,Day,oid,iy,ey = ltmp
        srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        srcPath = sort(glob.glob(srcDir + '/1C.GPM.GMI.XCAL2016-C.*.%06d.*.HDF5'%(oid)))[0]
    
        with h5py.File(srcPath, 'r') as h:
            a2lat = h['S1/Latitude'][:]             
            a2lon = h['S1/Longitude'][:]
    
            a2lat = a2lat[iy:ey+1,:]
            a2lon = a2lon[iy:ey+1,:]
    
            M.scatter(a2lon, a2lat, c=a2lat, cmap='jet', s=2)
            M.drawcoastlines()
    
    figPath = '/home/utsumi/temp/temp.us.png'
    plt.savefig(figPath)
    print figPath
    plt.clf()
    '''
