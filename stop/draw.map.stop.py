import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util


myhost = socket.gethostname()
if myhost == 'shui':
    workbaseDir = '/work'
    tankbaseDir = '/tank'
    figDir   = '/home/utsumi/temp/stop'

elif myhost == 'well':
    workbaseDir = '/home/utsumi/mnt/lab_work'
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()


reftype = 'dpr'
#reftype = 'mrms'


## Illinois (RNS)
#oid = 1456
#Year,Mon,Day = 2014,6,1
#iy, ey = -9999,-9999
##iy,ey = 1311, 1351
#clat    = 28.3
#clon    = 133.0
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'glb.v03.minrec1000.maxrec10000'
##expr = 'glb.stop-wgt-obs-01.minrec1000.maxrec10000'

## Amazon (Underestimation)
#oid = 1732
#Year,Mon,Day = 2014,6,18
#iy, ey = -9999,-9999
##iy,ey = 1311, 1351
#clat    = 2
#clon    = 303 -360.
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'glb.minrec1000.maxrec10000'


## West of Great lakes (Underestimation)
#oid = 1871
#Year,Mon,Day = 2014,6,27
#iy, ey = -9999,-9999
##iy,ey = 1311, 1351
#clat    = 37
#clon    = 270 -360.
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'glb.minrec1000.maxrec10000'


## West of Great lakes (Underestimation)
#oid = 1917
#Year,Mon,Day = 2014,6,30
##iy, ey = -9999,-9999
#iy,ey = 1826, 1886
#clat    = 42
#clon    = 269 -360.
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'test'
##expr = 'glb.minrec1000.maxrec10000'


## Alaska to check (Heavy rain missing)
#oid = 1700
#Year,Mon,Day = 2014,6,16
##iy, ey = -9999,-9999
#iy,ey = 1311, 1351
#clat    = 57
#clon    = 211 -360.
#DB_MAXREC = 10000
##DB_MAXREC = 20000
#DB_MINREC = 1000
##expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
##expr = 'glb.nprof'

## New York to check (good)
#oid = 1496
#Year,Mon,Day = 2014,6,3
##iy, ey = -9999,-9999
#iy, ey = 1054, 1134
#clat    = 42
#clon    = 286 -360.
#DB_MAXREC = 10000
##DB_MAXREC = 20000
#DB_MINREC = 1000
#expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

## Canada to check batch-version (Under est)
#oid = 1686
#Year,Mon,Day = 2014,6,16
##iy, ey = -9999,-9999
#iy, ey = 1720,1780
#clat    = 52.0
#clon    = 270.0 -360.
#DB_MAXREC = 10000
##DB_MAXREC = 20000
#DB_MINREC = 1000
##expr  = 'test.wtb'
#expr  = 'test.ntb'
##expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

## Yellow sea to check batch-version
#oid = 1693
#Year,Mon,Day = 2014,6,16
#iy, ey = -9999,-9999
#clat    = 32.0
#clon    = 123.0
#DB_MAXREC = 10000

## Myanmar to check batch-version
#oid = 1688
#Year,Mon,Day = 2014,6,16
#iy, ey = -9999,-9999
#clat    = 26.7
#clon    = 98.2
#DB_MAXREC = 10000

## Europe just to check batch-version
#oid = 1927
#Year,Mon,Day = 2014,7,1
#iy, ey = 1764, 1784
#clat    = 50.8
#clon    = 28.1
#DB_MAXREC = 10000
##DB_MAXREC = 20000


## West of Great Lakes (good)
#oid = 1748
#Year,Mon,Day = 2014,6,20
#iy, ey = -9999, -9999
#clat    = 42
#clon    = -92
#DB_MAXREC = 20000

## West of Great Lakes (Under est.)
#oid = 1574
#Year,Mon,Day = 2014,6,8
#iy, ey = -9999,-9999
#clat    = 41
#clon    = -(360-254)
#DB_MAXREC = 10000
##DB_MAXREC = 20000
#expr = 'glb.minrec1000.maxrec%d'%(DB_MAXREC)

## Southeast of Grate lakes (Fair)
#oid = 1573
#Year,Mon,Day = 2014,6,8
#iy, ey = -9999, -9999
#clat    = 41
#clon    = -(360-283)
#DB_MAXREC = 20000



#### Great Lake Snow
#oid = 4140
#Year,Mon,Day = 2014,11,20
#iy, ey = 1145,1205
#clat    = 43
#clon    = -79.5
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

## Great Lake Snow (some missing)
#oid = 4914
#Year,Mon,Day = 2015,1,9
#iy, ey = 1755, 1815
#clat    = 43
#clon    = -79.5
#DB_MAXREC = 10000
#DB_MINREC = 1000
#expr = 'glb.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

## SE.US case, oid=003556, 2014/10/14
#oid = 3556
#Year,Mon,Day = 2014,10,14
#iy, ey = -9999,-9999
##iy, ey = 917,1117
#clat    = 34    # SE.US case. oid = 003556
#clon    = -86   # 2014/10/14  05:42:03 UTC
#DB_MINREC = 1000
#DB_MAXREC = 10000
#expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#stopexpr = 'best01'
#stopact  = 'HTQZ'


## East of Japan case
#oid = 1963
#Year,Mon,Day = 2014,7, 3
#iy, ey = -9999,-9999
##iy, ey = 917,1117
#clat    = 34   
#clon    = 141  
#DB_MINREC = 1000
#DB_MAXREC = 10000
#expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#stopexpr = 'best01'
#stopact  = 'HTQZ'

# South of Japan case
oid = 2532
Year,Mon,Day = 2014,8,9
iy, ey = -9999,-9999
#iy, ey = 917,1117
clat    = 32 
clon    = 134
DB_MINREC = 1000
DB_MAXREC = 10000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
stopexpr = 'best01'
stopact  = 'HTQZ'




## SW.Japan typhoon case, oid=019015, 2017/07/03
#oid = 19015
#Year,Mon,Day = 2017,7,3
#iy, ey = 1793, 1973
#clat    = 33    # SE.US case. oid = 003556
#clon    = 130   # 2014/10/14  05:42:03 UTC
#DB_MAXREC = 20000

## SW.Japan typhoon case, oid=019015, 2017/07/03
#oid = 16166
#Year,Mon,Day = 2017,1,1
##iy, ey = -9999,-9999
#iy, ey = 2246, 2446
#clat    = -10    # SE.US case. oid = 003556
#clon    = -75   # 2014/10/14  05:42:03 UTC
#DB_MAXREC = 20000

## E of Hokkaido, Japan case oid=004305, 2014/12/1
#oid = 4305
#Year,Mon,Day = 2014,12,1
#iy, ey = 1536,2136
#clat    = 45    # SE.US case. oid = 003556
#clon    = 143   # 2014/10/14  05:42:03 UTC
#DB_MAXREC = 20000

## Japan sea case oid=004300, 2014/12/1
#oid = 4300
#Year,Mon,Day = 2014,12,1
#iy, ey = -9999,-9999 
#clat    = 41   
#clon    = 137
#DB_MAXREC = 20000

## Minnesota case oid=002081  2014/7/11
#oid = 2081
#Year,Mon,Day = 2014,7,11
#iy, ey = -9999,-9999 
#clat    = 46 
#clon    = -92
#DB_MAXREC = 20000

## QJRMS case, oid=012149, 2016/4/18
#oid = 12149
#iy, ey = 1038, 1098
#clat    = 32. # QJRMS case, oid=012149
#clon    = -94 # -180 - +180
#DB_MAXREC = 20000

thpr = 0.2

if clat !=-9999.:
    #dlatlon = 20
    #dlatlon = 10
    dlatlon = 4
    #dscan   = 100
    dscan   = 200  # used to extract GPROf
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
    #BBox    = [[clat-2, clon-4],[clat+2,clon+4]] # test
else:
    BBox    = [[-60,-180],[60,180]]

#BBox = [[41,-84], [45,-75]]   # test

[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.

#*****************************
# Lat & Lon of GMI scans
#-----------------------------
srcDir = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch = srcDir + '/1C.GPM.GMI.*.%06d.V05A.HDF5'%(oid)
srcPath= glob.glob(ssearch)[0]
with h5py.File(srcPath,'r') as h:
    a2latpmw = h['S1/Latitude'][:]
    a2lonpmw = h['S2/Longitude'][:]

if clat != -9999.:
    _, a2latpmw, a2lonpmw, iyTmp, eyTmp = epcfunc.extract_domain_2D(a2latpmw, a2latpmw, a2lonpmw, clat, clon, dlatlon, dscan, returnidx=True)


#*****************************
# Read ML storm top data
#-----------------------------
srcDir = tankbaseDir + '/utsumi/PMM/stop/orbit/%s-%s-ssn0/%04d/%02d/%02d'%(stopexpr, stopact, Year, Mon, Day)
srcPath= srcDir + '/stop.%06d.npy'%(oid)

a2mlstop = np.load(srcPath)
if clat != -9999.:
    a2mlstop = a2mlstop[iyTmp:eyTmp+1]

#**********************************
# Read EPC-top storm top data and surface precip
#**********************************
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

srcDir = tankbaseDir + '/utsumi/PMM/retepc/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
srcPath= srcDir + '/top-heightStormTopNS.%s.npy'%(stamp)
a2epcstop = np.load(srcPath)

srcPath= srcDir + '/nsurfNScmb.%s.npy'%(stamp)
a2epcprec = np.load(srcPath)

if clat != -9999.:
    a2epcstop = a2epcstop[iyTmp:eyTmp+1]
    a2epcprec = a2epcprec[iyTmp:eyTmp+1]


#**********************************
#- Read DPR/CMB data 
#**********************************
srcDir = workbaseDir + '/hk02/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch= srcDir + '/2B.GPM.DPRGMI.CORRA*.%06d.V06A.HDF5'%(oid)
dprPath= glob.glob(ssearch)[0]
with h5py.File(dprPath,'r') as h:
    a2dprprec = h['NS/surfPrecipTotRate'][:]
    a2dprstop = h['NS/Input/stormTopAltitude'][:] 
    a2latdpr  = h['NS/Latitude'][:]
    a2londpr  = h['NS/Longitude'][:]


#********************************
#-- Draw figure ---
print 'Draw figure'
fig   = plt.figure(figsize=(7,7))

for i in range(4):
    if i in [2,3]:
        ssize=0.2
    else:
        #ssize=0.5
        ssize=0.9

    if i == 3:
        mycm  = 'gist_ncar_r'
        vmin, vmax = 1,20

    else:
        mycm  = 'jet'
        vmin, vmax = 2,10

    if i==0:
        #-- ML stop -----
        ax = fig.add_axes([0.05,0.1,0.35,0.35])
        a2dat = ma.masked_where(a2epcprec<thpr, a2mlstop) /1000.  # [km]
        a2lat = a2latpmw
        a2lon = a2lonpmw
        stype = 'FNN' + ' (%s %s)'%(stopexpr, stopact)
        cax = fig.add_axes([0.41,0.13,0.01,0.3])

    elif i==1:
        #-- EPC stop -----
        ax = fig.add_axes([0.55,0.1,0.35,0.35])
        a2dat = ma.masked_where(a2epcprec<thpr, a2epcstop) / 1000. # [km]
        a2lat = a2latpmw
        a2lon = a2lonpmw
        stype = 'EPC Storm Top' 

        cax = fig.add_axes([0.91,0.13,0.01,0.3])

    elif i==2:
        ax = fig.add_axes([0.05,0.55,0.35,0.35])
        a2dat = ma.masked_where(a2dprprec<thpr, a2dprstop) / 1000. # [km]
        a2lat = a2latdpr
        a2lon = a2londpr
        stype = 'Ku StormTopHeight'
        cax = fig.add_axes([0.41,0.58,0.01,0.3])
    
    elif i==3:
        ax = fig.add_axes([0.55,0.55,0.35,0.35])
        a2dat = ma.masked_where(a2dprprec<thpr, a2dprprec)
        a2lat = a2latdpr
        a2lon = a2londpr
        stype = 'Precip. (CMB)'
        cax = fig.add_axes([0.91,0.58,0.01,0.3])
      
    M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    im    = M.scatter(a2lon, a2lat, c=a2dat, cmap=mycm, s=ssize, vmin=vmin, vmax=vmax)
    
    #M.bluemarble()
    M.drawcoastlines(color='k',linewidth=1)
    
    dgrid      = 5
    parallels  = arange(-90,90, dgrid)
    meridians  = arange(-180,180,dgrid)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d', color='0.8')
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d', color='0.8')
    
    ax.set_title(stype)

    #-- DPR edge ---- 
    a1lat0 = a2latdpr[:,0]
    a1lon0 = a2londpr[:,0]
    a1lat1 = a2latdpr[:,-1]
    a1lon1 = a2londpr[:,-1]
    ax.plot(a1lon0,a1lat0,'-',color='k', linewidth=0.8)
    ax.plot(a1lon1,a1lat1,'-',color='k', linewidth=0.8)




    #-- Colorbar (Shared, Storm top) ----
    plt.colorbar(im, orientation='vertical', cax=cax)
    


#-- Suptitle -------------
ssuptitle = '%04d/%02d/%02d #%06d'%(Year,Mon,Day,oid)
plt.suptitle(ssuptitle, fontsize=12)
##------------
util.mk_dir(figDir)
outPath  = figDir + '/stop.map.%s-%s.%06d.y%d-%d.png'%(stopexpr,stopact,oid,iyTmp,eyTmp)
plt.savefig(outPath)
print outPath
