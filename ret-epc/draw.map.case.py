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
    #srcDir = '/home/utsumi/temp/out'
    #srcDir = '/home/utsumi/temp/out/my'
    srcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    workbaseDir = '/work'
    mrmsDir  = '/work/hk01/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    mrmsDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()


#reftype = 'dpr'
reftype = 'mrms'

## Africa case
#oid = 2421
#iy, ey = 2029, 2089
#clat    = 14 # Africa case
#clon    = 2  # -180 - +180

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
#iy, ey = 1746, 1766
#clat    = 52.0
#clon    = 270.0 -360.
#DB_MAXREC = 10000
##DB_MAXREC = 20000
#DB_MINREC = 1000
#expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

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



### Great Lake Snow
oid = 4140
Year,Mon,Day = 2014,11,20
iy, ey = 1145,1205
clat    = 43
clon    = -79.5
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

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
#iy, ey = 987, 1047
##iy, ey = 1012, 1022
#clat    = 34    # SE.US case. oid = 003556
#clon    = -86   # 2014/10/14  05:42:03 UTC
#DB_MAXREC = 10000
#expr = 'test.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

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

thpr = 0.1  # Minimum threshold for EPC
#thpr = 0.01  # Minimum threshold for EPC
#thpr = 0.  # Minimum threshold for EPC

if clat !=-9999.:
    #dlatlon = 20
    dlatlon = 5
    #dscan   = 55
    dscan   = 200  # used to extract GPROf
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
    #BBox    = [[clat-2, clon-4],[clat+2,clon+4]] # test
else:
    BBox    = [[-60,-180],[60,180]]

#BBox = [[41,-84], [45,-75]]   # test

[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.


#stamp = '%06d.y%04d-%04d.nrec%d'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

srcDir         = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)

#srcDir = '/home/utsumi/temp/out/my'   # test
nsurfMSPath    = srcDir + '/nsurfMS.%s.npy'%(stamp)
nsurfNSPath    = srcDir + '/nsurfNS.%s.npy'%(stamp)
nsurfMScmbPath = srcDir + '/nsurfMScmb.%s.npy'%(stamp)
nsurfNScmbPath = srcDir + '/nsurfNScmb.%s.npy'%(stamp)
#topnsurfMScmbPath = srcDir + '/top-nsurfMScmb.%s.npy'%(stamp)
#topnsurfNScmbPath = srcDir + '/top-nsurfNScmb.%s.npy'%(stamp)

print nsurfMSPath

latPath   = srcDir + '/lat.%s.npy'%(stamp)
lonPath   = srcDir + '/lon.%s.npy'%(stamp)
#
#-- GPROF --
gprofDir = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid)
gprofPath= glob.glob(ssearch)[0]

#-- MRMS --
ssearch  = mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.%04d%02d%02d.%06d.*.npy'%(Year,Mon,Day,oid)
print ssearch
try:
    mrmsPath = glob.glob(ssearch)[0]
except:
    mrmsPath = 'xx'
#*****************
#- Read My data ----
print 'Read my data'
#a2MS    = np.load(nsurfMSPath)
#a2NS    = np.load(nsurfNSPath)
a2MS    = np.load(nsurfMSPath)
a2NS    = np.load(nsurfNSPath)
a2MScmb = np.load(nsurfMScmbPath)
a2NScmb = np.load(nsurfNScmbPath)
#a2topMScmb = np.load(topnsurfMScmbPath)
#a2topNScmb = np.load(topnsurfNScmbPath)

a2latMy = np.load(latPath)
a2lonMy = np.load(lonPath)

#-- Screen <thpr mm/h
a2MS= ma.masked_less(a2MS,thpr)
a2NS= ma.masked_less(a2NS,thpr)
a2MScmb= ma.masked_less(a2MScmb,thpr)
a2NScmb= ma.masked_less(a2NScmb,thpr)

#*****************
#- Read GPROF data ----
print 'Read GPROF'
with h5py.File(gprofPath) as h:
    a2esurfgp = h['S1/surfacePrecipitation'][:]
    a2latgpOrg= h['S1/Latitude'][:]
    a2longpOrg= h['S1/Longitude'][:]

if clat != -9999.:
    a2esurfgp, a2latgp, a2longp, iygp, eygp = epcfunc.extract_domain_2D(a2esurfgp, a2latgpOrg, a2longpOrg, clat, clon, dlatlon, dscan, returnidx=True)

else:
    a2latgp = a2latgpOrg
    a2longp = a2longpOrg

a2esurfgp = ma.masked_less_equal(a2esurfgp,0)

#*****************
#- Read MRMS-on-orbit data ----
if (reftype=='mrms')and(os.path.exists(mrmsPath)):
    print 'Read MRMS'
    iymr,eymr = map(int, mrmsPath.split('.')[-2].split('-'))
    a2mrms  = np.load(mrmsPath)
    a2latmr = a2latgpOrg[iymr:eymr+1]
    a2lonmr = a2longpOrg[iymr:eymr+1]

#*****************
#- Read DPR/CMB data ----
#if (reftype=='dpr')and(os.path.exists(mrmsPath)):
if (reftype=='dpr'):
    dprDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch= dprDir + '/2B.GPM.DPRGMI.CORRA*.%06d.V06A.HDF5'%(oid)
    dprPath= glob.glob(ssearch)[0]
    with h5py.File(dprPath,'r') as h:
        a2dpr = h['NS/surfPrecipTotRate'][:]
        a2latdpr = h['NS/Latitude'][:]
        a2londpr = h['NS/Longitude'][:]


#********************************
#-- Draw figure ---
print 'Draw figure'
fig   = plt.figure(figsize=(12,9))
vmin,vmax = 1,20
#vmin,vmax = 0,3   # test

#-- My retrieval --
for i in range(6):
    if (i==5)and(reftype=='dpr'):
        ssize=0.2
    else:
        ssize=1

    if i==0:
        #ax = fig.add_axes([0.1,0.66,0.35,0.3])
        ax = fig.add_axes([0.05,0.1,0.26,0.45])
        a2dat = ma.masked_less_equal(a2MS,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'MS' + ' (>%.2fmm/h)'%(thpr)

    elif i==1:
        #ax = fig.add_axes([0.5,0.66,0.35,0.3])
        ax = fig.add_axes([0.35,0.1,0.26,0.45])
        a2dat = ma.masked_less_equal(a2NS,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'NS' + ' (>%.2fmm/h)'%(thpr)

    elif i==2:
        #ax = fig.add_axes([0.1,0.33,0.35,0.3])
        ax = fig.add_axes([0.05,0.5,0.26,0.45])
        a2dat = ma.masked_less_equal(a2NS,0)
        a2dat = ma.masked_less_equal(a2MScmb,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'MScmb' + ' (>%.2fmm/h)'%(thpr)

    elif i==3:
        #ax = fig.add_axes([0.5,0.33,0.35,0.3])
        ax = fig.add_axes([0.35,0.5,0.26,0.45])
        a2dat = ma.masked_less_equal(a2NScmb,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'NScmb' + ' (>%.2fmm/h)'%(thpr)
    elif i==4:
        #ax = fig.add_axes([0.1,0.0,0.35,0.3])
        ax = fig.add_axes([0.65,0.1,0.26,0.45])
        a2dat = ma.masked_less_equal(a2esurfgp,0)
        a2lat = a2latgp
        a2lon = a2longp
        stype = 'GPROF'
    elif i==5:
        if reftype =='mrms':
            if not os.path.exists(mrmsPath):
                continue
            a2dat = ma.masked_less_equal(a2mrms,0)
            a2lat = a2latmr
            a2lon = a2lonmr
        elif reftype == 'dpr':
            a2dat = ma.masked_less_equal(a2dpr,0)
            a2lat = a2latdpr
            a2lon = a2londpr
        else:
            print 'check reftype',reftype
        #ax = fig.add_axes([0.5,0.0,0.35,0.3])
        ax = fig.add_axes([0.65,0.5,0.26,0.45])
        if reftype=='mrms':
            stype = 'MRMS'
        elif reftype=='dpr':
            stype = 'CMB'

        
    M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=vmin, vmax=vmax)
    M.drawcoastlines()
    
    dgrid      = 5
    parallels  = arange(-90,90, dgrid)
    meridians  = arange(-180,180,dgrid)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')

    plt.title(stype)
#-- Colorbar (Shared) ----
cax = fig.add_axes([0.93, 0.2, 0.02, 0.6])
plt.colorbar(im, orientation='vertical', cax=cax)

#-- Suptitle -------------
ssuptitle = '%04d/%02d/%02d #%06d (%s)'%(Year,Mon,Day,oid,expr)
plt.suptitle(ssuptitle, fontsize=12)
##------------
util.mk_dir(figDir)
outPath  = figDir + '/prcp.map.%s.%06d.y%d-%d.png'%(expr,oid,iy,ey)
plt.savefig(outPath)
print outPath
