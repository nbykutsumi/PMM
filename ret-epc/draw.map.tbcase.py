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
    mrmsDir  = '/work/hk01/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    workbaseDir = '/home/utsumi/mnt/lab_work'
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    mrmsDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()


## Africa case
#oid = 2421
#iy, ey = 2029, 2089
#clat    = 14 # Africa case
#clon    = 2  # -180 - +180


## West of Great Lakes (good)
#oid = 1748
#Year,Mon,Day = 2014,6,20
#iy, ey = -9999, -9999
#clat    = 42
#clon    = -92
#DB_MAXREC = 20000

# West of Great Lakes (Under est.)
oid = 1574
Year,Mon,Day = 2014,6,8
iy, ey = -9999, -9999
clat    = 41
clon    = -(360-254)
DB_MAXREC = 20000

## Southeast of Grate lakes (Fair)
#oid = 1573
#Year,Mon,Day = 2014,6,8
#iy, ey = -9999, -9999
#clat    = 41
#clon    = -(360-283)
#DB_MAXREC = 20000



## Great Lake Snow oid=003556, 2014/11/20
#oid = 4140
#Year,Mon,Day = 2014,11,20
#iy, ey = -9999, -9999
#clat    = 43
#clon    = -79.5
#DB_MAXREC = 20000

## Great Lake Snow oid=003556, 2015/1/9
#oid = 4914
#Year,Mon,Day = 2015,1,9
#iy, ey = -9999, -9999
#clat    = 43
#clon    = -79.5
#DB_MAXREC = 20000

## SE.US case, oid=003556, 2014/10/14
#oid = 3556
#Year,Mon,Day = 2014,10,14
#iy, ey = 927, 1107
##iy, ey = 1012, 1022
#clat    = 34    # SE.US case. oid = 003556
#clon    = -86   # 2014/10/14  05:42:03 UTC
#DB_MAXREC = 20000

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

#thpr = 0.1  # Minimum threshold for EPC
#thpr = 0.01  # Minimum threshold for EPC
thpr = 0.  # Minimum threshold for EPC

if clat !=-9999.:
    #dlatlon = 20
    dlatlon = 5
    #dscan   = 55
    dscan   = 200  # used to extract GPROf
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
    #BBox    = [[clat-2, clon-4],[clat+2,clon+4]] # test
else:
    BBox    = [[-60,-180],[60,180]]

[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.


#-- Read Tb data --
tbDir   = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)

ssearch = tbDir + '/1C.GPM.GMI.XCAL*.%06d.V05A.HDF5'%(oid)
tbPath  = glob.glob(ssearch)[0]


with h5py.File(tbPath, 'r') as h:
    a3tb1    = h['/S1/Tc'][:]
    a3tb2org = h['/S2/Tc'][:]
    a2lat    = h['/S1/Latitude'][:]
    a2lon    = h['/S1/Longitude'][:]
    a2lat2= h['/S2/Latitude'][:]
    a2lon2= h['/S2/Longitude'][:]

#-- Matchup and Joint S1 and S2 Tb --
s2xPath= tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid)
s2yPath= tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid)

a1x2  = np.load(s2xPath).flatten()
a1y2  = np.load(s2yPath).flatten()

a1mask= ma.masked_less(a1x2,0).mask
a1x2  = ma.masked_less(a1x2,0).filled(0)
a1y2  = ma.masked_less(a1y2,0).filled(0)

nytmp, nxtmp, ztmp = a3tb2org.shape
a2tb2 = a3tb2org[a1y2, a1x2]
a2tb2[a1mask] = miss
a3tb2 = a2tb2.reshape(nytmp,nxtmp,-1)
a3tb = concatenate([a3tb1, a3tb2],axis=2)

print a3tb.shape
##*****************
#
#if clat != -9999.:
#    a2esurfgp, a2latgp, a2longp, iygp, eygp = epcfunc.extract_domain_2D(a2esurfgp, a2latgpOrg, a2longpOrg, clat, clon, dlatlon, dscan, returnidx=True)
#
#else:
#    a2latgp = a2latgpOrg
#    a2longp = a2longpOrg
#
#a2esurfgp = ma.masked_less_equal(a2esurfgp,0)



#********************************
#-- Draw figure ---
print 'Draw figure'
ssize = 1
fig = plt.figure(figsize=(16,12))
for i in range(13):
    ltblabel = [
     r'10.7V'    
    ,r'10.7H'    
    ,r'18.7V'    
    ,r'18.7H'    
    ,r'23.8V'    
    ,r'36.5V'    
    ,r'36.5H'    
    ,r'89.0V'    
    ,r'89.0H'    
    ,r'166V'     
    ,r'166H'     
    ,r'183.31$\pm$3V'
    ,r'183.31$\pm$7V'
    ]

    a2dat = a3tb[:,:,i]

    if int(i/5)==0:
        y0 = 0.05 + 2*0.3
    elif int(i/5)==1:
        y0 = 0.05 + 1*0.3
    elif int(i/5)==2:
        y0 = 0.05 + 0*0.3

    x0 = 0.05 + 0.18*(i%5)
    h  = 0.26
    w  = 0.16
    ax  = fig.add_axes([x0,y0,w,h])

    print a2lat.shape, a2lon.shape, a2dat.shape
    vmin, vmax = 220,290 
    M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=vmin, vmax=vmax)
    M.drawcoastlines()
    
    dgrid      = 5
    parallels  = arange(-90,90, dgrid)
    meridians  = arange(-180,180,dgrid)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')

    stitle = ltblabel[i]
    plt.title(stitle)
#-- Colorbar (Shared) ----
cax = fig.add_axes([0.94, 0.2, 0.015, 0.6])
plt.colorbar(im, orientation='vertical', cax=cax)

#-- Suptitle -------------
ssuptitle = '%04d/%02d/%02d #%06d'%(Year,Mon,Day,oid)
plt.suptitle(ssuptitle, fontsize=12)
##------------
util.mk_dir(figDir)
outPath  = figDir + '/tb.map.%06d.y%d-%d.png'%(oid,iy,ey)
plt.savefig(outPath)
print outPath

##-- My retrieval --
#for i in range(6):
#    if i==0:
#        #ax = fig.add_axes([0.1,0.66,0.35,0.3])
#        ax = fig.add_axes([0.05,0.1,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2MS,0)
#        a2lat = a2latMy
#        a2lon = a2lonMy
#        stype = 'MS' + ' (>%.2fmm/h)'%(thpr)
#
#    elif i==1:
#        #ax = fig.add_axes([0.5,0.66,0.35,0.3])
#        ax = fig.add_axes([0.35,0.1,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2NS,0)
#        a2lat = a2latMy
#        a2lon = a2lonMy
#        stype = 'NS' + ' (>%.2fmm/h)'%(thpr)
#
#    elif i==2:
#        #ax = fig.add_axes([0.1,0.33,0.35,0.3])
#        ax = fig.add_axes([0.05,0.5,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2NS,0)
#        a2dat = ma.masked_less_equal(a2MScmb,0)
#        a2lat = a2latMy
#        a2lon = a2lonMy
#        stype = 'MScmb' + ' (>%.2fmm/h)'%(thpr)
#
#    elif i==3:
#        #ax = fig.add_axes([0.5,0.33,0.35,0.3])
#        ax = fig.add_axes([0.35,0.5,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2NScmb,0)
#        a2lat = a2latMy
#        a2lon = a2lonMy
#        stype = 'NScmb' + ' (>%.2fmm/h)'%(thpr)
#    elif i==4:
#        #ax = fig.add_axes([0.1,0.0,0.35,0.3])
#        ax = fig.add_axes([0.65,0.1,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2esurfgp,0)
#        a2lat = a2latgp
#        a2lon = a2longp
#        stype = 'GPROF'
#    elif i==5:
#
#        if not os.path.exists(mrmsPath):
#            continue
#        #ax = fig.add_axes([0.5,0.0,0.35,0.3])
#        ax = fig.add_axes([0.65,0.5,0.26,0.45])
#        a2dat = ma.masked_less_equal(a2mrms,0)
#        stype = 'MRMS'
#        a2lat = a2latmr
#        a2lon = a2lonmr
#
#    M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
#    im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=vmin, vmax=vmax)
#    M.drawcoastlines()
#    
#    dgrid      = 5
#    parallels  = arange(-90,90, dgrid)
#    meridians  = arange(-180,180,dgrid)
#    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
#    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
#
#    plt.title(stype)
##-- Colorbar (Shared) ----
#cax = fig.add_axes([0.93, 0.2, 0.02, 0.6])
#plt.colorbar(im, orientation='vertical', cax=cax)
#
##-- Suptitle -------------
#ssuptitle = '%04d/%02d/%02d id=%06d'%(Year,Mon,Day,oid)
#plt.suptitle(ssuptitle, fontsize=12)
###------------
#util.mk_dir(figDir)
#outPath  = figDir + '/prcp.map.%06d.y%d-%d.png'%(oid,iy,ey)
#plt.savefig(outPath)
#print outPath
