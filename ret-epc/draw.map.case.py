import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob

## Africa case
#oid = 2421
#iy, ey = 2029, 2089
#clat    = 14 # Africa case
#clon    = 2  # -180 - +180

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

# SW.Japan typhoon case, oid=019015, 2017/07/03
oid = 16166
Year,Mon,Day = 2017,1,1
#iy, ey = -9999,-9999
iy, ey = 2246, 2446
clat    = -10    # SE.US case. oid = 003556
clon    = -75   # 2014/10/14  05:42:03 UTC
DB_MAXREC = 20000



## QJRMS case, oid=012149, 2016/4/18
#oid = 12149
#iy, ey = 1038, 1098
#clat    = 32. # QJRMS case, oid=012149
#clon    = -94 # -180 - +180
#DB_MAXREC = 20000

if clat !=-9999.:
    dlatlon = 20
    #dscan   = 55
    dscan   = 200  # used to extract GPROf
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
else:
    BBox    = [[-60,-180],[60,180]]

[[lllat,lllon],[urlat,urlon]] = BBox
miss = -9999.


#srcDir = '/home/utsumi/temp/out'
#srcDir = '/home/utsumi/temp/out/my'
srcDir = '/tank/utsumi/PMM/retepc/glb.wprof/2017/01/01'
#stamp = '%06d.y%04d-%04d.nrec%d'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

nsurfMSPath    = srcDir + '/nsurfMS.%s.npy'%(stamp)
nsurfNSPath    = srcDir + '/nsurfNS.%s.npy'%(stamp)
nsurfMScmbPath = srcDir + '/nsurfMScmb.%s.npy'%(stamp)
nsurfNScmbPath = srcDir + '/nsurfNScmb.%s.npy'%(stamp)
#topnsurfMScmbPath = srcDir + '/top-nsurfMScmb.%s.npy'%(stamp)
#topnsurfNScmbPath = srcDir + '/top-nsurfNScmb.%s.npy'%(stamp)


latPath   = srcDir + '/lat.%s.npy'%(stamp)
lonPath   = srcDir + '/lon.%s.npy'%(stamp)
#
#-- GPROF --
gprofDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid)
gprofPath= glob.glob(ssearch)[0]

#-- MRMS --
mrmsDir  = '/work/hk01/PMM/MRMS/match-GMI-orbit'
ssearch  = mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.%04d%02d%02d.%06d.*.npy'%(Year,Mon,Day,oid)
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
if os.path.exists(mrmsPath):
    print 'Read MRMS'
    iymr,eymr = map(int, mrmsPath.split('.')[-2].split('-'))
    a2mrms  = np.load(mrmsPath)
    a2latmr = a2latgpOrg[iymr:eymr+1]
    a2lonmr = a2longpOrg[iymr:eymr+1]


#a2dpr, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2dpr, a2latjplOrg, a2lonjplOrg, clat, clon, dlatlon, dscan)
#a2dpr = ma.masked_less_equal(a2dpr,0)
#
#a2latjpl = a2latjplOrg[iyjpl:eyjpl+1,:]
#a2lonjpl = a2lonjplOrg[iyjpl:eyjpl+1,:]
##-- Read Tb data --
#with h5py.File(tbPath, 'r') as h:
#    a3tb1    = h['/S1/Tc'][:]
#    a3tb2org = h['/S2/Tc'][:]
#    a2lattb = h['/S1/Latitude'][:]
#    a2lontb = h['/S1/Longitude'][:]
#
#s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Xpy.1.002421.npy'
#s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Ypy.1.002421.npy'
#
#
##-- Matchup and Joint S1 and S2 Tb --
#a1x2  = np.load(s2xPath).flatten()
#a1y2  = np.load(s2yPath).flatten()
#
#a1mask= ma.masked_less(a1x2,0).mask
#a1x2  = ma.masked_less(a1x2,0).filled(0)
#a1y2  = ma.masked_less(a1y2,0).filled(0)
#
#nytmp, nxtmp, ztmp = a3tb2org.shape
#a2tb2 = a3tb2org[a1y2, a1x2]
#a2tb2[a1mask] = miss
#a3tb2 = a2tb2.reshape(nytmp,nxtmp,-1)
#a3tb = concatenate([a3tb1, a3tb2],axis=2)
#a2mask = np.any(ma.masked_outside(a3tb, 50, 350).mask, axis=2)
#
#print a2mask.shape, a2lattb.shape, a2lontb.shape
#a2mask, a2lattmp, a2lontmp = epcfunc.extract_domain_2D(a2mask, a2lattb, a2lontb, clat, clon, dlatlon, dscan)
#
##-- test ---
#plt.imshow(a2mask*1)
#plt.colorbar()
#plt.savefig('/home/utsumi/temp/out/temp.png')
#plt.clf()
#print 'check a2dat'
#print a2mask[50:55,105:110]




#********************************
#-- Draw figure ---
print 'Draw figure'
fig   = plt.figure(figsize=(8,8))
ssize = 1

#-- My retrieval --
#for i in range(6):
for i in range(5):
    if i==0:
        ax = fig.add_axes([0.1,0.66,0.35,0.3])
        a2dat = ma.masked_less_equal(a2MS,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'MS'

    elif i==1:
        ax = fig.add_axes([0.5,0.66,0.35,0.3])
        a2dat = ma.masked_less_equal(a2NS,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'NS'

    elif i==2:
        ax = fig.add_axes([0.1,0.33,0.35,0.3])
        a2dat = ma.masked_less_equal(a2MScmb,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'MScmb'

    elif i==3:
        ax = fig.add_axes([0.5,0.33,0.35,0.3])
        a2dat = ma.masked_less_equal(a2NScmb,0)
        a2lat = a2latMy
        a2lon = a2lonMy
        stype = 'NScmb'
    elif i==4:
        ax = fig.add_axes([0.1,0.0,0.35,0.3])
        a2dat = ma.masked_less_equal(a2esurfgp,0)
        a2lat = a2latgp
        a2lon = a2longp
        stype = 'GPROF'
    elif i==5:
        ax = fig.add_axes([0.5,0.0,0.35,0.3])
        a2dat = ma.masked_less_equal(a2mrms,0)
        stype = 'MRMS'
        a2lat = a2latmr
        a2lon = a2lonmr

    M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    im    = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=ssize, vmin=0, vmax=20)
    M.drawcoastlines()
    
    dgrid      = 5
    parallels  = arange(-90,90, dgrid)
    meridians  = arange(-180,180,dgrid)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
    plt.colorbar(im, orientation='horizontal')
    
    plt.title(stype)
##------------
figDir   = '/home/utsumi/temp/out/my'
outPath  = figDir + '/prcp.map.%06d.y%d-%d.png'%(oid,iy,ey)
plt.savefig(outPath)
print outPath
