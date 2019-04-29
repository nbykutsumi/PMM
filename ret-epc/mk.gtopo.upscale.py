import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from numpy import *
import numpy as np
import sys, os
import glob
import myfunc.regrid.Regrid as Regrid

procparts = False
#procparts = True
procjoint = True
#-- Read GTOPO ----
lullat = [90,40,-10,-60]
lullon = [-180,-140,-100,-60,-20,20,60,100,140]
#lullat = [90]
#lullon = [-180]

llatlon = [(ullat,ullon) for ullat in lullat for ullon in lullon]

dlatUp  = 0.1
dlonUp  = 0.1
orogDir  = "/data1/hjkim/GTOPO30"
#orogDir  = "/media/disk2/share/data/GTOPO30"
upDir  = '/work/hk01/utsumi/gtopo/0.1deg'

ddem    = {}
for (ullat,ullon) in llatlon:
    if procparts != True: continue

    #ullat = int( (lat - (-60))/50. )*50. + 50 -60.
    #ullon = int( (lon - (-180))/40.)*40. -180.

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<=ullon)&(ullon<0):
        WE = "W"
    else:
        WE = "E"

    #orogPath = orogDir + "/E060N40.DEM"
    orogPath = orogDir + "/%s%03d%s%02d.DEM"%(WE, abs(ullon), SN, abs(ullat))
    dmwPath  = orogDir + "/%s%03d%s%02d.DMW"%(WE, abs(ullon), SN, abs(ullat))

    if not os.path.exists(orogPath):
        continue

    print orogPath

    """
    BYTEORDER      M
    LAYOUT       BIL
    NROWS         6000
    NCOLS         4800
    NBANDS        1
    NBITS         16
    BANDROWBYTES         9600
    TOTALROWBYTES        9600
    BANDGAPBYTES         0
    NODATA        -9999
    ULXMAP        60.00416666666667
    ULYMAP        39.99583333333333
    XDIM          0.00833333333333
    YDIM          0.00833333333333
    """

    #ny = 6000
    nxOrg = 4800
    #ddem[(ullat,ullon)] = ma.masked_equal(flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nxOrg)),-9999).filled(0)
    a2demOrg = ma.masked_equal(flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nxOrg)),-9999).filled(0)
    nyOrg   = a2demOrg.shape[0]

    # load DMW file
    f=open(dmwPath, "r"); lines=f.readlines(); f.close()
    lonmin = float(lines[4])
    latmax = float(lines[5])
    lonmax = lonmin + 0.00833333333333*(nxOrg-1)
    latmin = latmax - 0.00833333333333*(nyOrg-1)

    dlatOrg = 0.00833333333333
    dlonOrg = 0.00833333333333
    a1latOrg= arange(latmin,latmax+0.5*dlatOrg, dlatOrg)
    a1lonOrg= arange(lonmin,lonmax+0.5*dlonOrg, dlonOrg)
   
    # Set upscale properties     
    lllat   = float(round(latmin))
    lllon   = float(round(lonmin))
    urlat   = float(round(latmax))
    urlon   = float(round(lonmax))

    a1latUp = arange(lllat+dlatUp*0.5, urlat-dlatUp*0.5+dlatUp*0.1, dlatUp)
    a1lonUp = arange(lllon+dlonUp*0.5, urlon-dlonUp*0.5+dlonUp*0.1, dlonUp)

    print 'nyUp, nxUp=',len(a1latUp),len(a1lonUp)

    us = Regrid.UpScale()
    us(a1latOrg, a1lonOrg, a1latUp, a1lonUp, globflag=False)
    a2demUp = us.upscale(a2demOrg, pergrid=False, miss_in=-9999, miss_out=-9999.).astype(int16)

    #-- Save ---
    stamp  = 'lat.%03d.%03d.lon.%04d.%04d'%(lllat,urlat,lllon,urlon)
    upPath = upDir + '/dem.%s.npy'%(stamp)
    latPath= upDir + '/lat.%s.npy'%(stamp)
    lonPath= upDir + '/lon.%s.npy'%(stamp)

    np.save(upPath, a2demUp)
    np.save(latPath, a1latUp)
    np.save(lonPath, a1lonUp)
    print upPath

#-- Joint maps ------
lat0 = -90
lon0 = -180
lat1 = 90
lon1 = 360+lon0

a1latGlob = arange(lat0+dlatUp*0.5, 90-dlatUp*0.5, dlatUp)
a1lonGlob = arange(lon0+dlonUp*0.5, 360+lon0-dlonUp*0.5, dlonUp)
ny = int(round((90-lat0)/dlatUp))
nx = int(round(360/dlonUp))
a2glob = zeros([ny,nx],int16)

for (ullat, ullon) in llatlon:
    lsrcPath = glob.glob(upDir + '/dem.lat.???.%03d.lon.%04d.????.npy'%(ullat, ullon))

    if len(lsrcPath)==0: continue

    srcPath = lsrcPath[0]
    lllat,urlat = map(int, srcPath.split('.')[-6:-6+2])
    lllon,urlon = map(int, srcPath.split('.')[-3:-3+2])
    print srcPath

    a2reg = np.load(srcPath)
    nyReg,nxReg = a2reg.shape
    iy = int((lllat - lat0)/dlatUp)
    ix = int((lllon - lon0)/dlonUp)

    print a2glob.shape, a2reg.shape
    a2glob[iy:iy+nyReg, ix:ix+nxReg] = a2reg

outPath = upDir + '/dem.lat.%03d.%03d.lon.%04d.%04d.%dx%d.npy'%(lat0,lat1,lon0,lon1,ny,nx)
np.save(outPath, a2glob)
print outPath

#--- Figure ----
a2glob = ma.masked_equal(a2glob,0)
fig = plt.figure(figsize=(8,8))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
M = Basemap(llcrnrlat=lat0,llcrnrlon=lon0, urcrnrlat=lat1, urcrnrlon=lon1, ax=ax, resolution='l')
im = M.imshow(a2glob, origin='lower', cmap='terrain')
M.drawcoastlines()
plt.colorbar(im, orientation='horizontal')

meridians = arange(-180,180+1,10)
parallels = arange(-90,90+1, 10)
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d', rotation=-90)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')

figPath = upDir + '/dem.lat.%03d.%03d.lon.%04d.%04d.%dx%d.png'%(lat0,lat1,lon0,lon1,ny,nx)
plt.savefig(figPath)
print figPath
plt.clf()
