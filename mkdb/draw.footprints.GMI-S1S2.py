import matplotlib
matplotlib.use('Agg')
from numpy import *
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import myfunc.IO.GPM.l2_dpr as l2_dpr
import myfunc.IO.GPM.l1_gmi as l1_gmi
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import sys, os
import myfunc.util as util
from f_match_fov import *


#DTime = datetime(2017,1,2)
DTime = datetime(2017,1,30)
Year,Mon,Day = DTime.timetuple()[:3]
#oid   = '016171' # 2017-1-2  SCori=180
oid   = '016608'  # 2017-1-30 SCori=0

gmi = l1_gmi.L1_GMI()

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
cx  = 110  # GMI center angle bin (py-idx)
cw  = ex0-ix0+1    # look at this width at center
#cw  = 221    # look at this width at center
w   = int(cw/2)


verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
#ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.%s.*.HDF5'%(oid)
lsrcPathGMI = glob.glob(ssearchGMI)
srcPathGMI  = lsrcPathGMI[0]

Lat1 = gmi.load_var_granule(srcPathGMI, 'S1/Latitude')
Lon1 = gmi.load_var_granule(srcPathGMI, 'S1/Longitude')
Lat2 = gmi.load_var_granule(srcPathGMI, 'S2/Latitude')
Lon2 = gmi.load_var_granule(srcPathGMI, 'S2/Longitude')

dtime1=gmi.load_dtime_granule(srcPathGMI,'S1')
dtime2=gmi.load_dtime_granule(srcPathGMI,'S2')

scori = gmi.load_var_granule(srcPathGMI, 'S1/SCstatus/SCorientation')


#for itmp,y in enumerate([100,500,750,1000,1500,2000,2500]):
#for itmp,y in enumerate([750,1000,1250,1500,1750,2000]):
#for itmp,y in enumerate([0]):
for itmp,y in enumerate([0,750,1000,1250,1500,2000]):
    #iy,ey = 1000,1001
    iy,ey = y,y+2
   

    print 'S1',dtime1[iy:ey+1]
    print 'S2',dtime2[iy:ey+1]
 
    a2lat1 = Lat1[iy:ey+1,cx-w:cx+w+1]
    a2lon1 = Lon1[iy:ey+1,cx-w:cx+w+1]

    a2lat2 = Lat2[iy:ey+1,cx-w:cx+w+1]
    a2lon2 = Lon2[iy:ey+1,cx-w:cx+w+1]


    #-- centers ----
    a2latC1 = Lat1[iy:ey+1,cx]
    a2lonC1 = Lon1[iy:ey+1,cx]

    a2latC2 = Lat2[iy:ey+1,cx]
    a2lonC2 = Lon2[iy:ey+1,cx]
   
    lllat = min(a2lat1.min(),a2lat2.min())-0.1
    urlat = max(a2lat1.max(),a2lat2.max())+0.1
    lllon = min(a2lon1.min(),a2lon2.min())-0.1
    urlon = max(a2lon1.max(),a2lon2.max())+0.1

    #lllat = min(a2lat1[:,:30].min(),a2lat2[:,:30].min())-0.1
    #urlat = max(a2lat1[:,:30].max(),a2lat2[:,:30].max())+0.1
    #lllon = min(a2lon1[:,:30].min(),a2lon2[:,:30].min())-0.1
    #urlon = max(a2lon1[:,:30].max(),a2lon2[:,:30].max())+0.1



    fig   = plt.figure()
    axmap = fig.add_subplot(1,1,1)
    M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
    
    M.scatter(a2lon1, a2lat1, marker='o',edgecolor='k',facecolor='none',s=100)
    M.scatter(a2lon2, a2lat2, marker='o',edgecolor='gray',facecolor='none',s=100)

    M.scatter(a2lonC1, a2latC1, marker='o',color='k',s=100, alpha=0.7)
    M.scatter(a2lonC2, a2latC2, marker='o',color='gray', s=100, alpha=0.7)
    
    #wgrid = 0.1
    wgrid = 1 
    meridians = arange(int(lllon)-3,int(urlon)+3,wgrid)
    parallels = arange(int(lllat)-3,int(urlat)+3,wgrid)

    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.7)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.7)
 
    ori = (scori.min() + scori.max())*0.5 
    stitle = '%04d/%02d/%02d id=%s SCori=%d scan=%d'%(Year,Mon,Day,oid,ori,y) 
    plt.title( stitle )
    figPath = '/home/utsumi/temp/GMI.S1S2.centers.%d.png'%(y)
    plt.savefig(figPath)
    print figPath
