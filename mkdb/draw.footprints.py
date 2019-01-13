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
dpr = l2_dpr.L2_DPR()
mwscan= 'S1'
radar = 'Ku'

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
cx  = 110  # GMI center angle bin (py-idx)
cw  = 3    # look at this width at center
w   = int(cw/2)


verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

verDPR = '06'
subverDPR = 'A'
fullverDPR = '%s%s'%(verDPR,subverDPR)

baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
baseDirDPR = '/work/hk01/PMM/NASA/GPM.Ku/2A/V%s'%(verDPR)
idxbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.IDX'%(fullverGMI, mwscan, ix0, ex0, radar, fullverDPR)


srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
#ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.%s.*.HDF5'%(oid)
lsrcPathGMI = glob.glob(ssearchGMI)
srcPathGMI  = lsrcPathGMI[0]

srcDirDPR   = baseDirDPR + '/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch     = srcDirDPR  + '/2A.GPM.%s.*.%s.V%s.HDF5'%(radar,oid,fullverDPR)
lsrcPathDPR = glob.glob(ssearch)
srcPathDPR  = lsrcPathDPR[0]

idxDir      = idxbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
irank  = 1
idxPathX    = idxDir + '/Xpy.%d.%s.npy'%(irank,oid)
idxPathY    = idxDir + '/Ypy.%d.%s.npy'%(irank,oid)


Lat0 = gmi.load_var_granule(srcPathGMI, '%s/Latitude'%(mwscan))
Lon0 = gmi.load_var_granule(srcPathGMI, '%s/Longitude'%(mwscan))
#dtime0= gmi.load_dtime_granule(srcPathGMI, mwscan)

Lat1 = dpr.load_var_granule(srcPathDPR, 'NS/Latitude')
Lon1 = dpr.load_var_granule(srcPathDPR, 'NS/Longitude')
#dtime1= dpr.load_dtime_granule(srcPathDPR, 'NS')

print srcPathDPR
scori = dpr.load_var_granule(srcPathDPR, 'NS/scanStatus/SCorientation')

print scori
print scori.min(),scori.max()

X    = np.load(idxPathX)
Y    = np.load(idxPathY)



#for itmp,y in enumerate([100,500,750,1000,1500,2000,2500]):
#for itmp,y in enumerate([750,1000,1250,1500,1750,2000]):
for itmp,y in enumerate([700]
    #iy,ey = 1000,1001
    iy,ey = y,y+1
    
    ax   = X[iy:ey+1,27-5:27+5+1]
    ay   = Y[iy:ey+1,27-5:27+5+1]
    
    a2latgmi = Lat0[iy:ey+1,cx-15:cx+15+1]
    a2longmi = Lon0[iy:ey+1,cx-15:cx+15+1]
    
    a2latdpr = Lat1[ay.min()-2:ay.max()+2+1, ax.min()-2:ax.max()+2+1]
    a2londpr = Lon1[ay.min()-2:ay.max()+2+1, ax.min()-2:ax.max()+2+1]
    
    #-- centers ----
    a2latgmiC = Lat0[iy:ey+1,110]
    a2longmiC = Lon0[iy:ey+1,110]
    a2latdprC = Lat1[ay.min()-2:ay.max()+2+1, 24]
    a2londprC = Lon1[ay.min()-2:ay.max()+2+1, 24]
   
   
    lllat = a2latdpr.min()-0.1
    urlat = a2latdpr.max()+0.1
    lllon = a2londpr.min()-0.1
    urlon = a2londpr.max()+0.1


    fig   = plt.figure()
    axmap = fig.add_subplot(1,1,1)
    M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
    
    M.scatter(a2londpr, a2latdpr, marker='o',edgecolor='red',facecolor='none',s=100)
    M.scatter(a2longmi, a2latgmi, marker='o',edgecolor='k',facecolor='none',s=100)

    M.scatter(a2londprC, a2latdprC, marker='o',color='red',s=100, alpha=0.7)
    M.scatter(a2longmiC, a2latgmiC, marker='o',color='k', s=100, alpha=0.7)
    
    
    meridians = arange(int(lllon)-3,int(urlon)+3,0.1)
    parallels = arange(int(lllat)-3,int(urlat)+3,0.1)

    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.7)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.7)
 
    ori = (scori.min() + scori.max())*0.5 
    stitle = '%04d/%02d/%02d id=%s SCori=%d scan=%d'%(Year,Mon,Day,oid,ori,y) 
    plt.title( stitle )
    figPath = '/home/utsumi/temp/GMI.DPR.centers.%d.png'%(y)
    plt.savefig(figPath)
    print figPath
