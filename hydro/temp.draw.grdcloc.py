import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import myfunc.IO.GRDC as GRDC
import numpy as np
from collections import deque

listDir = '/work/hk01/utsumi/PMM/hydro/list'
listPath= listDir + '/grdc_loc_rev_2000_2009.txt'

f=open(listPath,'r'); lines=f.readlines(); f.close()
lstnid = []
dyx    = {}
dlatlon= {}
llat   = []
llon   = []
for line in lines:
    line  = line.split()
    stnid = int(line[0])
    x       = int(line[6])  # 0,1,2..
    y       = int(line[7])  # 0,1,2..
    lon     = float(line[4])
    lat     = float(line[5])
    bias    = abs(float(line[10]))
    if bias >= 0.05: continue
    dyx[stnid]     = [y,x]
    dlatlon[stnid] = [lat,lon]
    lstnid.append(stnid)
    llat  .append(lat)
    llon  .append(lon)

#***************************************************
# Draw map
#***************************************************
BBox = [[-90,-180],[90,180]]
[[lllat,lllon],[urlat,urlon]] = BBox
fig  = plt.figure(figsize=(8,4))
fig  = plt.figure()
ax   = fig.add_subplot(1,1,1)
M    = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
im  = M.scatter(llon, llat, marker='o', s=60)

meridians = arange(-180,180,30)
parallels = arange(-90,90,30)
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.7)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.7)
M.drawcoastlines()
stitle  = 'GHCN Stations'
plt.show()
#-- save --
figDir  = '/work/hk01/utsumi/PMM/hydro/fig'
figPath = figDir + '/map.ghcn_loc.png'
plt.savefig(figPath)
print figPath
 
