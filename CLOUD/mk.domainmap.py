import matplotlib
matplotlib.use("Agg")

from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


BBox    = [[-0.1, 113.875],[52.1, 180.125]]
lllat = BBox[0][0]-10
lllon = BBox[0][1]-15
urlat = BBox[1][0]+20
urlon = BBox[1][1]+80

parallels =range(-90,90+1,30)
meridians =range(0,360+1,30)

#fig  = plt.figure(figsize= (3,1.8))
fig  = plt.figure(figsize= (4.5,2.7))
ax   = fig.add_axes([0.13,0.1,0.86,0.8])
M    = Basemap(resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon)

M.drawcoastlines(linewidth=0.6)
M.drawparallels(parallels, labels=[1,0,0,0], fontsize=13, linewidth=0.7)
M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=13, linewidth=0.7)

y0,x0  = M(BBox[0][0],BBox[0][1])
y1,x1  = M(BBox[1][0],BBox[1][1])

M.plot([x0,x1],[y0,y0],"-",color="k",linewidth=3)
M.plot([x0,x1],[y1,y1],"-",color="k",linewidth=3)
M.plot([x0,x0],[y0,y1],"-",color="k",linewidth=3)
M.plot([x1,x1],[y0,y1],"-",color="k",linewidth=3)

figPath = "/home/utsumi/temp/domainmap.png"
plt.savefig(figPath)
print figPath


