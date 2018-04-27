import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
from numpy import *
from mpl_toolkits.basemap import Basemap
import hcell_func


lregion = ["NAT","NAF","ASI","NPA","NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]

figDir  = "/tank/utsumi/PMM/HCELL/pict"

figmap = plt.figure()
axmap  = figmap.add_axes([0.1,0.1,0.8,0.8])
M = Basemap( resolution="l", llcrnrlat=-37, llcrnrlon=0, urcrnrlat=37, urcrnrlon=360, ax=axmap)


for region in lregion:
    BBox = hcell_func.ret_regionBBox(region)
    [[lllat,lllon],[urlat,urlon]]=BBox
    if lllon <0: lllon = 360+lllon
    if urlon <0: urlon = 360+urlon

    if lllon < urlon:
        M.plot([lllon,urlon],[lllat,lllat],color="r")
        M.plot([lllon,urlon],[urlat,urlat],color="r")
    else:
        M.plot([lllon,360],[lllat,lllat],color="r")
        M.plot([0,urlon],[lllat,lllat],color="r")
        M.plot([lllon,360],[urlat,urlat],color="r")
        M.plot([0,urlon],[urlat,urlat],color="r")

    M.plot([lllon,lllon],[lllat,urlat],color="r")
    M.plot([urlon,urlon],[lllat,urlat],color="r")

    M.drawcoastlines()
    #M.drawparallels(arange(-30,30+1,30))
    #M.drawmeridians(arange(0,330+1,30))


figPath = figDir + "/region.png"
plt.savefig(figPath)
print figPath
