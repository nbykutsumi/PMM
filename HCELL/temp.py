import matplotlib
matplotlib.use("Agg")
from numpy import *
import myfunc.util as util
import myfunc.grids as grids
import hcell_func
from myfunc.fig import Fig

#lregion = ["NAT","NAF","ASI"]
#lregion = ["NAT","NAF","ASI","NPA","NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]
#lregion = ["NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]
lregion = ["NPA"]

for region in lregion:

    BBox= hcell_func.ret_regionBBox(region)
    print BBox

    Lat = arange(-37.0+0.25, 37.0-0.25+0.01, 0.5)
    Lon = arange(0+1, 360-1+0.01, 2.0)
    a2region = grids.mk_mask_BBox(Lat,Lon,BBox)
    bnd = [0,1,2]
    figname = "/home/utsumi/temp/temp.%s.png"%(region)
    Fig.DrawMapSimple(a2region, a1lat=Lat, a1lon=Lon, bnd=bnd, figname=figname)
    print figname


