import matplotlib
matplotlib.use("Agg")
from numpy import *
import myfunc.IO.GSMaP as GSMaP
from datetime import datetime, timedelta
import myfunc.fig.Fig as Fig
import calendar
import myfunc.util as util

#calcFlag= False
calcFlag= True
iYM = [2014,1]
eYM = [2014,2]
iYear,iMon = iYM
eYear,eMon = eYM
lYM  = util.ret_lYM(iYM,eYM)
BBox = [[33,360-121],[37,360-115]]
parallels = arange(-90,90+0.1,2)
meridians = arange(-180,360+0.1,2)


gsmap = GSMaP.GSMaP(prj="reanalysis",ver="v6",compressed=False)
Lat    = gsmap.Lat
Lon    = gsmap.Lon
ny     = gsmap.ny
nx     = gsmap.nx
print lYM

vmin = 0
vmax = 100
cmap      = "gnuplot2_r"
#cmap      = "rainbow"

figDir = "/home/utsumi/members/hjkim/PMM/fig"
util.mk_dir(figDir)
cbarPath=figDir + "/cbar.GSMaP.png"
for YM in lYM:
    Year,Mon = YM
    print Year,Mon
    oDir   = "/home/utsumi/members/hjkim/PMM"
    oPath  = oDir + "/GSMaP.%04d.%02d.%dx%d"%(Year,Mon,ny,nx)
    if calcFlag == True:
        iDTime = datetime(Year,Mon,1,0)
        eDay   = calendar.monthrange(Year,Mon)[1]
        eDTime = datetime(eYear,Mon,eDay,23)
        dDTime = timedelta(hours=1)
     
        a2dat  = gsmap.time_sum_mmh(iDTime, eDTime, dDTime)
        a2dat.tofile(oPath)
        print oPath

    a2fig  = fromfile(oPath,float32).reshape(ny,nx) 
    figPath= figDir + "/GSMaP.%04d.%02d.png"%(Year,Mon)
#    Fig.DrawMapSimple(a2in=a2fig, a1lat=Lat, a1lon=Lon,BBox=BBox, figname=figPath,cbarname=cbarPath, cmap=cmap,vmin=vmin, vmax=vmax)

    Fig.DrawMapSimple(a2in=a2fig, a1lat=Lat, a1lon=Lon,BBox=BBox, figname=figPath,cbarname=cbarPath, cmap=cmap, parallels=parallels, meridians=meridians, vmin=vmin, vmax=vmax)

    print figPath

