from numpy import *
from datetime import datetime, timedelta
import ioPMM
import myfunc.util as util
import myfunc.fig.Fig as Fig
iopmm  = ioPMM.ioRAdom()
match  = ioPMM.matchRAdom()
match(prtype1="RA", prtype2="GPM.KuPR", Trc="ALL", nclass=4)

Year = 2014
Mon  = 6
lat = iopmm.Lat
lon = iopmm.Lon
BBox= iopmm.BBox


aall = iopmm.loadMonRAmmh(Year,Mon,var="ave",maskflag=True)
acol = match.ret_monPrc(prtype="RA", tag="plain", Year=Year,Mon=Mon)

sallname = "./temp.Mon.all.png"
scolname = "./temp.Mon.col.png"

stitle   = "%s-%02d"%(Year, Mon)
Fig.DrawMap(a2in=acol, a1lat=lat, a1lon=lon, BBox=BBox, figname=scolname, stitle=stitle, parallels=arange(-90,90+0.01,5), meridians=arange(0.0,360+0.01,5),vmax=4)
Fig.DrawMap(a2in=aall, a1lat=lat, a1lon=lon, BBox=BBox, figname=sallname, stitle=stitle, parallels=arange(-90,90+0.01,5), meridians=arange(0.0,360+0.01,5),vmax=4)
#



#iopmm  = ioPMM.ioRAdom()
#match  = ioPMM.matchRAdom()
#match(prtype1="RA", prtype2="GPM.KuPR", Trc="ALL", ncclass=4)
#iDTime = datetime(2014,6,22,0)
#eDTime = datetime(2014,6,22,0)
#iDTime = datetime(2014,6,1,0)
#eDTime = datetime(2014,6,30,18)
#dDTime = timedelta(hours=6)
#lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
#lat = iopmm.Lat
#lon = iopmm.Lon
#BBox= iopmm.BBox
#for DTime in lDTime:
#  try:
#    amask = iopmm.loadGPMmmh(DTime, prj="GPM.KuPR", maskflag=True)
#  except TabError:
#    continue
#  aall = iopmm.loadRAmmh(DTime, maskflag=True)
#  acol = ma.masked_where(amask.mask, aall)
#  print DTime, aall.mean(), acol.mean()
#
#  sallname = "./temp.all.png"
#  scolname = "./temp.col.png"
#
#  stitle   = "%s"%(DTime)
#  Fig.DrawMap(a2in=acol, a1lat=lat, a1lon=lon, BBox=BBox, figname=scolname, stitle=stitle, parallels=arange(-90,90+0.01,5), meridians=arange(0.0,360+0.01,5),vmax=20)
#  Fig.DrawMap(a2in=aall, a1lat=lat, a1lon=lon, BBox=BBox, figname=sallname, stitle=stitle, parallels=arange(-90,90+0.01,5), meridians=arange(0.0,360+0.01,5),vmax=20)
#

