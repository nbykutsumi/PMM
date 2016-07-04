from numpy import *
from pyhdf import SD
import gzip, os
from datetime import *
from collections import deque
import ctrack_func
import itertools
#------------------
def Load_gz(spath):

  print "Load_gz", spath.split("/")[-1],spath.split("/")[-1][8:8+2],spath.split("/")[-1][10:10+2],spath.split("/")[-1][12:12+5]
  f = gzip.open(spath, "rb")
  dat  = f.read()
  f.close()

  fname= spath.split("/")[-1][:-3]

  year = int(spath.split("/")[-3])
  mon  = int(spath.split("/")[-2])
  #odir = odir_root + "/%04d/%02d"%(year,mon)
  #tpath= odir + "/%s"%(fname)
  tpath= "/media/disk2/data/temp" + "/%s"%(fname)
  #ctrack_func.mk_dir(odir)
  f    = open(tpath, "wb")
  f.write(dat)
  f.close()

  h4   = SD.SD(tpath)
  os.remove(tpath)
  return h4
#------------------
def get_dat(prjname, spath):
  if prjname in ["PR.2A25"]:
    vname = "e_SurfRain"
  elif prjname in ["TMI.2A12"]:
    vname = "surfacePrecipitation"
  #-
  #h4       = SD.SD(spath)
  h4       = Load_gz(spath)
  adat        = h4.select(vname)[:]
  alat        = h4.select("Latitude")[:]
  alon        = h4.select("Longitude")[:]
  aY,aM,aD,aH,aMin,aSec = h4.select("Year")[:], h4.select("Month")[:], h4.select("DayOfMonth")[:], h4.select("Hour")[:], h4.select("Minute")[:], h4.select("Second")[:]
  adtime      = [datetime(Y,M,D,H,Min,Sec) for (Y,M,D,H,Min,Sec) in zip(aY,aM,aD,aH,aMin,aSec)]
  h4.end()
  latmax = alat.max()

  return adat, alat, alon, adtime
#------------------
def mk_a1dtime(iY,iM,iD,iH,iMin,iSec, eY,eM,eD,eH,eMin,eSec):
  idtime = datetime(iY,iM,iD,iH,iMin,iSec)
  edtime = datetime(eY,eM,eD,eH,eMin,eSec)
  dtime = idtime
  a1dtime = []
  while dtime <= edtime:
    a1dtime.append(dtime)
    dtime = dtime + timedelta(seconds=1)
  return a1dtime
#------------------

name1 = "/media/disk2/data/TRMM.PR/L2A25/2001/02/T1PR2001020118310_2A25F0007.01.gz"
name2 = "/media/disk2/data/TRMM.TMI/L2A12/2001/02/T1TMI2001020118310_2A12F0007.01.gz"
#
a2dat1, a2lat1, a2lon1, a1dtime1 = get_dat("PR.2A25", name1) 
a2dat2, a2lat2, a2lon2, a1dtime2 = get_dat("TMI.2A12",name2)

iY,iM,iD,iH,iMin,iSec= a1dtime2[0].year, a1dtime2[0].month, a1dtime2[0].day, a1dtime2[0].hour, a1dtime2[0].minute, a1dtime2[0].second

eY,eM,eD,eH,eMin,eSec= a1dtime1[-1].year, a1dtime1[-1].month, a1dtime1[-1].day, a1dtime1[-1].hour, a1dtime1[-1].minute, a1dtime1[-1].second

a1dtimeAll = mk_a1dtime(iY,iM,iD,iH,iMin,iSec,eY,eM,eD,eH,eMin,eSec)

a1out1 = array([])
a1out2 = array([])

olat1 = array([])
olat2 = array([])
olon1 = array([])
olon2 = array([])

thDeg  = 0.02**2.0
#thDeg  = 0.2**2.0
#for dtime in a1dtimeAll[1000:]:
for dtime in a1dtimeAll:
  ik1 = searchsorted(a1dtime1, dtime)
  ek1 = searchsorted(a1dtime1, dtime + timedelta(seconds=1))
  ik2 = searchsorted(a1dtime2, dtime - timedelta(seconds=60+4))
  ek2 = searchsorted(a1dtime2, dtime - timedelta(seconds=60-4))

  print a1dtime1[ik1], shape(a1out1), a2lat1[ik1][0]
  for k1 in range(ik1,ek1+1):
    a1lat1,a1lon1 = a2lat1[k1], a2lon1[k1]
    a1dat1        = a2dat1[k1]

    for k2 in range(ik2,ek2+1):
      a1lat2,a1lon2 = a2lat2[k2], a2lon2[k2]
      a1dat2        = a2dat2[k2]

      a1lat2,a1lon2 = a1lat2[50:-50], a1lon2[50:-50]
      a1dat2        = a1dat2[50:-50]
      
      a2prodDat = array(list(itertools.product(a1dat1,a1dat2)))
      a2prodLat = array(list(itertools.product(a1lat1,a1lat2)))
      a2prodLon = array(list(itertools.product(a1lon1,a1lon2)))
      

      a1dDeg    = square(a2prodLat[:,0] - a2prodLat[:,1]) + square(a2prodLon[:,0] - a2prodLon[:,1])

      a1prodDat1 = a2prodDat[:,0]
      a1prodDat2 = a2prodDat[:,1]
      a1Dat1    = ma.masked_where(a1dDeg > thDeg, a1prodDat1).compressed()
      a1Dat2    = ma.masked_where(a1dDeg > thDeg, a1prodDat2).compressed()

      #a1Lat1    = ma.masked_where(a1dDeg > thDeg, a2prodLat[:,0]).compressed()
      #a1Lat2    = ma.masked_where(a1dDeg > thDeg, a2prodLat[:,1]).compressed()
      #a1Lon1    = ma.masked_where(a1dDeg > thDeg, a2prodLon[:,0]).compressed()
      #a1Lon2    = ma.masked_where(a1dDeg > thDeg, a2prodLon[:,1]).compressed()

      if len(a1Dat1) !=0:
        a1mask    = ma.masked_where((a1Dat1==0.0)&(a1Dat2==0.0), a1Dat1).mask
        a1Dat1  = ma.masked_where(a1mask, a1Dat1).compressed()
        a1Dat2  = ma.masked_where(a1mask, a1Dat2).compressed()

        #a1Lat1  = ma.masked_where(a1mask, a1Lat1).compressed()
        #a1Lat2  = ma.masked_where(a1mask, a1Lat2).compressed()
        #a1Lon1  = ma.masked_where(a1mask, a1Lon1).compressed()
        #a1Lon2  = ma.masked_where(a1mask, a1Lon2).compressed()

      a1out1 =r_[a1out1, a1Dat1]
      a1out2 =r_[a1out2, a1Dat2]

      #olat1  =r_[olat1, a1Lat1]
      #olat2  =r_[olat2, a1Lat2]
      #olon1  =r_[olon1, a1Lon1]
      #olon2  =r_[olon2, a1Lon2]

      #print a1Lat1, a1Lat2,"***",a1Lon1, a1Lon2
      #print a1dtime1[k1],a1dtime2[k2],a1dDeg.min()
