from numpy import *
from datetime import *
from collections import deque
import os, gzip, calendar
import h5py
from pyhdf import SD
from glob import glob
from pmm_fsub import *
import ctrack_func
#------------------
# forward precip
#------------------
def ret_prjinfo(prjname):
  if   prjname == "PR.2A25":
    odir_root = "/media/disk2/data/TRMM/PR.2A25.sa.one.FWD"
    hdftype   = "HDF4" 
    vname     = "e_SurfRain"
    Grp       = False
  elif prjname == "TMI.2A12":
    odir_root = "/media/disk2/data/TRMM/TMI.2A12.sa.one.FWD"
    hdftype   = "HDF4" 
    vname     = "surfacePrecipitation"
    Grp       = False
  
  elif prjname == "DPR.L2.DPR":
    odir_root = "/media/disk2/data/DPR/DPR.L2.DPR.sa.one.FWD"
    hdftype   = "HDF5" 
    vname     = "SLV/precipRateESurface"
    Grp       = "NS"
  
  elif prjname == "KaPR.L2":
    odir_root = "/media/disk2/data/DPR/KaPr.L2.sa.one.FWD"
    hdftype   = "HDF5" 
    vname     = "SLV/precipRateESurface"
    Grp       = "MS"
  
  elif prj == "KuPR.L2":
    odir_root = "/media/disk2/data/DPR/KuPr.L2.sa.one.FWD"
    hdftype   = "HDF5" 
    vname     = "SLV/precipRateESurface"
    Grp       = "NS"
  
  elif prjname == "GMI.L2":
    odir_root = "/media/disk2/data/DPR/GMI.L2.sa.one.FWD"
    hdftype   = "HDF5" 
    vname     = "surfacePrecipitation"
    Grp       = "S1"
  return odir_root, hdftype, Grp, vname
#------------------
def Load_file(spath,gzflag=True):

  print "Load_file", spath.split("/")[-1],spath.split("/")[-1][8:8+2],spath.split("/")[-1][10:10+2],spath.split("/")[-1][12:12+5]
  if gzflag == False:
    h4   = SD.SD(spath)

  else:
    f = gzip.open(spath, "rb")
    dat  = f.read()
    f.close()
  
    fname= spath.split("/")[-1][:-3]
  
    year = int(spath.split("/")[-3])
    mon  = int(spath.split("/")[-2])
    odir = odir_root + "/%04d/%02d"%(year,mon)
    tpath= odir + "/%s"%(fname)
    ctrack_func.mk_dir(odir)
    f    = open(tpath, "wb")
    f.write(dat)
    f.close()
  
    h4   = SD.SD(tpath)
    os.remove(tpath)
  return h4
#------------------
def ret_lpath_day(year,mon,day):
  #print "ret_lpath_day"
  if gzflag == True:
    stail = "*.gz"
  else:
    stail = "*01"
  #------
  if prjname == "PR.2A25":
    #idir   = "/mnt/iis.data2/GPM/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
    idir   = "/media/disk2/data/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"T1PR%04d%02d%02d%s"%(year,mon,day,stail))))
  elif prjname == "TMI.2A12":
    idir   = "/media/disk2/data/TRMM.TMI/L2A12/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"T1TMI%04d%02d%02d%s"%(year,mon,day,stail))))
  elif prjname == "DPR.L2.DPR":
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/2014/03"
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_DPR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))
  elif prjname == "KaPR.L2":
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/2014/03"
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_KAR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))



  return lpath

def ret_lpath_days_raw(iyear,imon,iday, eyear,emon,eday):
  #print "ret_lpath_days_raw"
  idtime = datetime(iyear,imon,iday)
  edtime = datetime(eyear,emon,eday) + timedelta(days=1)

  ddtime = edtime - idtime
  dD     = ddtime.days
  ldtime = [idtime + timedelta(days=d) for d in range(dD)]
  lpath  = deque([])
  for dtime in ldtime:
    lpath.append(ret_lpath_day(dtime.year,dtime.month,dtime.day))

  lpath  = concatenate(lpath) 
  return lpath

def ret_lpath_days(iyear,imon,iday,eyear,emon,eday):
  #print "ret_lpath_days"
  idtime = datetime(iyear,imon,iday)
  edtime = datetime(eyear,emon,eday) + timedelta(days=1)

  predtime = idtime - timedelta(days=1)
  posdtime = edtime

  lpath     = deque(ret_lpath_days_raw(iyear,imon,iday,eyear,emon,eday))

  lpath_pre = ret_lpath_day(predtime.year, predtime.month, predtime.day)
  lpath_pos = ret_lpath_day(posdtime.year, posdtime.month, posdtime.day)

  #-- initial --   
  for spath in lpath_pre[::-1]:
    #-- HDF4 ----
    if hdftype == "HDF4":
      #h4       = SD.SD(spath)
      print "1st"
      print spath
      h4       = Load_file(spath, gzflag)
  
      eY,eM,eD = h4.select("Year")[-1], h4.select("Month")[-1], h4.select("DayOfMonth")[-1]
      iY,iM,iD = h4.select("Year")[0],  h4.select("Month")[0],  h4.select("DayOfMonth")[0]
      h4.end()
    #-- HDF5 ----
    if hdftype == "HDF5":
      h5   = h5py.File(spath, "r")
      aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
      aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
      aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
      h5.close()

      if len(aY)==0:
        continue

      eY, eM, eD = aY[-1], aM[-1], aD[-1]
      iY, iM, iD = aY[0],  aM[0],  aD[0]


    if datetime(eY,eM,eD) < idtime:
      break
    elif (datetime(iY,iM,iD) <= idtime and idtime <= datetime(eY,eM,eD)):
      lpath.appendleft(spath)
      break

  #-- last --   
  for spath in lpath_pos:
    print "2nd"
    #-- HDF4 ----
    if hdftype == "HDF4":
      #h4       = SD.SD(spath)
      h4       = Load_file(spath, gzflag)
      eY,eM,eD = h4.select("Year")[-1], h4.select("Month")[-1], h4.select("DayOfMonth")[-1]
      iY,iM,iD = h4.select("Year")[0],  h4.select("Month")[0],  h4.select("DayOfMonth")[0]
      h4.end()

    #-- HDF5 ----
    elif hdftype == "HDF5":
      print "XXXXXXXXXXX"
      print spath
      h5   = h5py.File(spath, "r")
      aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
      aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
      aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
      h5.close()

      if len(aY)==0:
        continue

      eY, eM, eD = aY[-1], aM[-1], aD[-1]
      iY, iM, iD = aY[0],  aM[0],  aD[0]

    print "iD,eD",iD,eD
    if edtime < datetime(iY,iM,iD):
      break
    elif (datetime(iY,iM,iD) <= edtime and edtime <= datetime(eY,eM,eD)):
      lpath.append(spath)
      break
  return lpath

def get_dat(spath):
  print "get_dat:",spath
  #---- HDF4 ----------
  if hdftype == "HDF4":
    #h4       = SD.SD(spath)
    h4       = Load_file(spath, gzflag)
    adat        = h4.select(vname)[:]
    alat        = h4.select("Latitude")[:]
    alon        = h4.select("Longitude")[:]
    aY,aM,aD,aH = h4.select("Year")[:], h4.select("Month")[:], h4.select("DayOfMonth")[:], h4.select("Hour")[:]
    h4.end()

  #---- HDF5 ----------
  elif hdftype == "HDF5":
    h5   = h5py.File(spath, "r")
    adat = h5["%s/%s"%(Grp,vname)][:]
    alat = h5["%s/Latitude"%(Grp)][:]
    alon = h5["%s/Longitude"%(Grp)][:]

    aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
    aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
    aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
    aH   = h5['%s/ScanTime/Hour'%(Grp)][:].astype('int')
    h5.close()

  
  adtime      = [datetime(Y,M,D,H) for (Y,M,D,H) in zip(aY,aM,aD,aH)] 
  return adat, alat, alon, adtime

def get_dat_days(iY,iM,iD,eY,eM,eD):

  lpath  = ret_lpath_days(iY,iM,iD, eY,eM,eD)
  lt =  [[fname.split("/")[-1][8:8+2],fname.split("/")[-1][10:10+2],fname.split("/")[-1][12:12+2]] for fname in lpath ]

  aDat   = deque([])
  aLat   = deque([])
  aLon   = deque([])
  aDtime = deque([])

  print "lpath="
  print lpath
  print "*****************************"

  #********** test ************
  #lpath = list(lpath)[0:2]
  #****************************

  for spath in lpath:
    print "PATH", spath
    adat, alat, alon, adtime = get_dat(spath)   
    aDat.append(adat)
    aLat.append(alat)
    aLon.append(alon)
    aDtime.append(adtime) 

  if len(aDat) == 0:
    aDat   = array([])
    aLat   = array([])
    aLon   = array([])
    aDtime = array([])
  else:
    aDat  = concatenate(aDat)
    aLat  = concatenate(aLat)
    aLon  = concatenate(aLon)
    aDtime= concatenate(aDtime)

  return aDat, aLat, aLon, aDtime
  
def obt2saone_days(iyear,imon,iday,eyear,emon,eday):
  miss            = -9999.
  aDat,aLat,aLon,aDtime  = get_dat_days(iyear,imon,iday,eyear,emon,eday)
  
  idtime = datetime(iyear,imon,iday,0)
  edtime = datetime(iyear,emon,eday,0) + timedelta(days=1)
  ddtime = edtime - idtime
  dH     = int(ddtime.total_seconds() / 3600.)
  ldtime = [idtime + timedelta(hours=H) for H in range(dH)]
  
  for dtime in ldtime:
    ii     = searchsorted(aDtime, dtime)
    ei     = searchsorted(aDtime, dtime+timedelta(hours=1))
  
    adat   = aDat[ii:ei]
    alat   = aLat[ii:ei]
    alon   = aLon[ii:ei]
    adtime = aDtime[ii:ei]
    if len(adtime) ==0:
      print "***************************************"
      print "No data for",dtime
      print "Make Dummy"
      print "***************************************"
      a2pr   = ones([ny,nx],float32)*(-9999.)
    else:
      #lout   = pmm_fsub.obt2saone(adat, alon, alat)
      lout   = pmm_fsub.obt2saone(adat.T, alon.T, alat.T)
      a2sum  = lout[0].T
      a2num  = lout[1].T
      a2pr   = ma.masked_where(a2num==0.0, a2sum)/a2num
      a2pr   = a2pr /60./60. # mm/h-->mm/s
      a2pr   = a2pr.filled(miss)
      a2pr   = a2pr[50:130]  # 40S-40N

    #--------------------
    # write file
    #--------------------
    Y = dtime.year
    M = dtime.month
    D = dtime.day
    H = dtime.hour
    #ofname = "PR.%04d.%02d.%02d.%02d.sa.one"%(Y,M,D,H)
    ofname = "%s.%04d.%02d.%02d.%02d.sa.one"%(prjname.split(".")[0],Y,M,D,H)
    ctrack_func.mk_dir( os.path.join(odir_root, "%04d"%(Y), "%02d"%(M)))
    opath  = os.path.join(odir_root,"%04d"%(Y),"%02d"%(M),ofname)
    a2pr.tofile(opath)
    print opath

#*************************************
ny,nx = 80, 360
gzflag= False
iY,eY = 2014,2014
iM,eM = 3,7
lY    = range(iY,eY+1)
lM    = range(iM,eM+1)
incD  = 3
#incD  = 1
#prjname = "TMI.2A12"
#prjname = "PR.2A25"
#prjname = "DPR.L2.DPR"
prjname = "KaPR.L2"
#prjname = "KuPR.L2"
#prjname = "GMI.L2"

odir_root, hdftype, Grp, vname = ret_prjinfo(prjname) 

for Y,M in [[Y,M] for Y in lY for M in lM]:
  ##---------
  #if Y==2008 and M in [1,2,3]:
  #  continue
  ##---------
  eDayOfMon = calendar.monthrange(Y,M)[1]
  lD        = range(1,eDayOfMon+1)
  for iD in lD[::incD]:
    if iD + incD-1 <= eDayOfMon: 
      eD = iD + incD-1
    else:
      eD = eDayOfMon
    ##---
    if Y==2014 and M==3 and iD <10:
      continue 
    #---
    print Y,M,"days=",iD,"to",eD
    obt2saone_days(Y,M,iD,Y,M,eD) 

