from numpy import *
from pyhdf import SD
import h5py
import os, calendar, sys
from datetime import *
from collections import deque
from glob import glob
from pmm_fsub import *
import pmm_para
import pmm_func
from cf2.io.GPM import GPM
#***************************************
#------------------
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#------------------
def ret_a2pr(prtype, year, mon, day, hour):
  # ret_a2backward_org returns saone like array (from south to north)
  if prtype in ["GSMaP.v5","GSMaP.gauge","RA"]:
    a2pr  = pmm_func.ret_a2backward_org(year,mon,day,hour,minute=0, prtype=prtype, maskmiss=False)
  #---------
  return a2pr
#------------------
def global2prmap(a2global, maptype):
  if maptype in ["GSMaP.v5", "GSMaP.gauge"]:
    a2prmap  = a2global[300:1500]    # 60S-60N
  elif maptype in ["RA"]:
    a2prmap  = a2global
  return a2prmap
##------------------
#def Load_file(spath,gzflag=True):
#
#  #print "Load_file", spath.split("/")[-1],spath.split("/")[-1][8:8+2],spath.split("/")[-1][10:10+2],spath.split("/")[-1][12:12+5]
#  if gzflag == False:
#    h4   = SD.SD(spath)
#
#  else:
#    f = gzip.open(spath, "rb")
#    dat  = f.read()
#    f.close()
#
#    fname= spath.split("/")[-1][:-3]
#
#    year = int(spath.split("/")[-3])
#    mon  = int(spath.split("/")[-2])
#    odir = odir_root + "/%04d/%02d"%(year,mon)
#    tpath= odir + "/%s"%(fname)
#    mk_dir(odir)
#    f    = open(tpath, "wb")
#    f.write(dat)
#    f.close()
#
#    h4   = SD.SD(tpath)
#    os.remove(tpath)
#  return h4
##------------------
#def ret_prjinfo(obttype):
#  if   obttype == "PR.2A25":
#    hdftype   = "HDF4"
#    vname     = "e_SurfRain"
#    Grp       = False
#  elif obttype == "TMI.2A12":
#    hdftype   = "HDF4"
#    vname     = "surfacePrecipitation"
#    Grp       = False
#
#  elif obttype == "DPR.L2.DPR":
#    hdftype   = "HDF5"
#    vname     = "SLV/precipRateESurface"
#    Grp       = "NS"
#
#  elif obttype == "KaPR.L2":
#    hdftype   = "HDF5"
#    vname     = "SLV/precipRateESurface"
#    Grp       = "MS"
#
#  elif obttype == "KuPR.L2":
#    hdftype   = "HDF5"
#    vname     = "SLV/precipRateESurface"
#    Grp       = "NS"
#
#  elif obttype == "GMI.L2":
#    hdftype   = "HDF5"
#    vname     = "surfacePrecipitation"
#    Grp       = "S1"
#  return hdftype, Grp, vname
##------------------
#def ret_lpath_days_raw(iyear,imon,iday, eyear,emon,eday):
#  #print "ret_lpath_days_raw"
#  idtime = datetime(iyear,imon,iday)
#  edtime = datetime(eyear,emon,eday) + timedelta(days=1)
#
#  ddtime = edtime - idtime
#  dD     = ddtime.days
#  ldtime = [idtime + timedelta(days=d) for d in range(dD)]
#  lpath  = deque([])
#  for dtime in ldtime:
#    lpath.append(ret_lpath_day(dtime.year,dtime.month,dtime.day))
#
#  lpath  = concatenate(lpath)
#  return lpath
#
##------------------
#def ret_lpath_days(iyear,imon,iday,eyear,emon,eday):
#  #print "ret_lpath_days"
#  idtime = datetime(iyear,imon,iday)
#  edtime = datetime(eyear,emon,eday) + timedelta(days=1)
#
#  predtime = idtime - timedelta(days=1)
#  posdtime = edtime
#
#  lpath     = deque(ret_lpath_days_raw(iyear,imon,iday,eyear,emon,eday))
#
#  lpath_pre = ret_lpath_day(predtime.year, predtime.month, predtime.day)
#  lpath_pos = ret_lpath_day(posdtime.year, posdtime.month, posdtime.day)
#
#  #-- initial --
#  for spath in lpath_pre[::-1]:
#    #-- HDF4 ----
#    if hdftype == "HDF4":
#      #h4       = SD.SD(spath)
#      print "1st"
#      print spath
#      h4       = Load_file(spath, gzflag)
#
#      eY,eM,eD = h4.select("Year")[-1], h4.select("Month")[-1], h4.select("DayOfMonth")[-1]
#      iY,iM,iD = h4.select("Year")[0],  h4.select("Month")[0],  h4.select("DayOfMonth")[0]
#      h4.end()
#    #-- HDF5 ----
#    if hdftype == "HDF5":
#      h5   = h5py.File(spath, "r")
#      aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
#      aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
#      aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
#      h5.close()
#
#      if len(aY)==0:
#        continue
#
#      eY, eM, eD = aY[-1], aM[-1], aD[-1]
#      iY, iM, iD = aY[0],  aM[0],  aD[0]
#
#
#    if datetime(eY,eM,eD) < idtime:
#      break
#    elif (datetime(iY,iM,iD) <= idtime and idtime <= datetime(eY,eM,eD)):
#      lpath.appendleft(spath)
#      break
#
#  #-- last --
#  for spath in lpath_pos:
#    print "2nd"
#    #-- HDF4 ----
#    if hdftype == "HDF4":
#      #h4       = SD.SD(spath)
#      h4       = Load_file(spath, gzflag)
#      eY,eM,eD = h4.select("Year")[-1], h4.select("Month")[-1], h4.select("DayOfMonth")[-1]
#      iY,iM,iD = h4.select("Year")[0],  h4.select("Month")[0],  h4.select("DayOfMonth")[0]
#      h4.end()
#
#    #-- HDF5 ----
#    elif hdftype == "HDF5":
#      print "XXXXXXXXXXX"
#      print spath
#      h5   = h5py.File(spath, "r")
#      aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
#      aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
#      aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
#      h5.close()
#
#      if len(aY)==0:
#        continue
#
#      eY, eM, eD = aY[-1], aM[-1], aD[-1]
#      iY, iM, iD = aY[0],  aM[0],  aD[0]
#
#    print "iD,eD",iD,eD
#    if edtime < datetime(iY,iM,iD):
#      break
#    elif (datetime(iY,iM,iD) <= edtime and edtime <= datetime(eY,eM,eD)):
#      lpath.append(spath)
#      break
#  return lpath
#
#
#
##------------------
#def ret_lpath_day(year,mon,day):
#  #print "ret_lpath_day"
#  if gzflag == True:
#    stail = "*.gz"
#  else:
#    stail = "*01"
#  #------
#  if obttype == "PR.2A25":
#    #idir   = "/mnt/iis.data2/GPM/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
#    idir   = "/tank/hjkim/GPM/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
#    lpath  = sorted(glob(os.path.join(idir,"T1PR%04d%02d%02d%s"%(year,mon,day,stail))))
#    #print idir
#    #print lpath
#  elif obttype == "TMI.2A12":
#    #idir   = "/media/disk2/data/TRMM.TMI/L2A12/%04d/%02d"%(year,mon)
#    idir   = "/tank/hjkim/GPM/TRMM.TMI/L2A12/%04d/%02d"%(year,mon)
#    lpath  = sorted(glob(os.path.join(idir,"T1TMI%04d%02d%02d%s"%(year,mon,day,stail))))
#  elif obttype == "DPR.L2.DPR":
#    #idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/%04d/%02d"%(year,mon)
#    idir   = "/tank/hjkim/GPM/GPM.DPR/L2.DPR/%04d/%02d"%(year,mon)
#    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_DPR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))
#  elif obttype == "KaPR.L2":
#    #idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(year,mon)
#    idir   = "/tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(year,mon)
#    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_KAR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))
#  elif obttype == "KuPR.L2":
#    #idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(year,mon)
#    idir   = "/tank/hjkim/GPM/GPM.KuPR/L2/02/%04d/%02d"%(year,mon)
#    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_KUR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))
#
#
#
#  return lpath
#
##------------------------------
#def get_dat(spath):
#  #print "get_dat:",spath
#  #---- HDF4 ----------
#  if hdftype == "HDF4":
#    #h4       = SD.SD(spath)
#    h4       = Load_file(spath, gzflag)
#    adat        = h4.select(vname)[:]
#    alat        = h4.select("Latitude")[:]
#    alon        = h4.select("Longitude")[:]
#    h4.end()
#
#  #---- HDF5 ----------
#  elif hdftype == "HDF5":
#    alat = h5["%s/Latitude"%(Grp)][:]
#    alon = h5["%s/Longitude"%(Grp)][:]
#
#    aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
#    aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
#    aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
#    aH   = h5['%s/ScanTime/Hour'%(Grp)][:].astype('int')
#    h5.close()
#
#
#  adtime      = [datetime(Y,M,D,H) for (Y,M,D,H) in zip(aY,aM,aD,aH)]
#  return adat, alat, alon, adtime
#
#def get_dat_days(iY,iM,iD,eY,eM,eD):
#
#  lpath  = ret_lpath_days(iY,iM,iD, eY,eM,eD)
#  lt =  [[fname.split("/")[-1][8:8+2],fname.split("/")[-1][10:10+2],fname.split("/")[-1][12:12+2]] for fname in lpath ]
#
#  aDat   = deque([])
#  aLat   = deque([])
#  aLon   = deque([])
#  aDtime = deque([])
#
#  #print "lpath="
#  #print lpath
#  #print "*****************************"
#
#  for spath in lpath:
#    print "PATH", spath
#    adat, alat, alon, adtime = get_dat(spath)
#    aDat.append(adat)
#    aLat.append(alat)
#    aLon.append(alon)
#    aDtime.append(adtime)
#
#  if len(aDat) == 0:
#    aDat   = array([])
#    aLat   = array([])
#    aLon   = array([])
#    aDtime = array([])
#  else:
#    aDat  = concatenate(aDat)
#    aLat  = concatenate(aLat)
#    aLon  = concatenate(aLon)
#    aDtime= concatenate(aDtime)
#
#  return aDat, aLat, aLon, aDtime
##------------------------------
#def get_dat(spath):
#  #print "get_dat:",spath
#  #---- HDF4 ----------
#  if hdftype == "HDF4":
#    #h4       = SD.SD(spath)
#    h4       = Load_file(spath, gzflag)
#    adat        = h4.select(vname)[:]
#    alat        = h4.select("Latitude")[:]
#    alon        = h4.select("Longitude")[:]
#    #aY,aM,aD,aH = h4.select("Year")[:], h4.select("Month")[:], h4.select("DayOfMonth")[:], h4.select("Hour")[:]
#    aY,aM,aD,aH,aMin = h4.select("Year")[:], h4.select("Month")[:], h4.select("DayOfMonth")[:], h4.select("Hour")[:], h4.select("Minute")
#    h4.end()
#
#  #---- HDF5 ----------
#  elif hdftype == "HDF5":
#    h5   = h5py.File(spath, "r")
#    adat = h5["%s/%s"%(Grp,vname)][:]
#    alat = h5["%s/Latitude"%(Grp)][:]
#    alon = h5["%s/Longitude"%(Grp)][:]
#
#    aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
#    aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
#    aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
#    aH   = h5['%s/ScanTime/Hour'%(Grp)][:].astype('int')
#    aMin = h5['%s/ScanTime/Minute'%(Grp)][:].astype('int')
#    aSec = h5['%s/ScanTime/Second'%(Grp)][:].astype('int')
#    h5.close()
#
#
#  adtime      = [datetime(Y,M,D,H,Min,Sec) for (Y,M,D,H,Min,Sec) in zip(aY,aM,aD,aH,aMin,aSec)]
#  return adat, alat, alon, adtime

#------------------------------
def ret_a2prmap_forward(maptype,Y,M,D,H):
  ttime   = datetime(Y,M,D,H) + timedelta(hours=1)
  tY,tM,tD,tH = ttime.year, ttime.month, ttime.day, ttime.hour
  a2prmap = pmm_func.ret_a2backward_org(tY,tM,tD,tH,minute=0,prtype=maptype,maskmiss=False)
  return a2prmap
#------------------------------
def ret_coord(maptype):
  if maptype == "RA":
    ny,nx     = 2800,3200
    dlat,dlon = 0.01, 0.01
    lat_first = 20.0  +0.005   # grid center
    lon_first = 118.0 +0.005  # grid center
  elif maptype in ["GSMaP.gauge", "GSMaP.v5"]:
    ny,nx     = 1200,3600
    dlat,dlon = 0.1, 0.1
    lat_first = -59.95   # grid center
    lon_first = 0.05  # grid center

  return lat_first, lon_first, dlat, dlon, ny, nx

#------------------------------
#  Main part
#------------------------------
iY,eY = 2014,2014
iM,eM = 3,3    # test
#iM,eM = 5,5    # test
lY    = range(iY,eY+1)
lM    = range(iM,eM+1)
iD    = 1
#iD    = 15      # test
#-- parameters -------------
#obttype = "PR.2A25"
#obttype = "DPR.L2.DPR"
#obttype = "KaPR.L2"
#obttype = "KuPR.L2"
#obttype = "GMI.L2"
#
#maptype = "GSMaP.gauge"
maptype = "RA"

#obttype  = "GPM.KuPR"
#prdLv    = "L2"
#prdVer   = "02"
#varname  = "NS/SLV/precipRateESurface"

#obttype  = "GPM.DPR"
#prdLv    = "L2"
#prdVer   = "03"
#varname  = "NS/SLV/precipRateESurface"
#
#obttype  = "GPM.GMI"
#prdLv    = "L2"
#prdVer   = "02"
#varname  = "S1/surfacePrecipitation"

obttype  = "TRMM.PR"
prdLv    = "L2A25"
prdVer   = "07"
varname  = "e_SurfRain"

#obttype  = "TRMM.TMI"
#prdLv    = "L2A12"
#prdVer   = "07"
#varname  = "surfacePrecipitation"


ny,nx      = pmm_para.ret_nynx(maptype)
#---- load orbit data ---------------
gpm = GPM(obttype, prdLv,prdVer)
BBox= [[20.0, 118.0], [48.0, 150.0]]

#----------------------------
a2one   = ones([ny,nx],float32)
a2zero  = zeros([ny,nx],float32)

lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(maptype)
gzflag  = False
miss    = -9999.
#hdftype, Grp, vname = ret_prjinfo(obttype)

for Y,M in [[Y,M] for Y in lY for M in lM]:
  eD = calendar.monthrange(Y,M)[1]
  #-- Initialize ---------------
  a2sum_obt  = zeros([ny,nx],float32)
  a2sum_map  = zeros([ny,nx],float32)
  a2num      = zeros([ny,nx],float32)
 
  da2pr_tmp  = {}
  for dattype in ["obt","map"]:
    da2pr_tmp[dattype]   = zeros([ny,nx],float32)

  a2num_tmp =  zeros([ny,nx],float32)

  #-----------------------------
#  for D in range(iD,eD+1):
  for D in range(iD,iD+1):   # test
    #aDat,aLat,aLon,aDtime  = get_dat_days(Y,M,D,Y,M,D)

    idtime = datetime(Y,M,D,0)
    edtime = datetime(Y,M,D,0) + timedelta(days=1)
    ddtime = edtime - idtime
    dH     = int(ddtime.total_seconds() / 3600.)
    ldtime = [idtime + timedelta(hours=H) for H in range(dH)]
    for dtime in ldtime:
#    for dtime in ldtime[4:23+1]:   # test
      print "***********************"
      print dtime
      print "***********************"
      Y = dtime.year
      M = dtime.month
      D = dtime.day
      H = dtime.hour
      #ii   = searchsorted(aDtime, dtime)
      #ei   = searchsorted(aDtime, dtime+timedelta(hours=1))

      #adat   = aDat[ii:ei]
      #alat   = aLat[ii:ei]
      #alon   = aLon[ii:ei]
      #adtime = aDtime[ii:ei]


      #if len(adtime) ==0:
      #  print "***************************************"
      #  print "No data for",dtime
      #  print "Make Dummy"
      #  print "***************************************"
      #  a2pr   = ones([ny,nx],float32)*(-9999.)
      #else:
      #  ##--- test ------------
      #  #def lon2x( a2lon):
      #  #  a2out = floor( (a2lon - 118.0)/0.01 )  
      #  #  return a2out
      #  #def lat2x( a2lat):
      #  #  a2out = floor( (a2lat - 20.0)/0.01 )  
      #  #  return a2out
      #  #alat = array([0.1,10.004, 20.005, 20.015, 20.025, 60.005], float32).reshape(2,-1)
      #  #alon = array([0.2,10.004, 118.005, 118.015, 118.025, 160.005], float32).reshape(2,-1)
      #  ##adat = lon2x( alon )
      #  #adat = lat2x( alat )
      #  
      #  #---------------------
      #  if maptype in ["GSMaP.v5","GSMaP.gauge"]:
      #    lout   = pmm_fsub.obt2dec(adat.T, alon.T, alat.T)
      #  elif maptype in ["RA"]:
      #    lout   = pmm_fsub.obt2jp2800x3200(adat.T, alon.T, alat.T)
      #  a2s    = lout[0].T
      #  a2n    = lout[1].T
      #  a2pr   = ma.masked_where(a2n==0.0, a2s)/a2n
      #  a2pr   = a2pr /60./60. # mm/h-->mm/s  # org
      #  #a2pr   = a2pr # mm/h-->mm/s    # test
      #  a2pr   = a2pr.filled(miss)
      #  a2pr   = global2prmap(a2pr, maptype)

      #------------------------

      #----- cf module ----------
      sDTime = dtime
      eDTime = sDTime + timedelta(hours=1)
      #-----------------------------
      print sDTime
      print eDTime
      #-----------------------------


      try:      
        #obt    = gpm(varname, sDTime, eDTime, [[-89.99,-180.0],[89.99,179.99]]) 
        obt    = gpm(varname, sDTime, eDTime, BBox) 
      except ValueError:
        print ""
        print "SKIP!! ValueError", sDTime
        print ""
        continue
      except IndexError :
        print ""
        print "SKIP!! IndexError", sDTime
        continue


      lout2   = pmm_fsub.obt2jp2800x3200(obt.data.T, obt.lon.T, obt.lat.T)
      a2s2    = lout2[0].T
      a2n2    = lout2[1].T
      a2pr    = ma.masked_where(a2n2==0.0, a2s2)/a2n2
      a2pr    = a2pr /60./60. # mm/h-->mm/s  # org
      a2pr    = a2pr.filled(miss)
      #------------------------------
      da2pr_tmp["obt"] = a2pr

      #---------------------
      # load map-style precipitation
      #---------------------
      #da2pr_tmp["map"]  = ret_a2pr(maptype, Y, M, D, H)
      da2pr_tmp["map"]  = ret_a2prmap_forward(maptype, Y, M, D, H)


      if a2pr.max() >0.0:
        print "obt.max", da2pr_tmp["obt"].max()
        print "map.max", da2pr_tmp["map"].max()
        sys.exit()
      #--------------------
      # mask out obt where map==miss
      #--------------------

      da2pr_tmp["obt"] = ma.masked_where(da2pr_tmp["map"]==miss, da2pr_tmp["obt"]).filled(miss)

      #---------------------
      # mask out map where obt==miss
      #---------------------
      da2pr_tmp["map"] = ma.masked_where(da2pr_tmp["obt"] ==miss, da2pr_tmp["map"]).filled(miss)

      #---------------------
      # mask out num where map==miss
      #---------------------
      #a2num_tmp        = ma.masked_where(da2pr_tmp["map"] ==miss, a2num_tmp).filled(0.0)
      a2num_tmp        = ma.masked_where(da2pr_tmp["map"] ==miss, a2one).filled(0.0)
      #---------------------
      # replace miss by zero
      #---------------------
      da2pr_tmp["obt"] = ma.masked_equal(da2pr_tmp["obt"], miss).filled(0.0)
      da2pr_tmp["map"] = ma.masked_equal(da2pr_tmp["map"], miss).filled(0.0)
      #--------------------
      a2sum_obt  = a2sum_obt + da2pr_tmp["obt"]
      a2sum_map  = a2sum_map + da2pr_tmp["map"]
      a2num      = a2num     + a2num_tmp
  #******************  
  # write data
  #******************  
  odir_base = "/tank/utsumi/PMM/COMP.MAP.OBT/%s.vs.%s"%(maptype, obttype)
  odir      = odir_base + "/%04d"%(Y)
  mk_dir(odir)
  print odir

  NY,NX     = shape(da2pr_tmp["map"])
  sumname_map = odir + "/sum.%s.%02d.%d.%d"%(maptype, M, NY, NX) 
  sumname_obt = odir + "/sum.%s.%02d.%d.%d"%(obttype, M, NY, NX) 
  numname     = odir + "/num.%02d.%d.%d"%(M, NY, NX) 

  a2sum_obt   = array(a2sum_obt, float32)
  a2sum_map   = array(a2sum_map, float32)
  a2num_map   = array(a2num, float32)

  a2sum_obt.tofile(sumname_obt) 
  a2sum_map.tofile(sumname_map)
  a2num.tofile(numname)
  print sumname_obt
