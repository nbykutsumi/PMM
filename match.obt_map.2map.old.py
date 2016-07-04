from numpy import *
from pyhdf import SD
import os, calendar
from datetime import *
from collections import deque
from glob import glob
from pmm_fsub import *
import pmm_para
import pmm_func
import tag_func
import ctrack_func
#***************************************
#------------------
def ret_a2pr(prtype, year, mon, day, hour):
  if prtype in ["GSMaP.v5","GSMaP.gauge"]:
    a2pr  = pmm_func.ret_a2backward_org(year,mon,day,hour,minute=0, prtype=prtype, maskmiss=False)
  #---------
  return a2pr
#------------------
def global2prmap(a2global, maptype):
  if maptype in ["GSMaP.v5", "GSMaP.gauge"]:
    a2prmap  = a2global[300:1500]    # 60S-60N
  return a2prmap
#------------------
def Load_file(spath,gzflag=True):

  #print "Load_file", spath.split("/")[-1],spath.split("/")[-1][8:8+2],spath.split("/")[-1][10:10+2],spath.split("/")[-1][12:12+5]
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
def ret_prjinfo(obttype):
  if   obttype == "PR.2A25":
    hdftype   = "HDF4"
    vname     = "e_SurfRain"
    Grp       = False
  elif obttype == "TMI.2A12":
    hdftype   = "HDF4"
    vname     = "surfacePrecipitation"
    Grp       = False

  elif obttype == "DPR.L2.DPR":
    hdftype   = "HDF5"
    vname     = "SLV/precipRateESurface"
    Grp       = "NS"

  elif obttype == "KaPR.L2":
    hdftype   = "HDF5"
    vname     = "SLV/precipRateESurface"
    Grp       = "MS"

  elif prj == "KuPR.L2":
    hdftype   = "HDF5"
    vname     = "SLV/precipRateESurface"
    Grp       = "NS"

  elif obttype == "GMI.L2":
    hdftype   = "HDF5"
    vname     = "surfacePrecipitation"
    Grp       = "S1"
  return hdftype, Grp, vname
#------------------
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

#------------------
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



#------------------
def ret_lpath_day(year,mon,day):
  #print "ret_lpath_day"
  if gzflag == True:
    stail = "*.gz"
  else:
    stail = "*01"
  #------
  if obttype == "PR.2A25":
    #idir   = "/mnt/iis.data2/GPM/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
    idir   = "/media/disk2/data/TRMM.PR/L2A25/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"T1PR%04d%02d%02d%s"%(year,mon,day,stail))))
    #print idir
    #print lpath
  elif obttype == "TMI.2A12":
    idir   = "/media/disk2/data/TRMM.TMI/L2A12/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"T1TMI%04d%02d%02d%s"%(year,mon,day,stail))))
  elif obttype == "DPR.L2.DPR":
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/2014/03"
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_DPR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))
  elif obttype == "KaPR.L2":
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/2014/03"
    idir   = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(year,mon)
    lpath  = sorted(glob(os.path.join(idir,"GPMCOR_KAR_%s%02d%02d*.h5"%( ("%04d"%(year))[2:], mon, day ))))



  return lpath

#------------------------------
def get_dat(spath):
  #print "get_dat:",spath
  #---- HDF4 ----------
  if hdftype == "HDF4":
    #h4       = SD.SD(spath)
    h4       = Load_file(spath, gzflag)
    adat        = h4.select(vname)[:]
    alat        = h4.select("Latitude")[:]
    alon        = h4.select("Longitude")[:]
    h4.end()

  #---- HDF5 ----------
  elif hdftype == "HDF5":
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

  #print "lpath="
  #print lpath
  #print "*****************************"

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
#------------------------------
def get_dat(spath):
  #print "get_dat:",spath
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

#------------------------------
def ret_a2prmap_forward(maptype,Y,M,D,H):
  ttime   = datetime(Y,M,D,H) + timedelta(hours=1)
  tY,tM,tD,tH = ttime.year, ttime.month, ttime.day, ttime.hour
  a2prmap = pmm_func.ret_a2backward_org(Y,M,D,H,minute=0,prtype=maptype,maskmiss=True).filled(miss)
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
iY,eY = 2002,2010
iM,eM = 1,12    # test
lY    = range(iY,eY+1)
lM    = range(iM,eM+1)
iD    = 1      # test
#-- parameters -------------
obttype   = "PR.2A25"

para_tag  = pmm_para.Para_tag()
maptype = "GSMaP.v5"
para_tag.prtype = maptype
para_tag.wnflag = False

ny,nx      = pmm_para.ret_nynx(maptype)
#-- tag !! Should be the same as ltag in pmm_func--
if para_tag.wnflag == False:
  #----------
  ltag      = ["tc","c","fbc","te","tf","ef","tef","ot"]
  ltag_elem = ["tc","c","fbc","te","tf","ef","tef"]  
else:
  print "check wnflag!!  wnflag=",wnflag
  sys.exit()

#----------------------------
a2one   = ones([ny,nx],float32)
a2zero  = zeros([ny,nx],float32)

lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(maptype)
gzflag  = False
miss    = -9999.
hdftype, Grp, vname = ret_prjinfo(obttype)

for Y,M in [[Y,M] for Y in lY for M in lM]:
  eD = calendar.monthrange(Y,M)[1]
#  eD = iD     # test
  #-- Initialize ---------------
  da2sum_obt   = {}
  da2sum_map   = {}
  da2num   = {}
  for tag in ltag+["plain"]:
    da2sum_obt[tag] = zeros([ny,nx],float32)
    da2sum_map[tag] = zeros([ny,nx],float32)
    da2num    [tag] = zeros([ny,nx],float32) 
  
  da2pr_tmp  = {}
  for dattype in ["obt","map"]:
    da2pr_tmp[dattype]   = zeros([ny,nx],float32)

  a2num_tmp =  zeros([ny,nx],float32)

  #-----------------------------
  for D in range(iD,eD+1):
    aDat,aLat,aLon,aDtime  = get_dat_days(Y,M,D,Y,M,D)

    idtime = datetime(Y,M,D,0)
    edtime = datetime(Y,M,D,0) + timedelta(days=1)
    ddtime = edtime - idtime
    dH     = int(ddtime.total_seconds() / 3600.)
    ldtime = [idtime + timedelta(hours=H) for H in range(dH)]
    for dtime in ldtime:
#    for dtime in ldtime[:3]:   # test
      #print "***********************"
      #print dtime
      #print "***********************"
      Y = dtime.year
      M = dtime.month
      D = dtime.day
      H = dtime.hour
      ii   = searchsorted(aDtime, dtime)
      ei   = searchsorted(aDtime, dtime+timedelta(hours=1))

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
        lout   = pmm_fsub.obt2dec(adat.T, alon.T, alat.T)
        a2sum  = lout[0].T
        a2num  = lout[1].T
        a2pr   = ma.masked_where(a2num==0.0, a2sum)/a2num
        a2pr   = a2pr /60./60. # mm/h-->mm/s
        a2pr   = a2pr.filled(miss)
        a2pr   = global2prmap(a2pr, maptype)
      #------------------------
      da2pr_tmp["obt"] = a2pr
      a2num_tmp = ma.masked_where(a2pr ==miss, a2one).filled(0.0)  # <-- shared by "obt" and "map"
      #---------------------
      # load map-style precipitation
      #---------------------
      da2pr_tmp["map"]  = ret_a2pr(maptype, Y, M, D, H)

      #--------------------
      # mask out orb where map==miss
      #--------------------
      da2pr_tmp["obt"] = ma.masked_where(da2pr_tmp["map"]==miss, da2pr_tmp["obt"]).filled(miss)

      #---------------------
      # mask out map where obt==miss
      #---------------------
      da2pr_tmp["map"] = ma.masked_where(da2pr_tmp["obt"] ==miss, da2pr_tmp["map"]).filled(miss)


      #---------------------
      # mask out num where map==miss
      #---------------------
      a2num_tmp        = ma.masked_where(da2pr_tmp["map"] ==miss, a2num_tmp).filled(0.0)
      #---------------------
      # replace miss by zero
      #---------------------
      da2pr_tmp["obt"] = ma.masked_equal(da2pr_tmp["obt"], miss).filled(0.0)
      da2pr_tmp["map"] = ma.masked_equal(da2pr_tmp["map"], miss).filled(0.0)
      #--------------------
      # load tag
      #--------------------
      da2tag   = pmm_func.ret_da2tag(Y,M,D,H,para_tag)
      ##*****************************
      # tag precipitation
      #-----------------------------
      da2sum_obt["plain"] = da2sum_obt["plain"] + da2pr_tmp["obt"]
      da2sum_map["plain"] = da2sum_map["plain"] + da2pr_tmp["map"]
      da2num    ["plain"] = da2num["plain"]     + a2num_tmp    # <-- "a2num_tmp" is same for "map" and "obt"

      for tag in ltag:
        a2tagpr_obt = da2pr_tmp["obt"] *da2tag[tag]
        a2tagpr_map = da2pr_tmp["map"] *da2tag[tag]
        a2tagnum    = a2num_tmp        *da2tag[tag]
        #
        da2sum_obt[tag] = da2sum_obt[tag] + a2tagpr_obt
        da2sum_map[tag] = da2sum_map[tag] + a2tagpr_map
        da2num    [tag] = da2num    [tag] + a2tagnum

  #******************  
  # write data
  #******************  
  para_tag.thpr = 0.
  para_tag.ny= ny
  para_tag.nx= nx

  for tag in ltag+["plain"]:
    sdir_root, sdir, sumname1, sumname2, numname \
        =tag_func.ret_name_tagthpr_match_sumnum_mon\
        (\
           maptype, obttype, tag, para_tag, Y, M\
        )

    ctrack_func.mk_dir(sdir)    
    da2sum_map[tag].tofile(sumname1)
    da2sum_obt[tag].tofile(sumname2)
    da2num    [tag].tofile(numname)
    print numname

