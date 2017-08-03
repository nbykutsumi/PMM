from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
from collections import deque
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import myfunc.util as util
import calendar
import sys
from mpl_toolkits.basemap import Basemap

#clVer = "JMA1"
#clVer = "MyWNP1"
clVer = "MyWNP3"

iYM    = [2014,4]
#eYM    = [2014,4]
eYM    = [2015,6]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
ldattype = ["KuPR","GMI"]

BBox    = [[-0.1, 113.875],[52.1, 180.125]]
ny,nx   = 261, 265
miss    = -9999.

vmin    = 0.1

#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = int(clVer[5:])
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%d"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%d"%(ver)



lcltype = cl.licl
ncltype = len(lcltype)
dclName = cl.dclName
dclShortName = cl.dclShortName
dclid   = cl.dclid
LatUp   = cl.Lat
LonUp   = cl.Lon

# RadarAMeDAS
ra    = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
LatOrgRA = ra.Lat
LonOrgRA = ra.Lon
usRa  = Regrid.UpScale()
usRa(LatOrgRA, LonOrgRA, LatUp, LonUp, globflag=False)

maskPath  = "/tank/utsumi/data/RadarAMeDAS/mask/RAmask.kubota.0.20x0.25WNP.261x265" 
a2ramask  = fromfile(maskPath, float32).reshape(261,265)
a2ramask  = ma.masked_equal(a2ramask, -9999.)

# Init Precip
#*****************
class MyIOException(Exception): pass


#*****************
def ret_dDTime(dattype):
  if dattype.split(".")[0]=="GSMaP":
    return timedelta(hours=1)
  elif dattype == "RA":
    return timedelta(hours=1)
  elif dattype.split(".")[0]=="IMERG":
    return timedelta(minutes=30)
  elif dattype in ["KuPR","GMI"]:
    return timedelta(minutes=30)
  else:
    print "check ret_dDTime"; sys.exit()

#*****************
def ret_eMinute(dattype):
  if   ret_dDTime(dattype)== timedelta(hours=1):
    eMinute = 0
  elif ret_dDTime(dattype)== timedelta(minutes=30):
    eMinute = 30 
  else:
    raise ValueError
  return eMinute
#*****************
def ret_mmh(DTime, dattype):
  if   dattype =="GSMaP":
    a2prOrg = ma.masked_less(gsmap.load_mmh(DTime),0.0).filled(miss)    # mm/h, forward
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)
  elif dattype =="GSMaP.MW":
    a2sateinfo = gsmap.load_sateinfo(DTime)
    a2prOrg = ma.masked_where(a2sateinfo <0,
               ma.masked_less(gsmap.load_mmh(DTime),0.0)
              ).filled(miss)    # mm/h, forward

    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="GSMaP.IR":
    a2sateinfo = gsmap.load_sateinfo(DTime)
    a2prOrg = ma.masked_where(a2sateinfo >=0,
                ma.masked_less(gsmap.load_mmh(DTime),0.0)
              ).filled(miss)    # mm/h, forward

    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime),0.0).filled(miss)    # mm/h, forward
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG.IR":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime, var="IRprecipitation"),0.0).filled(miss)    # mm/h, forward
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG.MW":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime, var="HQprecipitation"),0.0).filled(miss)    # mm/h, forward
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="RA":
    a2prOrg = ra.loadForward_mmh(DTime,mask=True).filled(miss)    # mm/h, forward
    a2pr    = usRa.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="KuPR":
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour
    Minute  = DTime.minute
    prj     = 'GPM.KuPR'
    prdLv   = 'L2'
    prdVer  = '03'
    baseDir = "/tank/utsumi/PMM/WNP.261x265/%s/%s/%s"%(prj,prdLv,prdVer)
    dataDir = baseDir + "/%04d/%02d"%(Year,Mon)
    dataPath= dataDir + "/pr.%02d.%02d.%02d.%dx%d"%(Day,Hour,Minute,ny,nx)
    try:
      a2pr    = fromfile(dataPath, float32).reshape(ny,nx)
      a2pr    = ma.masked_equal(a2pr,miss)*60.*60.   # mm/s -> mm/h
      print "Read",dattype,DTime
    except IOError:
      raise MyIOException()

  elif dattype =="GMI":
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour
    Minute  = DTime.minute
    prj     = 'GPM.GMI'
    prdLv   = 'L2'
    prdVer  = '03'
    baseDir = "/tank/utsumi/PMM/WNP.261x265/%s/%s/%s"%(prj,prdLv,prdVer)
    dataDir = baseDir + "/%04d/%02d"%(Year,Mon)
    dataPath= dataDir + "/pr.%02d.%02d.%02d.%dx%d"%(Day,Hour,Minute,ny,nx)
    try:
      a2pr    = fromfile(dataPath, float32).reshape(ny,nx)
      a2pr    = ma.masked_equal(a2pr,miss)*60.*60.   # mm/s -> mm/h
      print "Read",dattype,DTime
    except IOError:
      raise MyIOException()
    
  else:
    print "check ret_mmh"; sys.exit()
  return a2pr
#*****************
for dattype in ldattype:
  #---------
  if   dattype.split(".")[0] =="GSMaP":
    import myfunc.IO.GSMaP as GSMaP
    gsmap = GSMaP.GSMaP(prj="standard", ver="v6", BBox=BBox)
    LatOrg= gsmap.Lat
    LonOrg= gsmap.Lon
    us    = Regrid.UpScale()
    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
  elif dattype.split(".")[0] == "IMERG":
    import myfunc.IO.IMERG as IMERG
    imerg = IMERG.IMERG(PRD="PROD",VER="V03",crd="sa", BBox=BBox)
    LatOrg= imerg.Lat
    LonOrg= imerg.Lon
    us    = Regrid.UpScale()
    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
 
  elif dattype in ["KuPR","GMI"]:
    pass
  else:
    print "check dattype", dattype; sys.exit()

  #---------

  for YM in lYM:
    Year   = YM[0]
    Mon    = YM[1]
    iDay   = 1
    eDay   = calendar.monthrange(Year,Mon)[1]
    eMinute= ret_eMinute(dattype)
  
    iDTime = datetime(Year,Mon,iDay,0,0)
    eDTime = datetime(Year,Mon,eDay,23,eMinute)
 
    dDTime = ret_dDTime(dattype)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    #lDTime = lDTime[:5]   # test
  
    # No Data list for Cloud
    f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
    lNoData = [s.strip() for s in f.readlines()]
    f.close()
  
    # Initialize

    dpr = {icl:deque([]) for icl in lcltype} 
    dra = {icl:deque([]) for icl in lcltype} 
    for DTime in lDTime:
      print DTime

      # load Cloud Type
      if DTime.minute == 30:
        DTimeCL = DTime + dDTime
      else:
        DTimeCL = DTime
      try:
        if   clVer  == "JMA1":
          a2cl = cl.loadData(DTimeCL, DType="clc")
        elif clVer[:5] == "MyWNP":
          a2cl = cl.loadData(DTimeCL)

      except IOError:
        Year = DTimeCL.year
        Mon  = DTimeCL.month
        Day  = DTimeCL.day
        Hour = DTimeCL.hour
        if "%04d-%02d-%02d-%02d"%(Year,Mon,Day,Hour) in lNoData:
          print "skip", DTimeCL
        else:
          print "check cloud data"
          print DTimeCL
          print "stop"
          sys.exit()

      #-- load precipitation 
      try:
        a2pr    = ret_mmh(DTime, dattype=dattype)    # mm/h, forward
        a2pr    = ma.masked_equal(a2pr, miss)
      except MyIOException as error:
        print "MyIOException: Skip",DTime
        continue

      #-- load RA -----
      a2ra  = ret_mmh(DTime, dattype="RA") 
      a2ra  = ma.masked_equal(a2ra, miss) 

      #-- mask by cloud
      a2mask = a2ramask
      a2mask = ma.masked_where((a2ra< vmin)&(a2pr<vmin), a2mask)

      for icl in lcltype: 
        clid = dclid[icl]
        a2mask1 = ma.masked_where(a2cl !=clid, a2mask)
        a2prtmp = ma.masked_where(a2mask1.mask, a2pr)
        a2ratmp = ma.masked_where(a2mask1.mask, a2ra)

        dpr[icl].extend(a2prtmp.compressed())
        dra[icl].extend(a2ratmp.compressed())

    #-- Save --------- 
    for icl in lcltype:
      #baseDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265"
      baseDir  = ibaseDir
      sDir     = baseDir + "/VsRA.CL.%s"%(dattype)
      prPath   = sDir + "/%s.%04d.%02d.%s.bn"%(dattype,Year,Mon,dclShortName[icl])
      raPath   = sDir + "/RA.%04d.%02d.%s.bn"%(Year,Mon,dclShortName[icl])

      util.mk_dir(sDir)
      array(dpr[icl], float32).tofile(prPath)
      array(dra[icl], float32).tofile(raPath)
      print prPath

