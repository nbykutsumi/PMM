from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
from collections import deque
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import myfunc.util as util
import calendar
import sys, os
from mpl_toolkits.basemap import Basemap


#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP3"
clVer = "MyWNP.M.3"

#iYM    = [2014,4]
#eYM    = [2015,6]
iYM    = [2014,5]
eYM    = [2015,6]


lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
ldattype = ["IMERG.IR"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]


BBox    = [[-0.1, 113.875],[52.1, 180.125]]
ny,nx   = 261, 265
miss    = -9999.

rootDir = "/home/utsumi/mnt/wellshare"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)

if clVer[5:8]==".M.":
  MidFlag = True
else:
  MidFlag = False
print "*"*50+"\n"
print "MidFlag=",MidFlag
print "*"*50


lcltype = cl.licl
#lcltype = [3]
ncltype = len(lcltype)
dclName = cl.dclName
dclShortName = cl.dclShortName
dclid   = cl.dclid
LatUp   = cl.Lat
LonUp   = cl.Lon



landfracDir = "/home/utsumi/mnt/wellshare/PMM/WNP.261x265/MASK"
landfracPath= landfracDir + "/landfrac.%dx%d"%(cl.ny, cl.nx)
a2lndfrc    = fromfile(landfracPath, float32).reshape(261,265)

a2lndmask   = ma.masked_less(a2lndfrc,1.0).mask  # mask sea&coast
a2seamask   = ma.masked_greater(a2lndfrc,0.0).mask # mask land&coast
a2cstmask   = ma.masked_equal(a2lndfrc, 1.0).mask \
             +ma.masked_equal(a2lndfrc, 0.0).mask  # mask land&sea

dLSmask     = {"lnd":a2lndmask, "sea":a2seamask, "cst":a2cstmask}
llndsea     = ["lnd","sea","cst"]

# Lat and Lon array
a2lon, a2lat = meshgrid(LonUp, LatUp)

# Y and X
a2x, a2y     = meshgrid(range(nx),range(ny))


#*****************
class MyIOException(Exception): pass

#*****************
def check_nodata(lNoData, DTime):
    DTime = DTime
    Year = DTime.year
    Mon  = DTime.month
    Day  = DTime.day
    Hour = DTime.hour
    Minute = DTime.minute
    if "%04d-%02d-%02d-%02d-%02d"%(Year,Mon,Day,Hour,Minute) in lNoData:
      print "skip", DTime
      return "NoData"
    else:
      print DTime,"Not listed in"
      print nodataPath
      return "error"

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
    var    = "precipitationCal"
    DTime0 = DTime
    DTime1 = DTime + timedelta(minutes=30)
    a3prOrg    = empty([2,imerg.ny, imerg.nx])
    a3prOrg[0] = imerg.load_mmh(DTime0, var=var)  # mm/h, forward
    a3prOrg[1] = imerg.load_mmh(DTime1, var=var)  # mm/h, forward
    a2prOrg  = (ma.masked_less(a3prOrg, 0.0).mean(axis=0)).filled(miss)
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG.IR":
    var    = "IRprecipitation"
    DTime0 = DTime
    DTime1 = DTime + timedelta(minutes=30)
    a3prOrg    = empty([2,imerg.ny, imerg.nx])
    a3prOrg[0] = imerg.load_mmh(DTime0, var=var)  # mm/h, forward
    a3prOrg[1] = imerg.load_mmh(DTime1, var=var)  # mm/h, forward
    a2prOrg  = (ma.masked_less(a3prOrg, 0.0).mean(axis=0)).filled(miss)
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)


  elif dattype =="IMERG.MW":
    var    = "HQprecipitation"
    DTime0 = DTime
    DTime1 = DTime + timedelta(minutes=30)
    a3prOrg    = empty([2,imerg.ny, imerg.nx])
    a3prOrg[0] = imerg.load_mmh(DTime0, var=var)  # mm/h, forward
    a3prOrg[1] = imerg.load_mmh(DTime1, var=var)  # mm/h, forward
    a2prOrg  = (ma.masked_less(a3prOrg, 0.0).mean(axis=0)).filled(miss)
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="RA":
    a2prOrg = ra.loadForward_mmh(DTime,mask=True).filled(miss)    # mm/h, forward
    a2pr    = usRa.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="KuPR":
    DTime0   = DTime
    DTime1   = DTime + timedelta(minutes=30)

    Year0    = DTime0.year
    Mon0     = DTime0.month
    Day0     = DTime0.day
    Hour0    = DTime0.hour
    Minute0  = DTime0.minute

    Year1    = DTime1.year
    Mon1     = DTime1.month
    Day1     = DTime1.day
    Hour1    = DTime1.hour
    Minute1  = DTime1.minute

    prj     = 'GPM.KuPR'
    prdLv   = 'L2'
    prdVer  = '03'
    baseDir = "/tank/utsumi/PMM/WNP.261x265/%s/%s/%s"%(prj,prdLv,prdVer)

    dataDir0 = baseDir  + "/%04d/%02d"%(Year,Mon)
    dataPath0= dataDir0 + "/pr.%02d.%02d.%02d.%dx%d"%(Day0,Hour0,Minute0,ny,nx)

    dataDir1 = baseDir  + "/%04d/%02d"%(Year,Mon)
    dataPath1= dataDir1 + "/pr.%02d.%02d.%02d.%dx%d"%(Day1,Hour1,Minute1,ny,nx)

    exist0   = os.path.exists(dataPath0)
    exist1   = os.path.exists(dataPath1)

    if (exist0==True)&(exist1==True):
      a2pr0    = fromfile(dataPath0, float32).reshape(ny,nx)
      a2pr1    = fromfile(dataPath1, float32).reshape(ny,nx)
      a2pr    = ma.masked_less(
                array([a2pr0,a2pr1]),0.0
                ).mean(axis=0)*60.*60.  # mm/s --> mm/h
    elif (exist0==True)&(exist1==False):
      a2pr    = ma.masked_less(
                fromfile(dataPath0, float32).reshape(ny,nx), 0.0
                )*60.*60.  # mm/s --> mm/h
    elif (exist0==False)&(exist1==True):
      a2pr    = ma.masked_less(
                fromfile(dataPath1, float32).reshape(ny,nx), 0.0
                )*60.*60.  # mm/s --> mm/h
    else:
      raise MyIOException()

  elif dattype =="GMI":
    DTime0   = DTime
    DTime1   = DTime + timedelta(minutes=30)

    Year0    = DTime0.year
    Mon0     = DTime0.month
    Day0     = DTime0.day
    Hour0    = DTime0.hour
    Minute0  = DTime0.minute

    Year1    = DTime1.year
    Mon1     = DTime1.month
    Day1     = DTime1.day
    Hour1    = DTime1.hour
    Minute1  = DTime1.minute

    prj     = 'GPM.GMI'
    prdLv   = 'L2'
    prdVer  = '03'
    baseDir = "/tank/utsumi/PMM/WNP.261x265/%s/%s/%s"%(prj,prdLv,prdVer)

    dataDir0 = baseDir  + "/%04d/%02d"%(Year,Mon)
    dataPath0= dataDir0 + "/pr.%02d.%02d.%02d.%dx%d"%(Day0,Hour0,Minute0,ny,nx)

    dataDir1 = baseDir  + "/%04d/%02d"%(Year,Mon)
    dataPath1= dataDir1 + "/pr.%02d.%02d.%02d.%dx%d"%(Day1,Hour1,Minute1,ny,nx)

    exist0   = os.path.exists(dataPath0)
    exist1   = os.path.exists(dataPath1)

    if (exist0==True)&(exist1==True):
      a2pr0    = fromfile(dataPath0, float32).reshape(ny,nx)
      a2pr1    = fromfile(dataPath1, float32).reshape(ny,nx)
      a2pr    = ma.masked_less(
                array([a2pr0,a2pr1]), 0.0
                ).mean(axis=0)*60.*60.  # mm/s --> mm/h

    elif (exist0==True)&(exist1==False):
      a2pr    = ma.masked_less(
                fromfile(dataPath0, float32).reshape(ny,nx), 0.0
                )*60.*60.  # mm/s --> mm/h

    elif (exist0==False)&(exist1==True):
      a2pr    = ma.masked_less(
                fromfile(dataPath1, float32).reshape(ny,nx), 0.0
                )*60.*60.  # mm/s --> mm/h
    else:
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
 
    iDTime = datetime(Year,Mon,iDay,0,0)
    eDTime = datetime(Year,Mon,eDay,23,0)
 
    dDTime = timedelta(hours=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    #lDTime = lDTime[:24*2]   # test
  
    # No Data list for Cloud
    nodataPath = cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon)
    f = open(nodataPath, "r")
    lNoData = [s.strip() for s in f.readlines()]
    f.close()
  
    # Initialize
    apr  = deque([])
    aku  = deque([])
    acl  = deque([])
    als  = deque([])
    ay   = deque([])
    ax   = deque([])
    alat = deque([])
    alon = deque([])
    #atime= deque([])
    aday = deque([])
    ahour= deque([])
    amin = deque([])


    for DTime in lDTime:
    #for DTime in lDTime[-10:-1]:

      # load Cloud Type
      if MidFlag == True:
        DTimeCL = DTime + timedelta(minutes=30)
      else:
        DTimeCL = DTime

      try:
        if   clVer  == "JMA1":
          a2cl = cl.loadData(DTimeCL, DType="clc")
        elif clVer[:5] == "MyWNP":
          a2cl = cl.loadData(DTimeCL)

      except IOError:
        print "NoCloud",DTimeCL
        check = check_nodata(lNoData, DTimeCL)
        if check=="NoData": continue
        else:sys.exit()

      #-- load KuPR precipitation --
      try:
        a2ku   = ret_mmh(DTime, dattype="KuPR")    # mm/h, forward
        a2ku   = ma.masked_equal(a2ku, miss)
      except MyIOException as error:
        print "MyIOException: Skip","KuPR",DTime
        continue

      #-- load precipitation 
      try:
        a2pr    = ret_mmh(DTime, dattype=dattype)    # mm/h, forward
        #a2pr    = ma.masked_equal(a2pr, miss)
      except MyIOException as error:
        print "MyIOException: Skip",dattype,DTime
        continue

      print DTime
      #-- general mask
      a2mask = a2ku.mask

      #-- mask 
      a2prtmp = ma.masked_where(a2mask, a2pr)
      a2kutmp = ma.masked_where(a2mask, a2ku)
      a2cltmp = ma.masked_where(a2mask, a2cl)

      aku.extend(a2kutmp.compressed())
      apr.extend(a2prtmp.compressed())
      acl.extend(a2cltmp.compressed())

      # Land sea info
      als.extend(ma.masked_array(a2lndfrc, a2mask).compressed()) 
      
      ## X and Y
      #ax.extend(ma.masked_array(a2x, a2mask).compressed())
      #ay.extend(ma.masked_array(a2y, a2mask).compressed())

      # Lon and Lat
      alon.extend(ma.masked_array(a2lon, a2mask).compressed())
      alat.extend(ma.masked_array(a2lat, a2mask).compressed())

      # Time
      #stime= "%04d,%02d,%02d,%02d,%02d"\
      #  %(DTimeCL.year, DTimeCL.month, DTimeCL.day, DTimeCL.hour, DTimeCL.minute)
      #atime.extend(deque([stime]*(ma.count(a2mask)-sum(a2mask))))

      counttemp = ma.count(a2mask)-sum(a2mask)
      aday. extend(deque([DTimeCL.day]*counttemp))
      ahour.extend(deque([DTimeCL.hour]*counttemp))
      amin.extend(deque([DTimeCL.minute]*counttemp))

    #-- Save --------- 
    baseDir = ibaseDir
    sDir  = baseDir + "/VsKuPR.LocTime/%04d"%(Year)

    PrPath   = sDir + "/%s.%04d.%02d.bn"%(dattype,Year,Mon)
    KuPath   = sDir + "/KuPR.%04d.%02d.bn"%(Year,Mon)
    CLPath   = sDir + "/CL.%04d.%02d.bn"%(Year,Mon)
    LSPath   = sDir + "/LS.%04d.%02d.bn"%(Year,Mon)
    LonPath  = sDir + "/Lon.%04d.%02d.bn"%(Year,Mon)
    LatPath  = sDir + "/Lat.%04d.%02d.bn"%(Year,Mon)
    DayPath  = sDir + "/Day.%04d.%02d.bn"%(Year,Mon)
    HourPath = sDir + "/Hour.%04d.%02d.bn"%(Year,Mon)
    MinPath  = sDir + "/Min.%04d.%02d.bn"%(Year,Mon)

    util.mk_dir(sDir)

    asarray(apr, dtype="float32").tofile(PrPath)
    asarray(aku, dtype="float32").tofile(KuPath)
    asarray(acl, dtype="int32").tofile(CLPath)
    asarray(als, dtype="float32").tofile(LSPath)
    asarray(alon, dtype="float32").tofile(LonPath)
    asarray(alat, dtype="float32").tofile(LatPath)
    asarray(aday,  dtype="int32").tofile(DayPath)
    asarray(ahour, dtype="int32").tofile(HourPath)
    asarray(amin,  dtype="int32").tofile(MinPath)


    print "*"*50
    print "Save!!"
    print PrPath
    print "*"*50

