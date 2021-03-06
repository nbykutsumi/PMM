from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import myfunc.util as util
import calendar
import sys

from mpl_toolkits.basemap import Basemap


iYM    = [2014,4]
eYM    = [2014,11]
#eYM    = [2014,12]
lYM    = util.ret_lYM(iYM, eYM)
#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
ldattype = ["GSMaP.IR","GSMaP.MW","IMERG"]

#dattype= "RA"
#dattype= "GSMaP"
#dattype= "GSMaP.MW"
#dattype= "GSMaP.IR"
#dattype= "IMERG"
#dattype= "KuPR"
#dattype= "GMI"

BBox    = [[-0.1, 113.875],[52.1, 180.125]]
ny,nx   = 261, 265
miss    = -9999.

#lbinPr  = [0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,999]  # Maxs of bin range
lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1,1) + range(10,50+1,2) +[999]  # Maxs of bin range
#lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7]
#lbinPr  = [3,5,7,9]  # Maxs of bin range
#lbinPr  = [22,24,26,28,32,34,36,38,42,44,46,48,50]  # Maxs of bin range
nbinPr  = len(lbinPr)
lcltype = range(0,7+1)
ncltype = len(lcltype)
dclid   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}

# Cloud Type
cl    = CLOUDTYPE.CloudWNP()
LatUp = cl.Lat
LonUp = cl.Lon

# RadarAMeDAS
ra    = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
LatOrgRA = ra.Lat
LonOrgRA = ra.Lon
usRa  = Regrid.UpScale()
usRa(LatOrgRA, LonOrgRA, LatUp, LonUp, globflag=False)
 

# Init Precip
#*****************
class MyIOException(Exception): pass


#*****************
def ret_dDTime(dattype):
  if dattype.split(".")[0]=="GSMaP":
    return timedelta(hours=1)
  elif dattype == "RA":
    return timedelta(hours=1)
  elif dattype in ["IMERG","KuPR","GMI"]:
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
  
  elif dattype == "IMERG":
    import myfunc.IO.IMERG as IMERG
    imerg = IMERG.IMERG(PRD="PROD",VER="V03",crd="sa", BBox=BBox)
    LatOrg= imerg.Lat
    LonOrg= imerg.Lon
    us    = Regrid.UpScale()
    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
  #elif dattype == "RA":
  #  import myfunc.IO.RadarAMeDAS as RadarAMeDAS
  #  ra    = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
  #  LatOrg= ra.Lat
  #  LonOrg= ra.Lon
  #
  #  us    = Regrid.UpScale()
  #  us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
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
  
    # No Data list for Cloud
    f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
    lNoData = [s.strip() for s in f.readlines()]
    f.close()
  
    # Initialize
    da3sum = {binPr: zeros([ncltype,ny,nx],float32) for binPr in lbinPr}
    da3num = {binPr: zeros([ncltype,ny,nx],int32  ) for binPr in lbinPr} # Int32 !!
    
    a2oneint= ones([ny,nx],int32)
  
    for DTime in lDTime:
      print DTime
      try:
        a2pr    = ret_mmh(DTime, dattype=dattype)    # mm/h, forward
        a2pr    = ma.masked_equal(a2pr, miss)
      except MyIOException as error:
        print "MyIOException: Skip",DTime
        continue
      #-- load RA -----
      a2ra  = ret_mmh(DTime, dattype="RA") 
      a2ra  = ma.masked_equal(a2ra, miss) 
 
      # load Cloud Type 
      if DTime.minute == 30:
        DTimeCL = DTime + dDTime
      else:
        DTimeCL = DTime
      try:
        a2cl = cl.loadData(DTimeCL, DType="clc")
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

      #-----------------
   
      for binPr in lbinPr:
        #if binPr==0.0:
        #  a2tmp1  = ma.masked_greater(a2pr, binPr)
        #else:
        #  a2tmp1  = ma.masked_greater_equal(a2pr, binPr)

        if binPr==0.0:
          a2cltmp  = ma.masked_where(a2ra>binPr,  a2cl)
        else:
          a2cltmp  = ma.masked_where(a2ra>=binPr, a2cl)

        for cltype in lcltype:
          clid    = dclid[cltype]
          a2tmp = ma.masked_where(a2cltmp !=clid, a2pr)
    
          da3sum[binPr][cltype] = da3sum[binPr][cltype] + a2tmp.filled(0.0)
          da3num[binPr][cltype] = da3num[binPr][cltype] + ma.masked_where(a2tmp.mask, a2oneint).filled(0)
  
    # Save Monthly output
    #baseDir = "/tank/utsumi/PMM/WNP.261x265"
    baseDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265"
    oDir    = baseDir + "/ByRA.CL.Pr.%s"%(dattype)
    util.mk_dir(oDir)
  
    # Bin file
    binPath = oDir + "/CloudType.txt"
    sout    = "\n".join(["%d:%d"%(cltype, dclid[cltype]) for cltype in lcltype]).strip()
    f=open(binPath,"w"); f.write(sout); f.close()
  
    for binPr in lbinPr:
      sumPath   = oDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
      numPath   = oDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  
      da3sum[binPr].tofile(sumPath)
      da3num[binPr].tofile(numPath)  
      print sumPath
