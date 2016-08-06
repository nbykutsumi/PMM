from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
import sys

iYM    = [2014,4]
#eYM    = [2014,4]
eYM    = [2014,12]
lYM    = util.ret_lYM(iYM, eYM)
#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
ldattype = ["KuPR","GMI","GSMaP","GSMaP.IR"]
#ldattype = ["GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]

BBox    = [[-0.1, 113.875],[52.1, 180.125]]
miss    = -9999.
#*****************
def ret_dDTime(dattype):
  if dattype.split(".")[0]=="GSMaP":
    return timedelta(hours=1)
  elif dattype == "RA":
    return timedelta(hours=1)
  elif dattype in ["IMERG","IMERG.IR","IMERG.MW","KuPR","GMI"]:
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
    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)
  elif dattype =="GSMaP.MW":
    a2sateinfo = gsmap.load_sateinfo(DTime)
    a2prOrg = ma.masked_where(a2sateinfo <=0,
               ma.masked_less(gsmap.load_mmh(DTime),0.0)
              ).filled(miss)    # mm/h, forward

    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="GSMaP.IR":
    a2sateinfo = gsmap.load_sateinfo(DTime)
    a2prOrg = ma.masked_where(a2sateinfo >=0,
                ma.masked_less(gsmap.load_mmh(DTime),0.0)
              ).filled(miss)    # mm/h, forward

    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime),0.0).filled(miss)    # mm/h, forward
    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG.IR":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime, var="IRprecipitation"),0.0).filled(miss)    # mm/h, forward
    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="IMERG.MW":
    a2prOrg = ma.masked_less(imerg.load_mmh(DTime, var="HQprecipitation"),0.0).filled(miss)    # mm/h, forward
    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

  elif dattype =="RA":
    a2prOrg = ra.loadForward_mmh(DTime,mask=True).filled(miss)    # mm/h, forward
    #a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

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
      a2prOrg = a2pr
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
      a2prOrg = a2pr
      print "Read",dattype,DTime
    except IOError:
      raise MyIOException()

  else:
    print "check ret_mmh"; sys.exit()

  return a2prOrg


#*****************
for dattype in ldattype:
  #---------
  if   dattype.split(".")[0] =="GSMaP":
    import myfunc.IO.GSMaP as GSMaP
    gsmap = GSMaP.GSMaP(prj="standard", ver="v6", BBox=BBox)
    LatOrg= gsmap.Lat
    LonOrg= gsmap.Lon
    #us    = Regrid.UpScale()
    #us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)

  elif dattype.split(".")[0] == "IMERG":
    import myfunc.IO.IMERG as IMERG
    imerg = IMERG.IMERG(PRD="PROD",VER="V03",crd="sa", BBox=BBox)
    LatOrg= imerg.Lat
    LonOrg= imerg.Lon
    #us    = Regrid.UpScale()
    #us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)

  elif dattype == "RA":
    import myfunc.IO.RadarAMeDAS as RadarAMeDAS
    ra    = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
    LatOrg= ra.Lat
    LonOrg= ra.Lon

    #us    = Regrid.UpScale()
    #us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)

  elif dattype in ["KuPR","GMI"]:
    pass
  else:
    print "check dattype", dattype; sys.exit()
  #---------
  if dattype in ["KuPR","GMI"]:
    ny = 280
    nx = 320
  else:
    ny     = len(LatOrg)
    nx     = len(LonOrg)

  lYM   = util.ret_lYM(iYM,eYM)    
  for YM in lYM:
   
    Year  = YM[0]
    Mon   = YM[1]
    iDay   = 1
    eDay   = calendar.monthrange(Year,Mon)[1]
    eMinute= ret_eMinute(dattype)
  
    iDTime = datetime(Year,Mon,iDay,0,0)
    eDTime = datetime(Year,Mon,eDay,23,eMinute)
   
    dDTime = ret_dDTime(dattype) 
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
  
    a2sum  = zeros([ny,nx],float32)
    a2num  = zeros([ny,nx],int32)
    a2ones = ones([ny,nx],int32)
    
    for DTime in lDTime:
      a2tmp = ret_mmh(DTime, dattype)
      a2sum = a2sum + ma.masked_less(a2tmp, 0.0).filled(0.0)
      a2num = a2num + ma.masked_where(a2tmp<0.0, a2ones).filled(0)
      print DTime
  
  
    a2acc = (ma.masked_invalid(a2sum/a2num)*24.0*eDay).filled(miss)
 
    sDir  = "/home/utsumi/test/AccEst"
    util.mk_dir(sDir)
    accPath = sDir + "/accEst.%s.%04d.%02d.%dx%d"%(dattype,Year,Mon,ny,nx)
    sumPath = sDir + "/sum.%s.%04d.%02d.%dx%d"%(dattype,Year,Mon,ny,nx)
    numPath = sDir + "/num.%s.%04d.%02d.%dx%d"%(dattype,Year,Mon,ny,nx)

    a2acc.astype(float32).tofile(accPath) 
    a2sum.astype(float32).tofile(sumPath) 
    a2num.astype(int32).tofile(numPath) 
    print accPath
    print sumPath
    print numPath

