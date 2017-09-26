from numpy import *
from datetime import datetime, timedelta
from myfunc.regrid import Regrid
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import calendar
import sys, os

from mpl_toolkits.basemap import Basemap

#clVer  = "MyWNP1"
#clVer  = "MyWNP2"
#clVer  = "MyWNP3"
clVer  = "MyWNP.M.3"
#clVer  = "JMA1"

iYM    = [2014,4]
eYM    = [2015,6]
#iYM    = [2015,6]
#eYM    = [2015,6]


lYM    = util.ret_lYM(iYM, eYM)
lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
print lYM 

#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["KuPR","GMI"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
ldattype = ["IMERG","IMERG.IR","IMERG.MW"]
#ldattype = ["KuPR","GMI"]

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

#lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1,1) + range(10,50+1,2) +[999]  # Maxs of bin range
lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1,1) + range(10,18+1,2) + range(20, 40+1,4) +[999]  # Maxs of bin range
nbinPr  = len(lbinPr)

rootDir = "/home/utsumi/mnt/well.share"
if   clVer == "JMA1":
  cl      = CLOUDTYPE.CloudWNP()
  lcltype = range(0,7+1)
  dclid   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}

elif clVer[:5] == "MyWNP":
  ver     = clVer[5:]
  cl      = CLOUDTYPE.MyCloudWNP(ver=ver)
  lcltype = cl.licl
  dclid   = cl.dclid
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

if clVer[5:8]==".M.":
  MidFlag = True
else:
  MidFlag = False
print "*"*50+"\n"
print "MidFlag=",MidFlag
print "*"*50

ncltype = len(lcltype)

# Cloud Type
LatUp = cl.Lat
LonUp = cl.Lon


#*****************
class MyIOException(Exception): pass


#*****************
def ret_mmh(DTime, dattype):
  if   dattype =="GSMaP":
    a2prOrg = ma.masked_less(gsmap.load_mmh(DTime),0.0).filled(miss)    # mm/h, forward
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)
  elif dattype =="GSMaP.MW":
    a2sateinfo = gsmap.load_sateinfo(DTime)
    a2prOrg = ma.masked_where(a2sateinfo <=0,
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
    a2pr    = us.upscale(a2prOrg, pergrid=False, miss_in=miss, miss_out=miss)

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
    #prdVer  = '03'
    prdVer  = '05'
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
                #).mean(axis=0)*60.*60.  # mm/s --> mm/h
                ).mean(axis=0)  # mm/h
    elif (exist0==True)&(exist1==False):
      a2pr    = ma.masked_less(
                fromfile(dataPath0, float32).reshape(ny,nx), 0.0
                #)*60.*60.  # mm/s --> mm/h
                )           #  mm/h
    elif (exist0==False)&(exist1==True):
      a2pr    = ma.masked_less(
                fromfile(dataPath1, float32).reshape(ny,nx), 0.0
                #)*60.*60.  # mm/s --> mm/h
                )          # mm/h
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
    #prdVer  = '03'
    prdVer  = '05'
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
      #a2pr    = ma.masked_less(
      #          array([a2pr0,a2pr1]), 0.0
      #          ).mean(axis=0)*60.*60.  # mm/s --> mm/h
      a2pr    = ma.masked_less(
                array([a2pr0,a2pr1]), 0.0
                ).mean(axis=0)  # mm/h


    elif (exist0==True)&(exist1==False):
      a2pr    = ma.masked_less(
                fromfile(dataPath0, float32).reshape(ny,nx), 0.0
                #)*60.*60.  # mm/s --> mm/h
                )           # mm/h
    elif (exist0==False)&(exist1==True):
      a2pr    = ma.masked_less(
                fromfile(dataPath1, float32).reshape(ny,nx), 0.0
                #)*60.*60.  # mm/s --> mm/h
                )           # mm/h
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
    #gsmap = GSMaP.GSMaP(prj="standard", ver="v6", BBox=BBox)
    gsmap = GSMaP.GSMaP(prj="standard", ver="v7", BBox=BBox, compressed=True)
    LatOrg= gsmap.Lat
    LonOrg= gsmap.Lon
    us    = Regrid.UpScale()
    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
  elif dattype.split(".")[0] == "IMERG":
    import myfunc.IO.IMERG as IMERG
    #imerg = IMERG.IMERG(PRD="PROD",VER="V03",crd="sa", BBox=BBox)
    imerg = IMERG.IMERG(PRD="PROD",VER="V04",crd="sa", BBox=BBox)
    LatOrg= imerg.Lat
    LonOrg= imerg.Lon
    us    = Regrid.UpScale()
    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)
  
  elif dattype == "RA":
    import myfunc.IO.RadarAMeDAS as RadarAMeDAS
    ra    = RadarAMeDAS.RadarAMeDAS(prj="ra_0.01")
    LatOrg= ra.Lat
    LonOrg= ra.Lon
  
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
    #eMinute= ret_eMinute(dattype)
  
    iDTime = datetime(Year,Mon,iDay,0,0)
    eDTime = datetime(Year,Mon,eDay,23,0)
  
    dDTime = timedelta(hours=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    #lDTime = lDTime[:24]   # test
  
    # No Data list for Cloud
    f = open(cl.baseDir + "/%04d%02d/list.nodata.txt"%(Year,Mon), "r")
    lNoData = [s.strip() for s in f.readlines()]
    f.close()
  
    # Initialize
    da3sum = {binPr: zeros([ncltype,ny,nx],float32) for binPr in lbinPr}
    da3num = {binPr: zeros([ncltype,ny,nx],int32  ) for binPr in lbinPr} # Int32 !!
    
    a2oneint= ones([ny,nx],int32)
  
    for DTime in lDTime:
      print dattype,DTime
      try:
        a2pr    = ret_mmh(DTime, dattype=dattype)    # mm/h, forward
        a2pr    = ma.masked_equal(a2pr, miss)
      except MyIOException as error:
        print "MyIOException: Skip",DTime
        continue
  
  
      # load Cloud Type 
      if MidFlag == True:
        DTimeCL = DTime + timedelta(minutes=30)
      elif MidFlag == False:
        DTimeCL = DTime
      else:
        print "check MidFlag", MidFlag, sys.exit()
      try:
        if clVer   == "JMA1":
          a2cl = cl.loadData(DTimeCL, DType="clc")
        elif clVer[:5] == "MyWNP":
          a2cl = cl.loadData(DTimeCL)

      except IOError:
        Year = DTimeCL.year
        Mon  = DTimeCL.month
        Day  = DTimeCL.day
        Hour = DTimeCL.hour
        Minute = DTimeCL.minute
        if "%04d-%02d-%02d-%02d-%02d"%(Year,Mon,Day,Hour,Minute) in lNoData:
          print "skip", DTimeCL
          continue
        else:
          print "check cloud data"
          print DTimeCL
          print "stop"
          sys.exit() 
      #-----------------
      for binPr in lbinPr:
        if binPr==0.0:
          a2tmp1  = ma.masked_greater(a2pr, binPr)
        else:
          a2tmp1  = ma.masked_greater_equal(a2pr, binPr)
    
        for cltype in lcltype:
          clid   = dclid[cltype]
          a2tmp2 = ma.masked_where(a2cl !=clid, a2tmp1)
      
          da3sum[binPr][cltype] = da3sum[binPr][cltype] + a2tmp2.filled(0.0)
          da3num[binPr][cltype] = da3num[binPr][cltype] + ma.masked_where(a2tmp2.mask, a2oneint).filled(0)
  
    # Save Monthly output
    #rootDir = "/tank/utsumi"
    rootDir = "/home/utsumi/mnt/wellshare"
    if   clVer == "JMA1":
      baseDir = rootDir + "/PMM/WNP.261x265/CL.JMA"
    elif clVer[:5] == "MyWNP":
      baseDir = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)

    oDir    = baseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
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
  
