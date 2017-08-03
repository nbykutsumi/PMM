from numpy import *
from myfunc.IO import GPyM
from datetime import datetime, timedelta
from PMM.pmm_fsub import *
from collections import deque
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import os

prj     = 'GPM.KuPR'
prdLv   = 'L2'
prdVer  = '03'
#var     = 'NS/SLV/precipRateESurface'
var     = 'NS/SLV/precipRate'
BBox    = [[-0.1, 113.875],[52.1, 180.125]]
clVer = "MyWNP.M.3"
iYM    = [2014,4]
eYM    = [2014,4]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
miss_out= -9999.
nz      = 176 
nw      = 49   # Ku

rootDir = "/home/utsumi/mnt/well.share"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)

[[lllat,lllon],[urlat,urlon]] = cl.BBox
dlat, dlon   = cl.dLat, cl.dLon
gpm    = GPyM.GPM(prj, prdLv, prdVer)


for YM in lYM:
  Year   = YM[0]
  Mon    = YM[1]
#  eDay   = calendar.monthrange(Year,Mon)[1]
  eDay   = 1
  iDTime = datetime(Year,Mon,1,   0,0)
  #eDTime = datetime(Year,Mon,eDay,23,0) 
  eDTime = datetime(Year,Mon,eDay,5,0) 
  dDTime = timedelta(hours=1)
  lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

  # Initialize outputs
  aoLat  = deque([])
  aoLon  = deque([])
  aoTim  = deque([]) 
  aoPrc  = deque([])
  aoCld  = deque([])

  for DTime in lDTime:
    DTimeCL= DTime + timedelta(minutes=30)
    try:
      vout   = gpm('NS/SLV/precipRate',         DTime, DTime+timedelta(hours=1), BBox=BBox)
      sout   = gpm('NS/SLV/precipRateESurface', DTime, DTime+timedelta(hours=1), BBox=BBox)
    except IndexError:
      print "-"*50
      print DTime, "No obs at the BBox"
      print "-"*50
      continue
    a1Lat    = vout.lat.flatten()
    a1Lon    = vout.lon.flatten()
    a1DTime  = (resize(vout.dtime, (nw, len(vout.dtime))).T).flatten()
    a2Prc    = vout.data.reshape(-1, vout.data.shape[-1])

    a1msk    = ma.masked_less_equal(sout.data.flatten(), 0.0).mask
    a1lat    = ma.masked_where(a1msk, a1Lat).compressed()
    a1lon    = ma.masked_where(a1msk, a1Lon).compressed()
    a1DTime  = ma.masked_where(a1msk, a1DTime).compressed()

    a2prc    = ma.masked_where(resize( a1msk.reshape(1,len(a1msk)), a2Prc.shape) , a2Prc).compressed().reshape(-1,nz)

    # pickup cloud type
    a2CL    = cl.loadData(DTimeCL)
    a1cl    = (pmm_fsub.pickup_data(a2CL.T, lllat, lllon, urlat, urlon, dlat, dlon, a1lon, a1lat, miss_out).T)[0]

    # Stack
    aoLat.extend(a1lat)
    aoLon.extend(a1lon)
    aoTim.extend(a1time)
    aoPrc.extend(a1prc)
    aoCld.extend(a1cl)
  
  
