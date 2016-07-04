from numpy import *
from datetime import datetime, timedelta
import calendar
import os, sys
import ioPMM
import TagPMM
import myfunc.util as util

iYear = 2014
eYear = 2014
lYear = range(iYear,eYear+1)
#lMon  = range(1,12+1)
#lMon  = [4,5,6,7,8,9,10,11,12]
lMon  = [12]
lHour = range(0,23+1)
ltag  = ["tc","fbc","c"]
ltagAll= ltag+["ot"]
prtype1 = "RA"

#prtype2 = "GPM.KuPR" 
prtype2 = "GPM.KaPR" 
#prtype2 = "GSMaP" 
#prtype2 = "RA" 

Trc  = "ALL"
#Trc  = "GPM.KuPR"

IO   = ioPMM.ioRAdom()
TAG  = TagPMM.TagPMM()
if ((prtype1 == "GSMaP")or(prtype2=="GSMaP")):
  IO.init_GSMaP(prj="standard",ver="v6") 

ny,nx = IO.ny, IO.nx


#--------------------------------
def loadPrcp(prtype, DTime):
  if   prtype in ["GPM.KuPR","GPM.KaPR"]:
    try:
      return IO.loadGPMmmh(DTime, prj=prtype, maskflag=True, return_dummy=False)
    except TabError:
      return "NoObs"

  elif prtype == "GSMaP":
    return IO.loadGSMaPmmh(DTime, maskflag=True)

  elif prtype == "RA":
    return IO.loadRAmmh(DTime, maskflag=True)

def retTag(DTime, ltag):
  if   DTime.hour in [3,9,15,21]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=3))
  elif DTime.hour in [4,10,16,22]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=2))
  elif DTime.hour in [5,11,17,23]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=1))
  elif DTime.hour in [6,12,18,0]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=0))
  elif DTime.hour in [7,13,19,1]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=-1))
  elif DTime.hour in [8,14,20,2]:
    dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=-2))
  return dTag 

def ret_oPath(Year,Mon,tag):
  sDir     = os.path.join(baseDir, "%04d"%(Year))
  sumPath1 = os.path.join(sDir, "sum.%s.%s.%04d%02d.%dx%d"%(prtype1,tag,Year,Mon,ny,nx))
  sumPath2 = os.path.join(sDir, "sum.%s.%s.%04d%02d.%dx%d"%(prtype2,tag,Year,Mon,ny,nx))
  numPath  = os.path.join(sDir, "num.%s.%04d%02d.%dx%d"%(tag,Year,Mon,ny,nx))
  return sDir, sumPath1, sumPath2, numPath  
                              
#---------------------------  -----
#-- inititalze Tag --         
dTag = {}
#--------------------

for Year,Mon in [[Year,Mon] for Year in lYear for Mon in lMon]:
  iDay = 1
  eDay = calendar.monthrange(Year,Mon)[1]
#  eDay = iDay
  lDay = range(iDay,eDay+1)

  #-- Initialize ----
  da2sum1 = {}
  da2sum2 = {}
  da2num  = {}
  for tag in ltagAll+["plain"]:
    da2sum1[tag] = zeros([ny,nx],float32)
    da2sum2[tag] = zeros([ny,nx],float32)
    da2num [tag] = zeros([ny,nx],float32)
  #------------------
  for Day, Hour in [[Day,Hour] for Day in lDay for Hour in lHour]:
    DTime = datetime(Year,Mon,Day,Hour)
    #--- caution ----------------
    if DTime >= datetime(2014,12,31,21):
      continue
    #----------------------------
    print Year,Mon,Day,Hour
    #-- load track --
    if Trc != "ALL":
      a2trc = loadPrcp(Trc, DTime)  
      if type(a2trc)==str:
        print "No track in the domain",DTime
        continue

    #-- load Prcp ---
    a2pr2 = loadPrcp(prtype2, DTime)
    if type(a2pr2)==str:
      print "No GPM Obs in the domain",DTime
      continue
    a2pr1 = loadPrcp(prtype1, DTime)

    #-- matching mask --
    if Trc != "ALL":
      Mask = a2pr1.mask + a2pr2.mask + a2trc.mask
    else:
      Mask = a2pr1.mask + a2pr2.mask

    a2pr1 = ma.masked_array(a2pr1, Mask)
    a2pr2 = ma.masked_array(a2pr2, Mask)
    a2num = ma.masked_array(ones([ny,nx],float32), Mask)

    #-- load Tag --
    if dTag =={}:
      dTag = retTag(DTime, ltag)

    elif Hour in [3,9,15,21]:
      dTag = TAG.dictTagFrac(ltag, DTime+timedelta(hours=3))

    #-- calc plain -------
    da2sum1["plain"] = da2sum1["plain"] + a2pr1.filled(0.0)
    da2sum2["plain"] = da2sum2["plain"] + a2pr2.filled(0.0)
    da2num ["plain"] = da2num ["plain"] + a2num.filled(0.0)

    #-- calc tagged ------
    for tag in ltagAll:
      da2sum1[tag] = da2sum1[tag] + (a2pr1 *dTag[tag]).filled(0.0)
      da2sum2[tag] = da2sum2[tag] + (a2pr2 *dTag[tag]).filled(0.0)
      da2num [tag] = da2num [tag] + (a2num *dTag[tag]).filled(0.0)

  #-- Write to file ------ 
  baseDir = "/tank/utsumi/PMM/RAdom/%s.vs.%s.on.%s.%dclass"%(prtype1, prtype2, Trc, len(ltagAll))
    
  for tag in ltagAll+["plain"]:
    sDir     = ret_oPath(Year,Mon,tag)[0]
    sumPath1 = ret_oPath(Year,Mon,tag)[1]
    sumPath2 = ret_oPath(Year,Mon,tag)[2]
    numPath  = ret_oPath(Year,Mon,tag)[3]

    util.mk_dir(sDir)
   
    da2sum1[tag].tofile(sumPath1)
    da2sum2[tag].tofile(sumPath2)
    da2num [tag].tofile(numPath )
    print numPath     

