from ftplib import FTP
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket

hostname  = "arthurhou.pps.eosdis.nasa.gov"
irootDir = "/gpmallversions"
#irootDir = "/sm/730/gpmdata"

myhost = socket.gethostname()
if  myhost =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif myhost =="well":
    orootDir   = "/media/disk2/share/data/GPM"

#GPM/TRMM.TMI/L2A12/07/2014/

iYM       = [2017, 12]
eYM       = [2017, 12]

#iYM       = [2009, 6]
#eYM       = [2009, 6]
lYM       = util.ret_lYM(iYM, eYM)
#lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM




iDay      = 1
spec      = ["TRMM","TMI","2A-CLIM","gprof","V05","A"]
sate      = spec[0]
sensor    = spec[1]
prdName   = spec[2]
prj       = spec[3]
ver       = spec[4]
minorver  = spec[5]

#iDay      = 1
#spec      = ["GPM","GMI","1C","1C","V05","A"]
#sate      = spec[0]
#sensor    = spec[1]
#prdName   = spec[2]
#prj       = spec[3]
#ver       = spec[4]
#minorver  = spec[5]




myid      = "nbyk.utsumi@gmail.com"
mypass    = "nbyk.utsumi@gmail.com"

#----------------------------------
def mk_dir(sodir):
  try:
    os.makedirs(sodir)
  except OSError:
    pass


ftp = FTP(hostname)
ftp.login(myid, mypass)

#----------------------------------
for YM in lYM:
  Year, Mon = YM

  eDay = calendar.monthrange(Year,Mon)[1]
  lDay = range(iDay,eDay+1)
  for Day in lDay:
    #if (datetime(2014,7,23)<datetime(Year,Mon,Day))and(datetime(Year,Mon,Day) < datetime(2015,4,1)):continue

    #--- path and directory: Remote -----------
    iDir = irootDir + "/%s/%04d/%02d/%02d/%s"%(ver,Year,Mon,Day,prj)
    oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d"%(sate,sensor,prdName,ver,Year,Mon) 
    #GPM/TRMM.TMI/L2A12/07/2014/
    mk_dir(oDir)
    #--- list --------------
    lPath = ftp.nlst(iDir)
    print iDir
    print lPath
    for sPath in lPath:
      fName = os.path.basename(sPath)
      prdNameTmp, sateTmp, sensorTmp, algFullTmp, dtime, gNum, verFullTmp, sfx = fName.split('.')

      if prdNameTmp !=prdName: continue
      if sateTmp    !=sate   : continue
      if sensorTmp  !=sensor  : continue
      if verFullTmp !='%s%s'%(ver,minorver): continue
      oPath = oDir + "/" + sPath.split("/")[-1]
      ftp.retrbinary("RETR %s"%(sPath), open(oPath,"wb").write)
      print oPath
     
ftp.close()

