from ftplib import FTP
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket
from datetime import datetime, timedelta

hostname  = "arthurhou.pps.eosdis.nasa.gov"
irootDir = "/gpmallversions"
#irootDir = "/sm/730/gpmdata"

myhost = socket.gethostname()
if  myhost =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif myhost =="shui":
    orootDir   = "/work/hk01/PMM/NASA"
elif myhost =="well":
    orootDir   = "/media/disk2/share/data/GPM"

#GPM/TRMM.TMI/L2A12/07/2014/


iDTime  = datetime(2016,12,31)
#iDTime  = datetime(2017,2,25)
#eDTime  = datetime(2017,12,20)
eDTime  = datetime(2018,1,1)
dDTime  = timedelta(days=1)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)

#spec      = ["TRMM","TMI","2A-CLIM","gprof","V05","A"]
#spec      = ["GPM","GMI","1C","1C","V05","A"]
#spec      = ["GPM","GMI","2A-CLIM","gprof","V05","A"] # input=ECMWF
spec      = ["GPM","GMI","2A","gprof","V05","A"] # input=GANAL
#spec      = ["GPM","Ku","2A","radar","V06","A"]

sate      = spec[0]
sensor    = spec[1]
prdName   = spec[2]
prj       = spec[3]
ver       = spec[4]
minorver  = spec[5]


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
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    #--- path and directory: Remote -----------
    iDir = irootDir + "/%s/%04d/%02d/%02d/%s"%(ver,Year,Mon,Day,prj)
    #oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d"%(sate,sensor,prdName,ver,Year,Mon) 
    oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d/%02d"%(sate,sensor,prdName,ver,Year,Mon,Day) 
    #GPM/TRMM.TMI/L2A12/07/2014/
    mk_dir(oDir)
    #--- list --------------
    lPath = ftp.nlst(iDir)
    print iDir
    print lPath
    for sPath in lPath:
        fName = os.path.basename(sPath)
        prdNameTmp, sateTmp, sensorTmp, algFullTmp, dtime, gNum, verFullTmp, sfx = fName.split('.')
        ymd = dtime.split('-')[0]

        if ymd        !='%04d%02d%02d'%(Year,Mon,Day): continue
        if prdNameTmp !=prdName: continue
        if sateTmp    !=sate   : continue
        if sensorTmp  !=sensor  : continue
        if verFullTmp !='%s%s'%(ver,minorver): continue
        oPath = oDir + "/" + sPath.split("/")[-1]
        ftp.retrbinary("RETR %s"%(sPath), open(oPath,"wb").write)
        print oPath
     
ftp.close()

