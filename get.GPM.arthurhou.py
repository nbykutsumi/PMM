from ftplib import FTP
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket
from datetime import datetime, timedelta
import numpy as np

hostname  = "arthurhou.pps.eosdis.nasa.gov"
irootDir = "/gpmallversions"
#irootDir = "/sm/730/gpmdata"

myhost = socket.gethostname()
if  myhost =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif myhost =="shui":
    orootDir   = "/work/hk02/PMM/NASA"
    #orootDir   = "/tank/utsumi/data/PMM/NASA"
elif myhost =="well":
    #orootDir   = "/media/disk2/share/data/GPM"
    orootDir   = "/media/disk2/share/data/PMM/NASA"

#GPM/TRMM.TMI/L2A12/07/2014/


iDTime  = datetime(2015,5,1)
#iDTime  = datetime(2014,5,1)
#iDTime  = datetime(2017,1,1)
#iDTime  = datetime(2014,7,11)
#eDTime  = datetime(2014,7,11)
eDTime  = datetime(2015,5,31)
dDTime  = timedelta(days=1)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)

#spec      = ["TRMM","TMI","2A-CLIM","gprof","V05","A"]
#spec      = ["GPM","GMI","1C","1C","V05","A"]
#spec      = ["GPM","GMI","2A-CLIM","gprof","V05","A"] # input=ECMWF
#spec      = ["GPM","GMI","2A","gprof","V05","A"] # input=GANAL
spec      = ["GPM","Ku","2A","radar","V06","A"]
#spec      = ["GPM","Ka","2A","radar","V06","A"]
#spec      = ["GPM","DPR","2A","radar","V06","A"]
#spec      = ["GPM","DPRGMI","2B","radar","V06","A"]

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
    #--- Skip missing or problematic data -----
    ldataoutage=[
     [2014,5,20]  # 1274-1277
    ,[2014,10,22] # 3682-
    ,[2014,10,23] # all
    ,[2014,10,24] # -3718
    ,[2014,12,9]  # 004425, 004426(missing some lat&lon)
    ,[2014,12,10]  # 004447-004452
    ,[2017,9,26]
    ,[2017,9,27]
    ,[2017,9,28]
    ,[2017,9,29]
    ]


    if [Year,Mon,Day] in ldataoutage:
        print 'Skip data outage day'
        continue

    #if Mon==7: continue  # test

    #--- path and directory: Remote -----------
    iDir = irootDir + "/%s/%04d/%02d/%02d/%s"%(ver,Year,Mon,Day,prj)
    #oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d"%(sate,sensor,prdName,ver,Year,Mon) 
    oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d/%02d"%(sate,sensor,prdName,ver,Year,Mon,Day) 
    #GPM/TRMM.TMI/L2A12/07/2014/
    mk_dir(oDir)
    #--- list --------------
    lPath = np.sort(ftp.nlst(iDir))
    print iDir
    print lPath
    for sPath in lPath:
        oid = int(sPath.split('.')[-3])

        #if oid != 1924: continue # test

        fName = os.path.basename(sPath)
        prdNameTmp, sateTmp, sensorTmp, algFullTmp, dtime, gNum, verFullTmp, sfx = fName.split('.')
        ymd = dtime.split('-')[0]
        if fName.split('.')[3]=='GPM-SLH': continue
        if ymd        !='%04d%02d%02d'%(Year,Mon,Day): continue
        if prdNameTmp !=prdName: continue
        if sateTmp    !=sate   : continue
        if sensorTmp  !=sensor  : continue
        if verFullTmp !='%s%s'%(ver,minorver): continue
        oPath = oDir + "/" + sPath.split("/")[-1]
        ftp.retrbinary("RETR %s"%(sPath), open(oPath,"wb").write)
        print oPath
     
ftp.close()

