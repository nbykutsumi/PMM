import os
import socket
#import pysftp
from ftplib import FTP
import myfunc.util as util
from datetime import datetime, timedelta


hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif hostname=="well":
    orootDir   = "/media/disk2/share/data/GPM"


iDTime = datetime(2017,11,20)
eDTime = datetime(2017,12,20)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

#ow  = False # OverWrite
ow  = True # OverWrite
#prj	= ["TRMM","TRMM.PR","L3A25","07"]
#prj	= ["TRMM","TRMM.PR","L2A23","07"]
#prj	= ["TRMM","TRMM.PR","L2A25","07"]
#prj	= ["TRMM","TRMM.TMI","L2A12","07"]
#prj	= ["GPM","GPM.KuPR","L2","05"]
prj	= ["GPM","GPM.DPR","2A.DPR","05A"]
#prj	= ["GPM","GPM.GMI","L2","05"]
sate	= prj[0]
sensor	= prj[1]
prdName = prj[2]
version = prj[3]

host        = 'ftp.gportal.jaxa.jp'
username    = 'nbyk.utsumi'
passwd      = 'anonymous'
port        = 21
ibaseDir    = os.path.join("/" 'standard', sate, sensor, prdName, version)
obaseDir    = os.path.join(orootDir, sensor, prdName, version)

ftp	= FTP(host=host,user=username, passwd=passwd)
ftp.cwd(ibaseDir)
print ibaseDir
print ftp.nlst()


for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    srcDir = os.path.join(ibaseDir, "%04d"%(Year), "%02d"%(Mon), "%02d"%(Day))
    outDir = os.path.join(obaseDir, "%04d"%(Year), "%02d"%(Mon), "%02d"%(Day))

    try:
        ftp.cwd(srcDir)
    except IOError:
        print "No remote Directory: %s"%(srcDir)

    lfileName = ftp.nlst()
    print "srcDir=",srcDir
    util.mk_dir(outDir)

    for fileName in lfileName:
        srcPath = os.path.join(srcDir,fileName)
        outPath = os.path.join(outDir,fileName)

        if os.path.exists(outPath) and ow==False:
            print "Skip (ow:%5s) :: %s"%(ow, outPath)

            continue

        with open(outPath, 'wb') as f:
            try:
                ftp.retrbinary('RETR %s'%(fileName), f.write)
            except IOError:
                print "No remote File: %s"%(srcPath)
                continue
        print outPath



