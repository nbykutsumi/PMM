import os
import socket
#import pysftp
from ftplib import FTP
import myfunc.util as util
from datetime import datetime, timedelta


hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif hostname =="shui":
    orootDir   = "/work/hk02/PMM/NASA"
    #orootDir   = "/tank/utsumi/data/PMM/NASA"
elif hostname=="well":
    #orootDir   = "/media/disk2/share/data/GPM"
    orootDir   = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA'


lYear = range(2012,2015)

#ow  = False # OverWrite
ow  = True # OverWrite
prj	= ['TRMM_GPMFormat','TRMM.PR','3A.PR.M','06A']
#prj	= ["TRMM","TRMM.PR","L2A23","07"]
#prj	= ["TRMM","TRMM.PR","L2A25","07"]
#prj	= ["TRMM","TRMM.TMI","L2A12","07"]
#prj	= ["GPM","GPM.KuPR","L2","05"]
#prj	= ["GPM","GPM.DPR","2A.DPR","05A"]
#prj	= ["GPM","GPM.GMI","L2","05"]
sate	= prj[0]
sensor	= prj[1]
prdName = prj[2]
version = prj[3]

host        = 'ftp.gportal.jaxa.jp'
username    = 'nbyk.utsumi'
passwd      = 'anonymous'
port        = 21

ftp	= FTP(host=host,user=username, passwd=passwd)

for Year in lYear:
    '''/standard/TRMM_GPMFormat/TRMM.PR/3A.PR.D/06A/2014'''
    ibaseDir    = os.path.join("/" 'standard', sate, sensor, prdName, version)
    obaseDir    = os.path.join(orootDir, sensor, prdName, version)
    #obaseDir    = orootDir + "/%s/%s/%s/%04d/%02d/%02d"%(sate,sensor,prdName,ver,Year,Mon,Day)

    srcDir = os.path.join(ibaseDir, "%04d"%(Year))
    outDir = os.path.join(obaseDir, "%04d"%(Year))

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



