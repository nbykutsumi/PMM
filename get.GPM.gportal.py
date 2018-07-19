import os
import socket
import pysftp
import myfunc.util as util

hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif hostname=="well":
    orootDir   = "/media/disk2/share/data/GPM"


#iYM     = [2014,4]
#eYM     = [2015,6]
iYM     = [2014,10]
eYM     = [2014,10]

lYM	= util.ret_lYM(iYM,eYM)
lYM = [YM for YM in lYM if not YM[1] in [11,12,1,2,3]]
lYM = lYM[::-1]

ow  = False # OverWrite
#ow  = True # OverWrite
#prj	= ["TRMM","TRMM.PR","L3A25","07"]
#prj	= ["TRMM","TRMM.PR","L2A23","07"]
#prj	= ["TRMM","TRMM.PR","L2A25","07"]
#prj	= ["TRMM","TRMM.TMI","L2A12","07"]
prj	= ["GPM","GPM.KuPR","L2","05"]
#prj	= ["GPM","GPM.DPR","L2","05"]
#prj	= ["GPM","GPM.GMI","L2","05"]
sate	= prj[0]
sensor	= prj[1]
prdName = prj[2]
version = prj[3]

host        = 'sftp.gportal.jaxa.jp'
username    = 'nbyk.utsumi'
port        = 2051
ibaseDir    = os.path.join("/" 'standard', sate, sensor, prdName, version)
obaseDir    = os.path.join(orootDir, sensor, prdName, version)

sftp	= pysftp.Connection(host=host,username=username, port=port)
sftp.cwd(ibaseDir)
print ibaseDir
print sftp.listdir()


for Year, Mon in lYM:
    print Year,Mon
    srcDir = os.path.join(ibaseDir, "%04d"%(Year), "%02d"%(Mon))
    outDir = os.path.join(obaseDir, "%04d"%(Year), "%02d"%(Mon))

    try:
        sftp.cwd(srcDir)
    except IOError:
        print "No remote Directory: %s"%(srcDir)

    lfileName = sftp.listdir()
    print "srcDir=",srcDir
    util.mk_dir(outDir)

    for fileName in lfileName:
        srcPath = os.path.join(srcDir,fileName)
        outPath = os.path.join(outDir,fileName)

        if os.path.exists(outPath) and ow==False:
            print "Skip (ow:%5s) :: %s"%(ow, outPath)

            continue

   
        with pysftp.cd(outDir):
            try:
                sftp.get(srcPath)
            except IOError:
                print "No remote File: %s"%(srcPath)
                continue
            print outPath



