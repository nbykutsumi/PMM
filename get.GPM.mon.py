import os
import socket
import pysftp
import myfunc.util as util

hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif hostname=="well":
    orootDir   = "/media/disk2/share/data/GPM"


iYM     = [2012,1]
eYM     = [2012,2]
lYM	= util.ret_lYM(iYM,eYM)


prj	= ["TRMM","TRMM.PR","L3A25","07"]
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
#print sftp.listdir()


for Year, Mon in lYM:
    print Year,Mon
    srcDir = os.path.join(ibaseDir, "%04d"%(Year), "%02d"%(Mon))
    outDir = os.path.join(obaseDir, "%04d"%(Year), "%02d"%(Mon))

    sftp.cwd(srcDir)
    lfileName = sftp.listdir()
    print "srcDir=",srcDir
    util.mk_dir(outDir)

    for fileName in lfileName:
        srcPath = os.path.join(srcDir,fileName)
        outPath = os.path.join(outDir,fileName)

        if os.path.exists(outPath) and ow==False:
            print "Skip (ow:%5s) ::"%(ow, outPath)

            continue

   
        with pysftp.cd(outDir):
            sftp.get(srcPath)
            print outDir
