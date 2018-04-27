import os
from ftplib import FTP
import myfunc.util as util

fformat = "GMIN"
remotehost = "hector.gsfc.nasa.gov"
rootDir = "/gpmarchive/Gauge/Sitelist"

user = "anonymous"

outDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
util.mk_dir(outDir)
os.chdir(outDir)

ftp = FTP(remotehost)
ftp.login(user)

lfileName = []
ftp.cwd(rootDir)
lfileName = ftp.nlst(".")
for fileName in lfileName:
    print fileName
    with open(fileName, "w") as f:
        ftp.retrbinary("RETR %s"%(fileName), f.write)

        outPath = outDir + "/"+ fileName
        print outPath



ftp.quit()
