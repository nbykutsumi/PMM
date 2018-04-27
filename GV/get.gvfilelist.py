import os
from ftplib import FTP
import myfunc.util as util

#fformat = "GMIN"
fformat = "2A56"
remotehost = "hector.gsfc.nasa.gov"
rootDir = "/gpmarchive/Gauge/%s"%(fformat)

user = "anonymous"

ftp = FTP(remotehost)
ftp.login(user)

lfilePath = []
lnworkDir = ftp.nlst(rootDir)
for nworkDir in lnworkDir:
    lsiteDir  = ftp.nlst(nworkDir)
    for siteDir in lsiteDir:
        lYearDir = ftp.nlst(siteDir)
        for YearDir in lYearDir:
            ltmpPath = ftp.nlst(YearDir)
            print ltmpPath
            if len(ltmpPath)==0: continue

            lfilePath = lfilePath + ltmpPath


sout = ("\n".join(lfilePath)).strip()
listDir  = "/work/a01/utsumi/data/GPMGV"
listPath = listDir + "/filelist.%s.txt"%(fformat)
f = open(listPath, "w"); f.write(sout); f.close()

print listPath



ftp.quit()
