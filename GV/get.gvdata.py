import os
from ftplib import FTP
import myfunc.util as util

#fformat = "GMIN"
fformat = "2A56"
remotehost = "hector.gsfc.nasa.gov"
rootDir = "/gpmarchive/Gauge/%s"%(fformat)

user = "anonymous"

obaseDir  = "/work/a01/utsumi/data/GPMGV/%s"%(fformat)
util.mk_dir(obaseDir)

ftp = FTP(remotehost)
ftp.login(user)

lnworkDir = ftp.nlst(rootDir)
for nworkDir in lnworkDir:
    lsiteDir  = ftp.nlst(nworkDir)
    for siteDir in lsiteDir:
        lYearDir = ftp.nlst(siteDir)
        for YearDir in lYearDir:
            nwName, siteName, Year = YearDir.split("/")[-3:]
            if nwName in ["BRAZIL",  "CALIFORNIA",  "DARWIN",  "FLORIDA",  "FRANCE"]: continue
            outDir = obaseDir + "/%s"%(nwName)
            print outDir
            util.mk_dir(outDir)
            outDir = obaseDir + "/%s/%s"%(nwName,siteName)
            print outDir
            util.mk_dir(outDir)
            outDir = obaseDir + "/%s/%s/%s"%(nwName, siteName, Year)
            print outDir
            util.mk_dir(outDir)

            # change shell directory
            os.chdir(outDir)
            #------------------

            print nwName, siteName, Year
            ftp.cwd(YearDir)
            lfileName = ftp.nlst(".")
            print lfileName
            if len(lfileName)==0: continue
            for fileName in lfileName:
                with open(fileName, "w") as f:
                    ftp.retrbinary("RETR %s"%(fileName), f.write)

            


ftp.quit()
