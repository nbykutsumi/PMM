import os, sys
import glob
import subprocess

#fformat = "GMIN"
fformat = "2A56"
rootDir = "/work/a01/utsumi/data/GPMGV/%s"%(fformat)

lnwDir = glob.glob(rootDir+"/*")
#lsiteName = ["SFL"]
#lsiteName = ["GSFC"]
lsiteName = []

for nwDir in lnwDir:
    lsiteDir = glob.glob(nwDir + "/*")
    for siteDir in lsiteDir:
        lYearDir = glob.glob(siteDir + "/*")
        for YearDir in lYearDir:
            os.chdir(YearDir)
            ldatPath = glob.glob(YearDir + "/*")
            for datPath in ldatPath:
                datName = datPath.split("/")[-1]
                siteName = datName.split("_")[1]
                if len(lsiteName)==0:
                    pass
                else:
                    if siteName not in lsiteName:
                        continue 

                if   datPath.split(".")[-1] in ["tgz"]:
                    print datPath
                    cmd = ["tar","-xzvf",datPath]
                    subprocess.call(cmd)
                elif datPath[-6:] == "tar.gz":
                    print datPath
                    cmd = ["tar","-xzvf",datPath]
                    subprocess.call(cmd)

                elif datPath.split(".")[-1] in ["tar"]:
                    print datPath
                    cmd = ["tar","-xvf",datPath]
                    print cmd
                    subprocess.call(cmd)
                else:
                    continue 

               


