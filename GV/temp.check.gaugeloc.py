import os, sys
import glob
import myfunc.util as util
from numpy import *

srcDir = "/work/a01/utsumi/data/GPMGV/sitelist"
lfilePath = glob.glob(srcDir + "/allfiles.*.csv")

for filePath in lfilePath:
    fileName = os.path.basename(filePath)
    regionName = fileName.split(".")[1]
    nwName     = fileName.split(".")[2]
    f=open(filePath,"r"); lines=f.readlines(); f.close()

    if len(lines)==0:
        print regionName,nwName,"No data"

    dlatlon = {}
    for line in lines:
        line = line.strip().split(",")
        gName= line[4]
        lat  = float(line[7])
        lon  = float(line[8])
        if gName in dlatlon.keys():
            if [lat,lon] not in dlatlon[gName]:
                dlatlon[regionName,nwName,gName].append([lat,lon])
        else:
            dlatlon[regionName,nwName,gName] = [[lat,lon]]

    for key in dlatlon.keys():
        regionName, nwName, gName = key
        alatlon = array(dlatlon[key])
        #print regionName, nwName, gName, dlatlon[gName]
        if len(alatlon) == 1: continue
        #if alatlon.std(axis=0).max() < 1: continue
        print regionName, nwName, gName, alatlon.mean(axis=0), alatlon.std(axis=0)



