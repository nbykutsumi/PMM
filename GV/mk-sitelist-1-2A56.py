import os, sys
import glob
from collections import deque
import myfunc.util as util

fformat = "2A56"
#rootDir = "/work/a01/utsumi/data/GPMGV/2A56/BRAZIL/LBA/1998/2A56_LBA_1998.tar"

rootDir = "/work/a01/utsumi/data/GPMGV/%s"%(fformat)

lregionPath = glob.glob(rootDir + "/*")

allFlag = True
#allFlag = False
summaryFlag = True
#summaryFlag = False

for regionPath in lregionPath:
    if allFlag !=True:
        continue

    regionName= regionPath.split("/")[-1]
    lnwPath = glob.glob(regionPath + "/*")

    #for nwPath in ['/work/a01/utsumi/data/GPMGV/2A56/VIRGINIA/HFD']:
    for nwPath in lnwPath:
        nwName= nwPath.split("/")[-1]
        #if not nwName == "IPHEx_NASA": continue
        lYearPath = glob.glob(nwPath + "/*")
        lYearPath = sorted(lYearPath)

        lout = deque([])
        for YearPath in lYearPath:
            Year     = YearPath.split("/")[-1]
            lsrcPath = glob.glob(YearPath + "/*")
            lsrcPath = sorted(lsrcPath)

            print regionName, nwName, Year
            for srcPath in lsrcPath:
                fileName = srcPath.split("/")[-1]
                sfx      = fileName.split(".")[-1]
                if   sfx=="asc":
                    YearMon = fileName.split("_")[4]
                    Year    = YearMon[0:4]
                    Mon     = YearMon[4:6]
                    f=open(srcPath,"r"); line=f.readline(); f.close()
                    line     = line.split()
                    nwCode   = line[1] + "_" + line[2]
                    gaugeID  = line[3]
                    lat  = line[7]
                    lon  = line[8]

                    #print sfx,regionName, nwName, nwCode, gaugeID, lat, lon, Year, Mon

                elif sfx=="2a56":
                    fName = os.path.basename(srcPath)
                    gaugeID = fName.split('-')[1]
                    f=open(srcPath,"r")
                    line=f.readline()
                    line     = line.split()
                    nwCode   = line[1]
                    #gaugeID  = line[2]
                    line = f.readline()   # read labels
                    line = f.readline()   # read first data line                    
                    line = line.split()
                    Year = line[0]
                    Mon  = line[1]
                    lat  = line[7]
                    lon  = line[8]

                    f.close()

                    #print sfx,regionName, nwName, nwCode, gaugeID, lat, lon, Year, Mon

                else: continue
                #----
                ltmp = [fileName,regionName,nwName,nwCode,gaugeID,Year,Mon,lat,lon]
                stmp = ",".join(ltmp)
                lout.append(stmp)


        #--- output ---
        lout    = list(lout)
        sout    = "\n".join(lout).strip()
        outDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
        outPath = outDir + "/allfiles.%s.%s.csv"%(regionName,nwName)
        f = open(outPath, "w"); f.write(sout); f.close()
        print outPath


#--- make summary -----
srcDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
lsrcPath= glob.glob(srcDir + "/allfiles.*.csv")

dlatlon = {}
dnwCode = {}
dsfx    = {}
dYM     = {}

for srcPath in lsrcPath:
    f=open(srcPath, "r"); lines=f.readlines(); f.close()
    print srcPath
    for line in lines:
        line     = line.strip().split(",")
        fileName = fileName = line[0]
        regionName  = line[1]
        nwName      = line[2]
        nwCode      = line[3]
        gName       = line[4]
        Year        = int(line[5])
        Mon         = int(line[6])
        Lat         = line[7]
        Lon         = line[8]
        sfx         = fileName.split(".")[-1]

        if (regionName, nwName, gName) in dlatlon.keys():
            dYM    [regionName, nwName, gName].append([Year,Mon])
        else:
            dlatlon[regionName, nwName, gName] = [Lat,Lon]
            dsfx   [regionName, nwName, gName] = sfx
            dnwCode[regionName, nwName, gName] = nwCode
            dYM    [regionName, nwName, gName] = [[Year,Mon]]
        

lkey = dlatlon.keys()
lkey = sorted(lkey)

lout = []
for key in lkey:

    #if key != ('BRAZIL', 'LBA', '0001'): continue

    lNoYM = []
    iYM = min(dYM[key])
    eYM = max(dYM[key])
    lYM = util.ret_lYM(iYM,eYM)
    for YM in lYM:
        if not YM in dYM[key]:
            s= "%04d%02d"%(YM[0],YM[1])
            lNoYM.append(s)

    regionName, nwName, gName = key
    nwCode = dnwCode[key]
    iYear, iMon = map(str, iYM)
    eYear, eMon = map(str, eYM)
    lat, lon    = dlatlon[key]

    ltmp = [regionName, nwName, nwCode, gName, lat, lon, iYear, iMon, eYear, eMon]
    if len(lNoYM) !=0:
        ltmp = ltmp + ["%s"%len(lNoYM)] + lNoYM
    else:
        ltmp = ltmp + ["0"]


    stmp = ",".join(ltmp)
    lout.append(stmp)


slabel= ",".join(["regionName","nwName","nwCode","gName","lat","lon","iYear","iMon","eYear","eMon","NoYM"])

sout = "\n".join(lout)
sout = slabel + "\n" + sout

outPath = srcDir + "/sitelist.csv"
f = open(outPath, "w")
f.write(sout)
f.close()
print outPath


 
