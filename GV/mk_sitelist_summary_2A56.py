import myfunc.util as util
from numpy import *

srcDir  =  "/work/a01/utsumi/data/GPMGV/sitelist"
srcPath = srcDir + "/sitelist_reclassified.csv"

f=open(srcPath,"r"); lines=f.readlines(); f.close()

Y0, Y1 = 1995,2020

#---------------------------
def ret_kYM(YM,Y0):
    Year = YM[0]
    Mon  = YM[1]
    return Year - Y0 + Mon -1

def YYYYMM2YM(s):
    Year = int(s[:4])
    Mon  = int(s[4:])
    return [Year,Mon]

def YM2YYYYMM(YM):
    return "%04d%02d"%(YM[0],YM[1])


#---------------------------


lkey = []
dYM  = {}
dnwCode={}
dlat = {}
dlon = {}
for line in lines[1:]:
    line = line.strip().split(",")
    region  = line[0]
    nwName  = line[1]
    nwCode  = line[2]
    domain  = line[3]
    gCode   = line[4]
    lat     = float(line[5])
    lon     = float(line[6])

    sYear   = int(line[7])
    sMon    = int(line[8])
    eYear   = int(line[9])
    eMon    = int(line[10])

    key = domain

    lYM     = util.ret_lYM([sYear,sMon],[eYear,eMon])
    if len(line)>11:
        lnoYM   = map(YYYYMM2YM, line[12:])
        ltmp    = [YM for YM in lYM if YM not in lnoYM ]
        lYM     = ltmp

    # initialize
    if key not in lkey:
        lkey.append(key)
        dYM[key] = lYM
        dnwCode[key] = nwCode
        dlat[key]= [lat]
        dlon[key]= [lon]


    else:
        dlat[key].append(lat)
        dlon[key].append(lon)

        for YM in lYM:
            if YM not in dYM[key]:
                dYM[key].append(YM)

#print dlat[("BRAZIL","INP")]


#---
lYM = util.ret_lYM([Y0,1],[Y1,12])
lk  = arange(len(lYM))
sout = ""
for key in lkey:
    domain = key
    region, nwName = key.split("-")[:2]
    nwCode = dnwCode[key]
    dYM[key] = sorted(dYM[key])
    sYM = min(dYM[key])
    eYM = max(dYM[key])

    lllat= min(dlat[key])
    urlat= max(dlat[key])
    lllon= min(dlon[key])
    urlon= max(dlon[key])

    ltmp = [region,nwName,nwCode,domain,lllat,urlat,lllon,urlon,sYM[0],sYM[1],eYM[0],eYM[1]]
    lyyyymm = map(YM2YYYYMM, dYM[key])
    ltmp = ltmp + lyyyymm
    ltmp = map(str, ltmp)
    sout = sout + ",".join(ltmp) + "\n"


slabel = ",".join(["region","nwName","nwCode","domain","lllat","urlat","lllon","urlon","sYear","sMon","eYear","eMon","ObsYYYYMM"])

sout    = slabel + "\n" + sout.strip()
outPath = srcDir + "/sitelist_summary.csv"
f = open(outPath, "w"); f.write(sout); f.close()
print outPath
