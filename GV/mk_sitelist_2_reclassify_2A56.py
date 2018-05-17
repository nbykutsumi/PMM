import os
import myfunc.util as util

srcDir  =  "/work/a01/utsumi/data/GPMGV/sitelist"
srcPath = srcDir + "/sitelist.csv"

f= open(srcPath,"r"); lines = f.readlines(); f.close()

sout = ""
for line in lines[1:]:
    line = line.strip().split(",")
    region  = line[0]
    nwName  = line[1]
    nwCode  = line[2]
    gCode   = line[3]
    lat     = float(line[4])
    lon     = float(line[5])

    sYear   = int(line[6])
    sMon    = int(line[7])
    eYear   = int(line[8])
    eMon    = int(line[9])


    regCode = "%s-%s"%(region,nwName)

    key = (region,nwName)
    if key == ("VIRGINIA","NASA"):
        if lon < -89.775:
            regCode = "%s-%s-W"%(region,nwName)
        elif lat < 36.315:
            regCode = "%s-%s-C"%(region,nwName)
        elif (36.315 <lat) & (lat < 37.68):
            regCode = "%s-%s-NE"%(region,nwName)
        else: 
            regCode = "%s-%s-SE"%(region,nwName)

    if key == ("FLORIDA","KAM"):
        if lon < -81.252:
            regCode = "%s-%s-W"%(region,nwName)
        else:
            regCode = "%s-%s-E"%(region,nwName)

    if key == ("FLORIDA","SFL"):
        if lat < 24.81:
            regCode = "%s-%s-S"%(region,nwName)
        else:
            regCode = "%s-%s-N"%(region,nwName)

    if key == ("FRANCE","HyMeX"):
        if lon < 5.85: 
            regCode = "%s-%s-W"%(region,nwName)
        else:
            regCode = "%s-%s-E"%(region,nwName)

    if key == ("MARYLAND","PCMK"):
        if lat < 37.584: 
            regCode = "%s-%s-S"%(region,nwName)
        else:
            regCode = "%s-%s-N"%(region,nwName)

    if key == ("VIRGINIA","NSWD"):
        if lat < 37.585: 
            regCode = "%s-%s-S"%(region,nwName)
        else:
            regCode = "%s-%s-N"%(region,nwName)

 

    ltmp = line[:3] + [regCode] + line[3:]
    stmp = ",".join(ltmp)
    sout = sout + stmp + "\n"


            
outPath = srcDir + "/sitelist_reclassified.csv"
f = open(outPath, "w"); f.write(sout); f.close()
print outPath




