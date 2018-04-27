import os, sys
import myfunc.util as util

def load_sitelist(siteName=None):
    srcDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
    srcPath = srcDir +  "/%s_loc.dat"%(siteName)
    f = open(srcPath,"r"); lines=f.readlines(); f.close()

    dloc = {}

    if siteName in ["CEMADEN"]:
        for line in lines[1:]:
            line    = line.split()
            gauge   = line[0]
            lat     = float(line[1])
            lon     = float(line[2])
            dloc[gauge] = [lat, lon]





    else:
        for line in lines:
            line = line.split()

            print len(line),line
            if len(line)==0: break

            gauge = line[0]
            lon   = float(line[2]) + float(line[3])/60.0 + float(line[4])/3600. 
            lat   = float(line[5]) + float(line[6])/60.0 + float(line[7])/3600. 
            dloc[gauge] = [lat, lon]

    return dloc
    


    #if siteName in ["APU","BDF","CAL","CSC","ERK","FIN","GSFC","HAR","INP","KAM","KAP","KP2","KSC","KWA","LBA","NNN","OMK","RMI","SFL","SFW","SPG","STJ","TAM","TFB","WAL",



if __name__ == "__main__":
    siteName = "GSFC"
    #siteName = "CEMADEN"
    a= load_sitelist(siteName)
    print a.keys()
    print a.values()

