from numpy import *
import Image, shutil, socket
import os

#var  = "Count"
var  = "PDF"
#Norm = True
Norm = False


lielat= [(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35)]  # degree
lielat= lielat[::-1]

lseason = ["DJF","MAM","JJA","SON"]
#lseason = ["DJF","MAM","JJA"]
lhemi   = ["N","S"]

rtype = "all"
for lndsea in ["lnd","sea"]:
    i = -1
    da2dat = {}
    for hemi in lhemi:
        for season in lseason:
            i = i+1
            figDir = "/tank/utsumi/PMM/HCELL/pict"
            figPath = os.path.join(figDir, "plot.ts.stormH.%s.%s.%s.%s.png"%(hemi,season,"all",lndsea))
    
            a2png  = Image.open(figPath)
            a2array= asarray(a2png)
            da2dat[i] = a2array
    
    a2row1 =hstack([da2dat[0],da2dat[1],da2dat[2],da2dat[3]])    
    a2row2 =hstack([da2dat[4],da2dat[5],da2dat[6],da2dat[7]])
    
    a2out  = vstack([a2row1,a2row2])
    oimg   = Image.fromarray(uint8(a2out))
    
    oPath  = os.path.join(figDir,"join.plot.ts.stormH.%s.png"%(lndsea))
    
    oimg.save(oPath)
    print oPath



