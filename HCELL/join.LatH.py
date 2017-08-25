from numpy import *
import Image, shutil, socket
import os

hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

#var  = "Count"
var  = "PDF"
Norm = True
#Norm = False

baseDir = os.path.join(rootDir, "HCELL/PDF/Type")


lseason = ["DJF","MAM","JJA","SON"]
#lseason = ["DJF","MAM","JJA"]

for landsea in ["lnd","sea"]:
    i = -1
    da2dat = {}
    for season in lseason:
        i = i+1 
        figDir = "/tank/utsumi/PMM/HCELL/pict"
        figPath= os.path.join(figDir,"plot.LatH.%s.%s.png"%(landsea,season))

        a2png  = Image.open(figPath)
        a2array= asarray(a2png)
        da2dat[i] = a2array

    a2line1 =hstack([da2dat[0],da2dat[1]])
    a2line2 =hstack([da2dat[2],da2dat[3]]) 

    a2out  = vstack([a2line1,a2line2])
    oimg   = Image.fromarray(uint8(a2out))
    
    oPath  = os.path.join(figDir,"join.LatH.%s.png"%(landsea))

    oimg.save(oPath)
    print oPath



