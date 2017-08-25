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


lielat= [(0,5),(5,10),(10,15),(15,20),(25,30),(30,35)]  # degree
lielat= lielat[::-1]

lseason = ["DJF","MAM","JJA","SON"]
#lseason = ["DJF","MAM","JJA"]

rtype = "all"
for landsea in ["land","sea"]:
    i = -1
    da2dat = {}
    for season in lseason:
        for ielat in lielat:
            ilat, elat = ielat
            i = i+1 
            figDir = "/tank/utsumi/PMM/HCELL/pict"
            if Norm==False:
                figPath= os.path.join(figDir,"TS.%s.%s.%s.%s.Lat.%.1f.%.1f.png"%(var,landsea,rtype,season,ilat,elat))
            elif Norm==True:
                figPath= os.path.join(figDir,"TS.%s.Norm.%s.%s.%s.Lat.%.1f.%.1f.png"%(var,landsea,rtype,season,ilat,elat))

            a2png  = Image.open(figPath)
            a2array= asarray(a2png)
            da2dat[i] = a2array

    a2col1 =vstack([da2dat[0],da2dat[1],da2dat[2],da2dat[3],da2dat[4],da2dat[5]])    
    a2col2 =vstack([da2dat[6],da2dat[7],da2dat[8],da2dat[9],da2dat[10],da2dat[11]])    
    a2col3 =vstack([da2dat[12],da2dat[13],da2dat[14],da2dat[15],da2dat[16],da2dat[17]])    
    a2col4 =vstack([da2dat[18],da2dat[19],da2dat[20],da2dat[21],da2dat[22],da2dat[23]])    

    a2out  = hstack([a2col1,a2col2,a2col3,a2col4])
#    a2out  = hstack([a2col1,a2col2,a2col3])
    oimg   = Image.fromarray(uint8(a2out))
    
    if Norm==False:
        oPath  = os.path.join(figDir,"join.TS.%s.%s.%s.png"%(var,landsea,rtype))
    elif Norm==True:
        oPath  = os.path.join(figDir,"join.TS.%s.Norm.%s.%s.png"%(var,landsea,rtype))

    oimg.save(oPath)
    print oPath



