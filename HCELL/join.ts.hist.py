from numpy import *
import Image, shutil, socket
import os

hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

#lvar = ["Count","PDF"]
lvar = ["Count"]
lNorm= [True,False]

baseDir = os.path.join(rootDir, "HCELL/PDF/Type")


lielat= [(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35)]  # degree
lielat= lielat[::-1]

lseason = ["DJF","MAM","JJA","SON"]
lregion = ["NAT","NAF","ASI","NPA","NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]
#lregion = ["SAT","SAF","SIN","OCE","SWP","SEP","SAM"]

dhemi = {}
for region in ["NAT","NAF","ASI","NPA","NAM"]:
    dhemi[region] = "N"
for region in ["SAT","SAF","SIN","OCE","SWP","SEP","SAM"]:
    dhemi[region] = "S"


rtype = "all"

lkey = [[var,Norm] for var in lvar for Norm in lNorm]
for [var,Norm] in lkey:
    for region in lregion:
        for landsea in ["lnd","sea"]:
            i = -1
            da2dat = {}
            for season in lseason:

                if dhemi[region]=="N":
                    lielat_tmp=lielat
                elif dhemi[region]=="S":
                    lielat_tmp=lielat[::-1]

                for ielat in lielat_tmp:
                    ilat, elat = ielat
                    i = i+1 
                    figDir = "/tank/utsumi/PMM/HCELL/pict"
                    if Norm==False:
                        figPath= os.path.join(figDir,"TS.%s.%s.%s.%s.%s.Lat.%.1f.%.1f.png"%(var,landsea,region,rtype,season,ilat,elat))
                    elif Norm==True:
                        figPath= os.path.join(figDir,"TS.%s.Norm.%s.%s.%s.%s.Lat.%.1f.%.1f.png"%(var,landsea,region,rtype,season,ilat,elat))
        
                    a2png  = Image.open(figPath)
                    a2array= asarray(a2png)
                    da2dat[i] = a2array
        
            a2col1 =vstack([da2dat[0],da2dat[1],da2dat[2],da2dat[3],da2dat[4],da2dat[5],da2dat[6]])    
            a2col2 =vstack([da2dat[7],da2dat[8],da2dat[9],da2dat[10],da2dat[11],da2dat[12],da2dat[13]])    
            a2col3 =vstack([da2dat[14],da2dat[15],da2dat[16],da2dat[17],da2dat[18],da2dat[19],da2dat[20]])    
            a2col4 =vstack([da2dat[21],da2dat[22],da2dat[23],da2dat[24],da2dat[25],da2dat[26],da2dat[27]])    
        
            a2out  = hstack([a2col1,a2col2,a2col3,a2col4])
        #    a2out  = hstack([a2col1,a2col2,a2col3])
            oimg   = Image.fromarray(uint8(a2out))
            
            if Norm==False:
                oPath  = os.path.join(figDir,"join.TS.%s.%s.%s.%s.png"%(var,landsea,region,rtype))
            elif Norm==True:
                oPath  = os.path.join(figDir,"join.TS.%s.Norm.%s.%s.%s.png"%(var,region,landsea,rtype))
        
            oimg.save(oPath)
            print oPath
    
    

