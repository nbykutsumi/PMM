from numpy import *
import Image


baseDir = "/home/utsumi/mnt/well.share"
figDir  = baseDir + "/PMM/WNP.261x265/pict"

#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["KuPR","IMERG","IMERG.IR","IMERG.MW","GMI","GSMaP","GSMaP.IR","GSMaP.MW"]

lcltype = range(0,7+1)
#lcltype = [1]
ncltype = len(lcltype)
dclid   = {0:0, 1:1, 2:201, 3:202, 4:4, 5:3, 6:204, 7:200}
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy", 99:"All"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw", 99:"All"}


iy  = 108  # top
ey  = -1  # bottom
ix  = 1   # left
ex  = -108 # right

for icl in lcltype:

  da2dat = {}
  for i,dattype in enumerate(ldattype):
    figPath   = figDir  + "/scatter.%s.vsRA.%s.png"%(dattype,dclShortName[icl])
    iimg      = Image.open(figPath)
    a2array   = asarray(iimg)
    print shape(a2array)
    da2dat[i] = a2array[iy:ey, ix:ex]

  a2line1  = hstack([da2dat[0], da2dat[1], da2dat[2], da2dat[3]])
  a2line2  = hstack([da2dat[4], da2dat[5], da2dat[6], da2dat[7]])
  a2oarray = vstack([a2line1, a2line2])
  oimg     = Image.fromarray(a2oarray)

  oPath    = figDir + "/join.pdf.%s.png"%(dclShortName[icl])
  oimg.save(oPath)
  print oPath


