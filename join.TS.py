from numpy import *
import Image, shutil

nclass = 4
ltag    = ["plain","c","fbc","tc","ot"]
lvtype  = ["Bias.Prc","Bias.Rte","R.Bias.Rte","TR.Bias.Prc","Corr"]
#lvtype  = ["Prcp.RA"]
figDir  = "/tank/utsumi/PMM/RAdom/pict/%dclass"%(nclass)

for vtype in lvtype:
  da2dat  = {}
  for i,tag in enumerate(ltag):
    figPath   = figDir + "/TS.%s.%s.png"%(vtype, tag)
    a2png     = Image.open(figPath)
    a2array   = asarray(a2png)
    da2dat[i] = a2array
  
  a2dummy  = ones(a2array.shape)*255
  a2line1  = hstack([a2dummy,   da2dat[0]])
  a2line2  = hstack([da2dat[1], da2dat[2]])
  a2line3  = hstack([da2dat[3], da2dat[4]])
  a2oarray = vstack([a2line1, a2line2, a2line3])
  oimg     = Image.fromarray(uint8(a2oarray))
  
  oPath    = figDir + "/join.TS.%s.png"%(vtype)
  oimg.save(oPath)
  print oPath
