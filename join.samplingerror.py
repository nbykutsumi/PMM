from numpy import *
import Image, shutil

nclass = 4
prtype = "RA"
lvtype  = ["Prcp","SamplingError"]
figDir  = "/tank/utsumi/PMM/RAdom/pict/%dclass"%(nclass)
ltag   = ["plain","c","fbc","tc","ot"]
da2dat  = {}
for itype, vtype in enumerate(lvtype):
  for itag,tag in enumerate(ltag):
    figPath   = figDir + "/TS.%s.%s.%s.png"%(vtype, prtype, tag)
    a2png     = Image.open(figPath)
    a2array   = asarray(a2png)
    da2dat[itype,itag] = a2array

da2comb = {} 
for itag, tag in enumerate(ltag):
  da2comb[itag] = vstack([da2dat[0,itag], da2dat[1,itag]]) 

a2dummy  = ones(da2comb[0].shape)*255
a2line1  = hstack([a2dummy,    da2comb[0]])
a2line2  = hstack([da2comb[1], da2comb[2]])
a2line3  = hstack([da2comb[3], da2comb[4]])
a2oarray = vstack([a2line1, a2line2, a2line3])
oimg     = Image.fromarray(uint8(a2oarray))

oPath    = figDir + "/join.TS.SamplingError.%s.png"%(prtype)
oimg.save(oPath)
print oPath
