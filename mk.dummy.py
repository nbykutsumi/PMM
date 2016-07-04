from numpy import *

Y,M  = 2001,8
idir = "/media/disk2/data/TRMM/PR.2A25.sa.one.FWD/%04d/%02d"%(Y,M)
for D in range(7,22+1):
  for H in range(0,23+1):
    #-------------
    if D==7 and H in [0,1,2,3,4,5,6,7,8,9,10,11,12,13]:
      continue
    #-------------
    a2dat = ones([80,360],float32)*(-9999.)
    sname = idir + "/PR.%04d.%02d.%02d.%02d.sa.one"%(Y,M,D,H)
    #a2dat.tofile(sname)
    #print sname
    

