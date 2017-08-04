import subprocess, os
from glob import glob
#-------------------
def ret_lpath(prj,year,mon):
  if prj[0] == "TRMM.PR":
    #idir  = os.path.join("/media/disk2/data",prj[0],prj[1],"%04d"%(year),"%02d"%(mon))
    #idir  = os.path.join("/data2/GPM",prj[0],prj[1],prj[2],"%04d"%(year),"%02d"%(mon))
    idir  = os.path.join("/home/utsumi/mnt/wellshare/data/GPM",prj[0],prj[1],prj[2],"%04d"%(year),"%02d"%(mon))
    lpath = glob( os.path.join(idir, "*.gz"))

  print idir
  return lpath

#prj          = ["TRMM.PR","L2A25","07"]
prj          = ["TRMM.PR","L2A25","07"]
#prj          = ["TRMM.PR","L3A25"]
#prj          = ["TRMM.TMI","L2A12"]
iyear, eyear = 2001,2001
lyear        = range(iyear,eyear+1)
imon,emon    = 1,1
lmon         = range(imon,emon+1)

for year in lyear:
  for mon in lmon:
    lpath  = ret_lpath(prj,year,mon)
    print lpath
    for path in lpath:
      cmd = "gunzip -dv %s"%(path)
      subprocess.call(cmd, shell=True)



