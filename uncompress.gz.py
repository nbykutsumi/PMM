import subprocess, os
from glob import glob
#-------------------
def ret_lpath(prj,year,mon):
  if prj == ["TRMM.PR","L2A25"] or prj == ["TRMM.TMI","L2A12"] or prj == ["TRMM.PR","L3A25"]:
    idir  = os.path.join("/media/disk2/data",prj[0],prj[1],"%04d"%(year),"%02d"%(mon))
    lpath = glob( os.path.join(idir, "*.gz"))

  print idir
  return lpath

#prj          = ["TRMM.PR","L2A25"]
prj          = ["TRMM.PR","L3A25"]
#prj          = ["TRMM.TMI","L2A12"]
iyear, eyear = 2001,2010
lyear        = range(iyear,eyear+1)
imon,emon    = 1,12
lmon         = range(imon,emon+1)

for year in lyear:
  for mon in lmon:
    lpath  = ret_lpath(prj,year,mon)
    for path in lpath:
      cmd = "gunzip -dv %s"%(path)
      subprocess.call(cmd, shell=True)



