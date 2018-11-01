import subprocess, os
import socket
from glob import glob
#-------------------
hostname = socket.gethostname()
if   hostname == "mizu":
    srcDir =  "/home/utsumi/mnt/wellshare/data/GPM"
elif hostname == "well":
    srcDir =  "/media/disk2/share/data/GPM"

#-------------------
def ret_lpath(prj,year,mon):
  #idir  = os.path.join("/media/disk2/data",prj[0],prj[1],"%04d"%(year),"%02d"%(mon))
  #idir  = os.path.join("/data2/GPM",prj[0],prj[1],prj[2],"%04d"%(year),"%02d"%(mon))
  #idir  = os.path.join("/home/utsumi/mnt/wellshare/data/GPM",prj[0],prj[1],prj[2],"%04d"%(year),"%02d"%(mon))
  idir  = os.path.join(srcDir,prj[0],prj[1],prj[2],"%04d"%(year),"%02d"%(mon))
  lpath = glob( os.path.join(idir, "*.gz"))

  print idir
  return lpath

#prj          = ["TRMM.PR","L2A25","07"]
prj          = ["TRMM.PR","L2A23","07"]
#prj          = ["TRMM.PR","L2A25","07"]
#prj          = ["TRMM.PR","L3A25","07"]
#prj          = ["TRMM.TMI","L2A12"]
#iyear, eyear = 2009,2013
iyear, eyear = 2014,2014
lyear        = range(iyear,eyear+1)
imon,emon    = 1,12
lmon         = range(imon,emon+1)
#lmon         = [4,5,6,7,8,9,10]
#lmon         = [10,9,8,7,6,5,4]

for year in lyear[::-1]:
  for mon in lmon:
    lpath  = ret_lpath(prj,year,mon)
    print lpath
    for path in lpath:
      #print os.path.exists(path), path
      cmd = "gunzip -dv %s"%(path)
      subprocess.call(cmd, shell=True)



