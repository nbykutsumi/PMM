import subprocess

iyear = 2001
eyear = 2011
imon  = 1
emon  = 12
lyear = range(iyear,eyear+1)
lmon  = range(imon,emon+1)
print lyear
print lmon
print "********************"
prog = "./get.GPM.py"
for year in lyear:
  for mon in lmon:
    cmd = "%s %s %s"%(prog, year, mon)
    #if year==2002 and mon in [1,2,3,4,5,6]:
    #  continue
    subprocess.call(cmd, shell=True)
