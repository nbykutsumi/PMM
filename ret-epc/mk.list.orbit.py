import os
import socket
import myfunc.util as util
import random
import glob

iYM = [2014,6]
#eYM = [2014,6]
eYM = [2015,5]
lYM = util.ret_lYM(iYM,eYM)
#nsample = 1000
nsample = 2

sensor  = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    workbaseDir = '/work'
    tankbaseDir = '/tank'
elif myhost =="well":
    workbaseDir = '/home/utsumi/mnt/lab_work'
    tankbaseDir = '/home/utsumi/mnt/lab_tank'

#*******************
# Function
#*******************
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass
#*******************
srcbaseDir = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05'
lsrcPath = []
for Year,Mon in lYM:
    ssearch = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/??/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon)
    lsrcPathTmp = sorted(glob.glob(ssearch))
    lsrcPath = lsrcPath + lsrcPathTmp


random.seed(0)
lsrcPath = sorted(random.sample(lsrcPath, nsample))
sout = '\n'.join(lsrcPath)

#-- Save -----
listDir = tankbaseDir + '/utsumi/PMM/retepc/list-orbit'
util.mk_dir(listDir)
outPath = listDir + '/list.%s.%s.%04d%02d-%04d%02d.%dobts.txt'%(sensor,myhost,iYM[0],iYM[1],eYM[0],eYM[1],nsample)
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
