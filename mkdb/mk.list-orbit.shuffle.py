import numpy as np
import myfunc.util as util
import random

prod       = '1C'
ver        = 'V05'
gmirootDir = '/work/hk01/PMM/NASA/GPM.GMI/%s/%s'%(prod,ver)

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)


listDir = '/work/hk01/utsumi/PMM/EPCDB/list'
util.mk_dir(listDir)
#-- Functions ---
def list2csv(a):
  if type(a[0]) !=list:
    sout = "\n".join(map(str,a))
  elif type(a[0]) ==list:
    lline = [",".join( map(str,line)) for line in a]
    sout  = "\n".join(lline).strip()
  return sout

#----------------

for (Year,Mon) in lYM:
    srcPath = listDir + '/list.%s.%s.%04d%02d.csv'%(prod,ver,Year,Mon)
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lines = [line.strip() for line in lines]
    random.seed(0)
    random.shuffle(lines)
    sout = util.list2csv(lines)

    outPath = listDir + '/list.shuffle.%s.%s.%04d%02d.csv'%(prod,ver,Year,Mon)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath

