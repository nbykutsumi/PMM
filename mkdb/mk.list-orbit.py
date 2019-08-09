import os, sys
import calendar
import glob
from collections import deque
import myfunc.util as util

prod       = '1C'
ver        = 'V05'
gmirootDir = '/work/hk01/PMM/NASA/GPM.GMI/%s/%s'%(prod,ver)

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)

ldataoutage=[
 [2014,5,20]  # 1274-1277
,[2014,10,22] # 3682-
,[2014,10,23] # all
,[2014,10,24] # -3718
,[2014,12,9]  # 004425, 004426(missing some lat&lon)
,[2014,12,10]  # 004447-004452
,[2017,9,26]
,[2017,9,27]
,[2017,9,28]
,[2017,9,29]
]


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
    iDay = 1
    eDay = calendar.monthrange(Year,Mon)[1]

    l = deque([])
    for Day in range(iDay,eDay+1):
        if [Year,Mon,Day] in ldataoutage:
            print 'skip',Year,Mon,Day
            continue
        srcDir  = gmirootDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = srcDir + '/*.HDF5'
        lpath = glob.glob(ssearch)
        print lpath

        for spath in lpath:
            sdate,obtnum = spath.split('/')[-1].split('.')[4:6]
            ymd, time0, time1 = sdate.split('-')
            time0 = time0[1:]
            time1 = time1[1:]
            l.append([obtnum, '%04d'%(Year),'%02d'%(Mon),'%02d'%(Day),time0,time1])

    l = sorted(l, key=lambda x: x[0]) 
    sout = list2csv(l)
    listPath = listDir + '/list.%s.%s.%04d%02d.csv'%(prod, ver, Year,Mon)
    f=open(listPath,'w'); f.write(sout); f.close()
    print listPath
