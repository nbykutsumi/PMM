from ftplib import FTP
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket
from datetime import datetime, timedelta
import numpy as np

hostname  = "ftp.eorc.jaxa.jp"
irootDir = "/pub/TRMM/hidden/GSMaP_L2"

myhost = socket.gethostname()
if myhost =="shui":
    orootDir   = "/work/hk02/PMM/GSMaP_L2"
elif myhost =="well":
    orootDir   = "/home/utsumi/mnt/lab_work/hk02/PMM/GSMaP_L2"


iDTime  = datetime(2018,1,1)
eDTime  = datetime(2018,12,31)

#lver = ['V7','V6']
lver = ['V6']
lname = [
#'amsr2',
#'gmi',
#'metop_A',
#'metop_B',
#'metop_N18',
#'metop_N19',
'ssmis_F16',
'ssmis_F17',
'ssmis_F18',
]

myid =''
mypass = ''
#----------------------------------
def mk_dir(sodir):
  try:
    os.makedirs(sodir)
  except OSError:
    pass

ftp = FTP(hostname)
ftp.login(myid, mypass)


#----------------------------------
for [name,ver] in [[name,ver] for ver in lver for name in lname]:
    #--- path and directory: Remote -----------
    iDir     = irootDir + '/STD_%s/%s'%(ver, name)
    obaseDir = orootDir + '/STD_%s/%s'%(ver, name) 
    #GPM/TRMM.TMI/L2A12/07/2014/
    #--- list --------------
    lPath = np.sort(ftp.nlst(iDir))
    print iDir
    #print lPath
    for sPath in lPath:
        fname = os.path.basename(sPath)
        stime = fname.split('_')[2]
        yy = int(stime[:2])
        if yy>=97:
            yyyy = 1900 + yy
        else:
            yyyy = 2000 + yy
        mm = int(stime[2:4])
        dd = int(stime[4:6])
        hh = int(stime[6:8])
        mn = int(stime[8:10])
        dtime = datetime(yyyy,mm,dd,hh,mn)

        if ((dtime < iDTime)or(eDTime+timedelta(days=1) < dtime)): continue

        #--- test -------
        #if (name=='amsr2'):continue
        #if (name=='gmi')&(dtime<datetime(2018,1,17)):continue
        #if (name=='metop_A')&(dtime<datetime(2018,10,12)):continue
        if (name=='ssmis_F16')&(dtime<datetime(2018,1,20)):continue
        #----------------

        print dtime

        oDir  = obaseDir + '/%04d/%02d/%02d'%(yyyy,mm,dd)
        oPath = oDir + "/" + fname
        mk_dir(oDir)
        ftp.retrbinary("RETR %s"%(sPath), open(oPath,"wb").write)
        print oPath

ftp.close()

