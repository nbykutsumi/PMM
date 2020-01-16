import os
from datetime import datetime, timedelta
import socket
import subprocess
import myfunc.util as util
import os, sys

myhost = socket.gethostname()
if myhost =='shui':
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    figDir   = '/home/utsumi/temp/geo'

elif myhost == 'well':
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()


DTime = datetime(2017,7,1,2,20)


## -------------------------------------------------------
## Quick reference list (Released gridded data and Himawari8 band)
## [EXT] 01:Band03
## [VIS] 01:Band01 02:Band02 03:Band04
## [SIR] 01:Band05 02:Band06
## [TIR] 01:Band13 02:Band14 03:Band15 04:Band16 05:Band07
##       06:Band08 07:Band09 08:Band10 09:Band11 10:Band12
## -------------------------------------------------------


#      wget ${FTP}/${YYYY}${MM}/${CHN}/${YYYY}${MM}${DD}${HH}${MN}.${CHN,,}.${NUM}.fld.geoss.bz2
yyyy,mm,dd,hh,mn = DTime.timetuple()[:5]
baseftp = 'ftp://hmwr829gr.cr.chiba-u.ac.jp/gridded/FD/V20190123/archived/archived2'

ch = 'tir'
chnum = 1
timestamp = '%04d%02d%02d%02d%02d'%(yyyy,mm,dd,hh,mn)
bzName  = '%s.%s.%02d.fld.geoss.bz2'%(timestamp,ch,chnum)

print '-----Download-------'
iPath  = baseftp + '/%04d%02d/%s/'%(yyyy,mm,str.upper(ch)) + bzName

tmpDir = tankbaseDir + '/utsumi/PMM/himawari/temp/temp.%s.%s.%02d'%(timestamp,ch,chnum)
util.mk_dir(tmpDir)
bzPath = tmpDir + '/' + bzName
if os.path.exists(bzPath):
    os.remove(bzPath)

scmd = 'wget %s -P %s'%(iPath, tmpDir)
lcmd = scmd.split(' ')
print ''
print scmd
subprocess.call(lcmd)


print '----- Uncompress -------'
bzPath = tmpDir + '/' + bzName
scmd = 'dd if=%s of=%s/little.geoss conv=swab'%(bzPath, tmpDir)
lcmd = scmd.split(' ')
print ''
print scmd
subprocess.call(lcmd)
print ''
print '----- Convert count to tbb -----'
print ''
para = '%s.%02d'%(ch, chnum)

os.chdir(tmpDir)

if ch=='tir':
    scmd = './tir.x %s/little.geoss %s'%(tmpDir, para)
    print scmd
    subprocess.call(scmd.split(' '))
    resolution="0.02"

elif ch=='vis':
    scmd = './vis.x %s/little.geoss %s'%(tmpDir,para)
    print scmd
    subprocess.call(scmd.split(' '))
    resolution="0.01"

elif ch=='ext':
    scmd = 'dd if=%s/little.geoss of=%s/01.geoss bs=576000000 count=1'%(tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = './ext.x %s/01.geoss %s && mv %s/grid05.dat %s/grid05_1.dat'%(tmpDir, para, tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'dd if=%s/little.geoss of=%s/02.geoss bs=576000000 skip=1'%(tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = './ext.x %s/02.geoss %s && mv %s/grid05.dat %s/grid05_2.dat'%(tmpDir, para, tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'cat %s/grid05_1.dat %s/grid05_2.dat > %s/grid05.dat'%(tmpDir, tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'cat %s/grid05_1.dat %s/grid05_2.dat > %s/grid05.dat'%(tmpDir, tmpDir, tmpDir)
    print scmd
    subprocess.call(scmd.split(' '))





