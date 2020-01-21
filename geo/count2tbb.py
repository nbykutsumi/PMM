import os
from datetime import datetime, timedelta
import socket
import subprocess
import myfunc.util as util
import os, sys
import shutil

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


argv = sys.argv
if len(argv) != 7+1:
    print 'cmd Year Mon Day Hour Min ch chnum'
    print argv
    sys.exit()

print 'argv[1:]=',argv[1:]
Year,Mon,Day,Hour,Min,ch,chnum = argv[1:]
Year,Mon,Day,Hour,Min, chnum = map(int, [Year,Mon,Day,Hour,Min, chnum])
DTime = datetime(Year,Mon,Day,Hour,Min)

baseftp = 'ftp://hmwr829gr.cr.chiba-u.ac.jp/gridded/FD/V20190123/archived/archived2'


#ch = 'tir'
##ch = 'vis'
##ch = 'ext'
#chnum = 1


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
print scmd
print ''
subprocess.call(lcmd)

if not os.path.exists(bzPath):
    sys.exit()
print '----- Uncompress -------'
scmd = 'bunzip2 %s'%(bzPath)
print scmd
print ''
subprocess.call(scmd.split(' '))



print '----- Swap endian (to little) ------'
bigPath = tmpDir + '/' + bzName[:-4]
scmd = 'dd if=%s of=%s/little.geoss conv=swab'%(bigPath, tmpDir)
lcmd = scmd.split(' ')
print scmd
print ''
subprocess.call(lcmd)

print ''
print '----- Copy parameter files to temporay directry -----'
print ''
para = '%s.%02d'%(ch, chnum)

shutil.copy(para, tmpDir + '/')
print 'copy',para
print 'to: ',tmpDir +'/'


print ''
print '----- Convert count to tbb -----'
print ''

os.chdir(tmpDir)
print os.getcwd()
if ch in ['tir','sir']:
    #scmd = './tir.x %s/little.geoss %s'%(tmpDir, para)
    scmd = '/home/utsumi/bin/PMM/geo/tir_my.x little.geoss %s'%(para)
    print scmd
    subprocess.call(scmd.split(' '))
    resolution="0.02"

elif ch=='vis':
    scmd = '/home/utsumi/bin/PMM/geo/vis_my.x little.geoss %s'%(para)
    print scmd
    subprocess.call(scmd.split(' '))
    resolution="0.01"

elif ch=='ext':
    scmd = 'dd if=little.geoss of=01.geoss bs=576000000 count=1'
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = '/home/utsumi/bin/PMM/geo/ext_my.x 01.geoss %s && mv grid05.dat grid05_1.dat'%(para)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'dd if=little.geoss of=02.geoss bs=576000000 skip=1'
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = '/home/utsumi/bin/PMM/geo/ext_my.x 02.geoss %s && mv grid05.dat grid05_2.dat'%(para)
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'cat grid05_1.dat grid05_2.dat > grid05.dat'
    print scmd
    subprocess.call(scmd.split(' '))

    scmd = 'cat grid05_1.dat grid05_2.dat > grid05.dat'
    print scmd
    subprocess.call(scmd.split(' '))





