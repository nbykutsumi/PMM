import subprocess
import os, sys
import myfunc.util as util

iYM = [2017,8]
eYM = [2017,9]
lYM = util.ret_lYM(iYM,eYM)
#lvar = ['2d','2t','cape','cin','mvimd','sp','tcwv','vo','w']
lvar = ['q','r','t','z']

srcbaseDir = '/tank/utsumi/era5'
outbaseDir = '/media/disk2/share/data/era5'

for (Year,Mon) in lYM:
    for var in lvar:
        srcPath= srcbaseDir + '/%s/%04d%02d'%(var,Year,Mon)
        outDir = outbaseDir + '/%s'%(var)

        scmd = 'rsync -avr utsumi@shui.iis.u-tokyo.ac.jp:%s %s'%(srcPath, outDir)
        print scmd

        subprocess.call(scmd.split())


