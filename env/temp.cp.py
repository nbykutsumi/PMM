import os, sys
import myfunc.util as util
import shutil
from datetime import  datetime, timedelta
import glob

baseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    
    lsrcPath = glob.glob(srcDir + '/*.npy')

    for srcPath in lsrcPath:
   
        filename = os.path.basename(srcPath)
        var    = filename.split('.')[0]
        if var=='full':
            var = filename.split('.')[1]
            
        #outDir = baseDir + '/%s/%04d/%02d/%02d'%(var,Year,Mon,Day)
        #util.mk_dir(outDir)
        #outPath = outDir + '/' + filename
        #shutil.copy(srcPath, outPath)
        #print outPath
