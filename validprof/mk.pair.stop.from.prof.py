import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket
import sys
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'

else:
    print 'check myhost'
    sys.exit()
#*******************************
iYM  = [2014,6]
eYM  = [2015,5]
lYM  = util.ret_lYM(iYM,eYM)
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
#lrettype = ['epc','gprof']
#lrettype = ['epc']
lrettype = ['gprof']
lvar = ['profpmw','profrad']
#lvar = ['profrad']
#lvar = ['profpmw','profrad','top-profpmw']
thwat = 0.01  # g/m3 for defining stop
nz    = 25
#********************************
def mk_stop(a2prof, thwat):
    a2bin  = np.ones(a2prof.shape).astype(int32) * np.arange(25)
    a2bin  = ma.masked_where(a2prof <=thwat, a2bin).argmax(axis=1)
    return (ma.masked_less_equal(a2bin,0) * 500 + 250).filled(-9999)
    
#********************************
for Year,Mon in lYM:
    for rettype in lrettype:
        for var in lvar:
            if (rettype=='gprof')and(var=='top-profpmw'):
                continue

      
            eDay   = calendar.monthrange(Year,Mon)[1]
            iDTime = datetime(Year,Mon,1)
            eDTime = datetime(Year,Mon,eDay)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
        

            for DTime in lDTime:
                print rettype, var, DTime
                Year,Mon,Day = DTime.timetuple()[:3]
                if [Year,Mon,Day] in lskipdates:
                    continue
                if rettype == 'epc':
                    srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
                elif rettype=='gprof':
                    srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/gprof/%04d/%02d/%02d'%(Year,Mon,Day)
                else:
                    print 'check rettype',rettype
                    sys.exit()
        
                ssearch = srcDir + '/%s.??????.npy'%(var)
                lprofPath = np.sort(glob.glob(ssearch))
    
                if len(lprofPath)==0:
                    print 'No files'
                    print ssearch
                    sys.exit() 
                for profPath in lprofPath:
                    oid = int(profPath.split('.')[-2])    
        
                    a2prof = np.load(srcDir + '/%s.%06d.npy'%(var,oid))[:,:nz] # Stored in bottom to top order.

                    a2stop = mk_stop(a2prof, thwat)

                    outPath = srcDir + '/stop-%s.%06d.npy'%(var, oid)
                    np.save(outPath, a2stop)
                    print outPath



