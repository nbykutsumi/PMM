import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket

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
#lrettype = ['gprof']
#lrettype = ['epc']
lrettype = ['epc','gprof']
#lvar = ['profpmw','profrad','top-profpmw']
lvar = ['profpmw','profrad']
#lvar = ['profrad']
for Year,Mon in lYM:
    for rettype in lrettype:
        for var in lvar:
            if (rettype=='gprof')and(var=='top-profpmw'):
                continue

            if var in ['profpmw','profrad','top-profpmw']:
                nz = 25
       
            eDay   = calendar.monthrange(Year,Mon)[1]
            iDTime = datetime(Year,Mon,1)
            eDTime = datetime(Year,Mon,eDay)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
        

            for DTime in lDTime:
                print rettype, var, DTime
                Year,Mon,Day = DTime.timetuple()[:3]
                #if [Year,Mon,Day] in lskipdates:
                #    continue
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
                    #sys.exit() 
                    continue
                for profPath in lprofPath:
                    oid = int(profPath.split('.')[-2])    
        
                    a2var  = np.load(srcDir + '/%s.%06d.npy'%(var,oid))[:,:nz] # Stored in bottom to top order. (500m resolution)
                    #a1prec = np.load(srcDir + '/precrad.%06d.npy'%(oid))
       
                    
                    #---- Mask if q(2-2.5km) or q(4-4.5km) is missing ---
                    a1mask1 = ma.masked_less(a2var[:,4],0).mask
                    a1mask2 = ma.masked_less(a2var[:,8],0).mask
                    a1mask  = a1mask1 + a1mask2
    
                    #---- Slope q(2-2.5km) - q(4-4.5km) ----
                    a1slope = ma.masked_where(a1mask, (a2var[:,4] - a2var[:,8])/2.0).filled(miss_out)   # g/m3 / km
    
                    #---- Save -----------------------------
                    outvar= {'profrad':'sloperad', 'profpmw':'slopepmw', 'top-profpmw':'top-slopepmw'}[var]
                    outPath = srcDir + '/%s.%06d.npy'%(outvar, oid)
                    print outPath
    
                    np.save(outPath, a1slope.astype('float32'))            
    
    
    
