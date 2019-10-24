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
iYM  = [2015,2]
eYM  = [2015,2]
lYM  = util.ret_lYM(iYM,eYM)
lYM  = [ym for ym in lYM if ym[1] not in [9,10,11]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
#lrettype = ['gprof']
lrettype = ['epc']
#lrettype = ['epc','gprof']
lvar = ['profpmw','profrad','top-profpmw']
#lvar = ['top-profpmw']
lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360
lthpr = [0.5, 10]
#lthpr = [0.5]
#lthpr = [10]
lskipdates = [[2014,12,9],[2014,12,10]]
for Year,Mon in lYM:
    for rettype in lrettype:
        for var in lvar:
            if (rettype=='gprof')and(var=='top-profpmw'):
                continue

            for thpr in lthpr:
                if var in ['profpmw','profrad','top-profpmw']:
                    nz = 25
                else:
                    nz = 1
        
                eDay   = calendar.monthrange(Year,Mon)[1]
                iDTime = datetime(Year,Mon,1)
                eDTime = datetime(Year,Mon,eDay)
                lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
                #lDTime = lDTime[:1]  # test
         
                a3ss  = np.zeros([ny,nx,nz],float32)
                a3sum = np.zeros([ny,nx,nz],float32)
                a2num = np.zeros([ny,nx],int32)
                a3numprof=np.zeros([ny,nx,nz],int32)
                for DTime in lDTime:
                    print rettype, DTime
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
        
                        a2var  = np.load(srcDir + '/%s.%06d.npy'%(var,oid))[:,:nz] # Stored in bottom to top order.
                        a1lat  = np.load(srcDir + '/Latitude.%06d.npy'%(oid)).astype(float64)
                        a1lon  = np.load(srcDir + '/Longitude.%06d.npy'%(oid)).astype(float64)
                        a1prec = np.load(srcDir + '/precrad.%06d.npy'%(oid))
        
        
                        #-- Screen invalid (nan) ----- 
                        a2var  = ma.masked_invalid(a2var).filled(miss_out) 
                        #-- Projection over grid map --
                        a1y = np.floor((a1lat - lat0)/dlatlon).astype(int32)
                        a1x = np.floor((a1lon - lon0)/dlatlon).astype(int32)
                        a1flagP = ma.masked_greater(a1prec,thpr).mask
        
                        #---------- 
                        if rettype =='gprof':
                            a1flagQ = ma.masked_equal(np.load(srcDir + '/qualityFlag.%06d.npy'%(oid)) ,0).mask   # Good quality (flag=0)
       
                            a1flagS = ma.masked_equal(np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid)) ,2)
                            a1flagS = ma.masked_inside(a1flagS,8,11)
                            a1flagS = ma.masked_equal(a1flagS,14)
                            a1flagS = ~(a1flagS.mask)
 
                            a1flag = a1flagP * a1flagQ * a1flagS
        
                        else:
                            a1flag = a1flagP
                        #---------- 
                        a1flagy = ma.masked_inside(a1y,0,ny-1).mask
                        a1flag  = a1flag * a1flagy
                        #----------
            
                        a1y = a1y[a1flag]
                        a1x = a1x[a1flag]
                        a2var=a2var[a1flag,:]
                        a1lon=a1lon[a1flag]
                        #print rettype, thpr,oid,a1flagP.sum()
               
                        a2bit= ma.masked_greater_equal(a2var,0).mask.astype(int32) 
                        a2var= ma.masked_less(a2var,0).filled(0.0)

                        for i in range(a2var.shape[0]):
                            y = a1y[i]
                            x = a1x[i]
            
                            a3sum[y,x,:] = a3sum[y,x,:] + a2var[i]
                            a3ss [y,x,:] = a3ss[y,x,:]  + a2var[i]**2
                            a2num[y,x] = a2num[y,x] + 1
                            a3numprof[y,x,:] = a3numprof[y,x,:] + a2bit[i]
       
                #-- Save ---
                a3ave = ma.masked_invalid( a3sum/a2num.reshape(ny,nx,1) ).filled(-9999.)
           
                if rettype =='epc':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
                elif rettype =='gprof':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
                else:
                    print 'check rettype (in outDir)',rettype
                    sys.exit()
        
                util.mk_dir(outDir)
                stamp  = 'pr%.1f.%04d%02d'%(thpr,Year,Mon)
                sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                numprofPath= outDir  + '/%s.numprof.%s.sp.one.npy' %(var,stamp)
                sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
                np.save(sumPath, a3sum)
                np.save(sum2Path, a3ss)
                np.save(numPath, a2num)
                np.save(numprofPath, a3numprof)
                print sumPath
