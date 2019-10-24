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
eYM  = [2015,2]
lYM  = util.ret_lYM(iYM,eYM)
lYM  = [ym for ym in lYM if ym[1] not in [9,10,11]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
lrettype = ['epc']
#lvar= ['stoprad','top-stoppmw']
#lvar= ['stoprad','top-stoppmw']
#lvar = ['precrad','precpmw']
#lvar = ['convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['precrad','precpmw','stoprad','top-stoppmw','convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['convfreqpmw','stratfreqpmw']
lvar = ['zeroDegAltituderad']

lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360
#lthpr = [0.1, 10]
lthpr = [0.5, 10]
lskipdates = [[2014,12,9],[2014,12,10]]

#-----------------
def ret_var_filename(var):
    if var in ['convfreqrad','stratfreqrad']:
        var_filename= 'typePreciprad'
    elif var in ['convfreqpmw','stratfreqpmw']:
        var_filename= 'top-typePrecippmw'
    else:
        var_filename= var
    return var_filename

#-----------------

for Year,Mon in lYM:
    for rettype in lrettype:
        for var in lvar:
            var_filename = ret_var_filename(var)
            for thpr in lthpr:
                eDay   = calendar.monthrange(Year,Mon)[1]
                iDTime = datetime(Year,Mon,1)
                eDTime = datetime(Year,Mon,eDay)
                lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
                #lDTime = lDTime[:1]  # test
         
                a2ss  = np.zeros([ny,nx],float32)
                a2sum = np.zeros([ny,nx],float32)
                a2num = np.zeros([ny,nx],int32)
                for DTime in lDTime:
                    print var,DTime
                    Year,Mon,Day = DTime.timetuple()[:3]
                    if [Year,Mon,Day] in lskipdates: continue
        
                    if rettype =='epc':
                        srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
                    elif rettype=='gprof':
                        srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/gprof/%04d/%02d/%02d'%(Year,Mon,Day)
                    else:
                        print 'check var (in srcDir)',var
                        sys.exit()
        
                    ssearch = srcDir + '/%s.??????.npy'%(var_filename)
                    lprofPath = np.sort(glob.glob(ssearch))
                    for profPath in lprofPath:
                        oid = int(profPath.split('.')[-2])    
        
                        avar  = np.load(srcDir + '/%s.%06d.npy'%(var_filename,oid))
                        a1lat  = np.load(srcDir + '/Latitude.%06d.npy'%(oid)).astype(float64)
                        a1lon  = np.load(srcDir + '/Longitude.%06d.npy'%(oid)).astype(float64)
                        if var[-3:] =='rad':
                            a1prec = np.load(srcDir + '/precrad.%06d.npy'%(oid))
                        elif var[-3:]=='pmw':
                            a1prec = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
                        else:
                            print 'check var',var
                            sys.exit()
       
                        #-- typePrecip ----
                        if var in ['convfreqrad','convfreqpmw','stratfreqrad','stratfreqpmw']:
                            avar = ma.masked_less(avar,0).astype('int16')
                            strat = (avar %10).astype('int16')
                            conv  = (avar%100-strat).astype('int16')/10
                            other = (avar/100).astype('int16')

                            if var in ['convfreqrad','convfreqpmw']:
                                avar = conv
                            elif var in ['stratfreqrad','stratfreqpmw']:
                                avar = strat
                            
                            avar = avar/9.0  # fraction in 9-grids
                        #-- Screen invalid (nan) ----- 
                        avar  = ma.masked_invalid(avar).filled(miss_out) 
        
                        #-- Projection over grid map --
                        a1y = np.floor((a1lat - lat0)/dlatlon).astype(int32)
                        a1x = np.floor((a1lon - lon0)/dlatlon).astype(int32)
                        a1flagP = ma.masked_greater_equal(a1prec,thpr).mask
                        a1flag = a1flagP
                        #---------- 
                        a1flagy = ma.masked_inside(a1y,0,ny-1).mask
                        a1flag  = a1flag * a1flagy
                        #----------
            
                        a1y = a1y[a1flag]
                        a1x = a1x[a1flag]
                        avar= avar[a1flag]
                        a1lon=a1lon[a1flag]
               
                        avar= ma.masked_less(avar,0).filled(0.0)
                        for i in range(avar.shape[0]):
                            y = a1y[i]
                            x = a1x[i]
           
                             
                            a2sum[y,x] = a2sum[y,x] + avar[i]
                            a2ss [y,x] = a2ss[y,x]  + avar[i]**2
                            a2num[y,x] = a2num[y,x] + 1
        
                #-- Save ---
                if rettype=='epc':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
                elif rettype=='gprof':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
                else:
                    print 'check var (in outDir)',var
                    sys.exit()
        
                util.mk_dir(outDir)
                stamp  = 'pr%.1f.%04d%02d'%(thpr,Year,Mon)
                sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
                np.save(sumPath, a2sum)
                np.save(sum2Path, a2ss)
                np.save(numPath, a2num)
                print sumPath
