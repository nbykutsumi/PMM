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
eYM  = [2014,8]
lYM  = util.ret_lYM(iYM,eYM)

nrec   = 20000
miss_out= -9999.
#lrettype= ['epc-top']
#lrettype= ['gprof']
#lrettype= ['epc','rad','gprof']
lrettype= ['rad-epc','rad-gprof']
#lrettype= ['rad','gprof','epc']
#lrettype= ['gprof','epc']
lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360
thpr = 0.5
for Year,Mon in lYM:
    for rettype in lrettype:
        if rettype =='rad':
            #dirtype = 'gprof'
            dirtype = 'epc'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype =='rad-epc':
            dirtype = 'epc'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype =='rad-gprof':
            dirtype = 'gprof'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype=='gprof':
            dirtype = 'gprof'
            filetype= 'pmw'
            profname= 'prof'
            nz = 36
        elif rettype=='epc':
            dirtype = 'epc'
            filetype= 'pmw'
            profname= 'prof'
            nz = 25
        elif rettype=='epc-top':
            dirtype = 'epc'
            filetype= 'pmw'
            profname= 'top-prof'
            nz = 25

    
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
            srcDir  = tankbaseDir + '/utsumi/validprof/pair.%s/%04d/%02d/%02d'%(dirtype,Year,Mon,Day)
            ssearch = srcDir + '/prof%s.??????.npy'%(filetype)
            lprofPath = np.sort(glob.glob(ssearch))
            for profPath in lprofPath:
                oid = int(profPath.split('.')[-2])    

                a2var  = np.load(srcDir + '/%s%s.%06d.npy'%(profname,filetype,oid))
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
                if dirtype =='gprof':
                    a1flagQ = ma.masked_equal(np.load(srcDir + '/qualityFlag.%06d.npy'%(oid)) ,0).mask   # Good quality (flag=0)
                    a1flag = a1flagP * a1flagQ

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
   
        outDir = tankbaseDir + '/utsumi/validprof/mapprof/%s'%(rettype)
        util.mk_dir(outDir)
    
        sumPath= outDir  + '/prof.sum.%04d%02d.sp.one.npy' %(Year,Mon)
        numPath= outDir  + '/prof.num.%04d%02d.sp.one.npy' %(Year,Mon)
        numprofPath= outDir  + '/prof.numprof.%04d%02d.sp.one.npy' %(Year,Mon)
        sum2Path= outDir + '/prof.sum2.%04d%02d.sp.one.npy'%(Year,Mon)
        avePath= outDir  + '/prof.ave.%04d%02d.sp.one.npy' %(Year,Mon)
        np.save(sumPath, a3sum)
        np.save(sum2Path, a3ss)
        np.save(numPath, a2num)
        np.save(numprofPath, a3numprof)
        np.save(avePath, a3ave)
        print avePath
