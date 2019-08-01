from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, os, sys
import numpy as np
import calendar

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
lat0   = -60
lon0   = -180
dlatlon= 1.0
ny,nx  = 120,360
#lvar   = ['profS','profrmse']
#lvar   = ['profS','precrad','precpmw']
#lvar   = ['precrad','precpmw','stoprad','stoppmw','peakhrad','peakhpmw']
lvar   = ['profS','profrmse']+['precrad','precpmw','stoprad','stoppmw','peakhrad','peakhpmw']
#lvar   = ['stoppmw']
#lvar   = ['peakhrad','peakhpmw']

def taylor_index(a2ref, a2dat, miss):
    a2mask1 = ma.masked_less_equal(a2ref, miss).mask
    a2mask2 = ma.masked_less_equal(a2dat, miss).mask
    a2mask  = a2mask1 + a2mask2
    a2ref = ma.masked_where(a2mask, a2ref)
    a2dat = ma.masked_where(a2mask, a2dat)

    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    corrmax= 1.0

    S = 4*(1.0+a1corr)**4 /((a1stdref/a1stddat + a1stddat/a1stdref)**2) / (1.0+corrmax)**4
    S = ma.masked_invalid(S).filled(miss)
    return S

def ret_stormtop(a2prof, miss_out=-9999.):
    thmin = 0.1 # g/m3
    ny,nz = a2prof.shape
    a1h = arange(nz).reshape(1,-1)*0.5 + 0.25
    a2h = np.repeat(a1h, ny, axis=0)
    a1stop = ma.masked_where(a2prof< thmin, a2h).max(axis=1).filled(miss_out)
    return a1stop 

def ret_peakheight(a2prof, miss_out=-9999.):
    thmin = 0.1 # g/m3
    a2prof = ma.masked_less(a2prof, thmin)
    a1mask = a2prof.mask.all(axis=1)
    a1h = a2prof.argmax(axis=1)*0.5 + 0.25
    a1h = ma.masked_where(a1mask, a1h).filled(miss_out)
    return a1h

#**** Month loop ********************
for Year,Mon in lYM:
    eDay = calendar.monthrange(Year,Mon)[1]
    #eDay = 2
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    
    for var in lvar:
        a2var = np.zeros([ny,nx],float32)
        a2num = np.zeros([ny,nx],int32)

        for DTime in lDTime:
            print var, DTime
            Year, Mon, Day = DTime.timetuple()[:3]
            srcDir = '/tank/utsumi/validprof/pair.gprof/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch= srcDir + '/Latitude.*.npy'
            llatPath= sort(glob.glob(ssearch))
        
            for latPath in llatPath:
                #print latPath
                oid = int(latPath.split('.')[-2])
                a1lat = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
                a1lon = np.load(srcDir + '/Longitude.%06d.npy'%(oid))
                a1qflag= np.load(srcDir + '/qualityFlag.%06d.npy'%(oid))        
                a1surftype= np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))        

                if var in ['precbias','precbrat']:
                    a1precrad =np.load(srcDir + '/precrad.%06d.npy'%(oid))
                    a1precpmw=np.load(srcDir + '/precpmw.%06d.npy'%(oid))
        
                elif var in ['profrmse', 'profS']:
                    a2profrad=np.load(srcDir + '/profrad.%06d.npy'%(oid))
                    a2profpmw=np.load(srcDir + '/profpmw.%06d.npy'%(oid))
                    a2profpmw=ma.masked_greater(a2profpmw,30).filled(-9999.) # screen extremely large values from gprof (cmb product ranges 0-18g/m3)

                elif var in ['precrad']:
                    a1precrad =np.load(srcDir + '/precrad.%06d.npy'%(oid))
 
                elif var in ['precpmw']:
                    a1precpmw=np.load(srcDir + '/precpmw.%06d.npy'%(oid))

                elif var in ['stoprad','peakhrad']:
                    a2prof = np.load(srcDir + '/profrad.%06d.npy'%(oid))
                elif var in ['stoppmw','peakhpmw']:
                    a2prof = np.load(srcDir + '/profpmw.%06d.npy'%(oid))

                #print a1precpmw.min(), a1precrad.min(), a1precpmw.max(), a1precrad.max()
                #-- Error ----
                if var=='profrmse':
                    a2profpmw = ma.masked_less(a2profpmw,0)
                    a2profrad = ma.masked_less(a2profrad,0)
                    a1var = np.sqrt(np.square(a2profpmw - a2profrad).mean(axis=1)).filled(-9999.)

                    print a1var.max(), a2profpmw.max(), a2profrad.max()

                elif var=='profS':
                    a1var = taylor_index(a2profrad, a2profpmw, miss=-9999.)

                elif var=='precrad':
                    a1var = ma.masked_less(a1precrad,0).filled(-9999.)
                elif var=='precpmw':
                    a1var = ma.masked_less(a1precpmw,0).filled(-9999.)

                elif var in ['stoprad','stoppmw']:
                    a1var = ret_stormtop(a2prof, miss_out=-9999.)

                elif var in ['peakhrad','peakhpmw']:
                    a1var = ret_peakheight(a2prof, miss_out=-9999.)


                else:
                    print 'check var',var
                    sys.exit()
        
                #-- Projection over grid map --
                a1y = ((a1lat - lat0)/dlatlon).astype(int32)
                a1x = ((a1lon - lon0)/dlatlon).astype(int32)
       
                a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)
                a1flagS1 = ma.masked_equal(a1surftype,1).mask   # Ocean
                a1flagS2 = ma.masked_inside(a1surftype,3,7).mask # Vegetation
                a1flagS3 = ma.masked_inside(a1surftype,12,13).mask # Standing water and rivers & Water/Coast boundary 
                a1flagS  = a1flagS1 + a1flagS2 + a1flagS3

                a1flag1 = ma.masked_inside(a1y,0,119).mask
                a1flag2 = ma.masked_inside(a1x,0,359).mask
                a1flag3 = ma.masked_not_equal(a1var, -9999.).mask 
                a1flag = a1flagQ * a1flagS *a1flag1 * a1flag2 * a1flag3

                a1y = a1y[a1flag]
                a1x = a1x[a1flag]
                a1var= a1var[a1flag]
                for i in range(len(a1var)):
                    y = a1y[i]
                    x = a1x[i]
                    if a1var[i] != -9999.:
                        a2var[y, x] = a2var[y,x] + a1var[i]
                        a2num[y, x] += 1
        imax = np.argmax(a2var)
        #print a2var.max(), imax, a2num.flatten()[imax]
        a2var = a2var / ma.masked_equal(a2num,0)
        a2var = a2var.filled(-9999.)

        imax = np.argmax(a2var)
        #print a2var.max(), imax, a2num.flatten()[imax]

    
        mapDir  = '/tank/utsumi/validprof/maperror'
        util.mk_dir(mapDir)

        mapPath = mapDir + '/ave.%s.%04d%02d.npy'%(var, Year,Mon)
        np.save(mapPath, a2var)

        numPath = mapDir + '/num.%s.%04d%02d.npy'%(var, Year,Mon)
        np.save(numPath, a2num)
        print mapPath


        #-- lat and lon ---
        if ([Year,Mon]==lYM[0])&(var==lvar[0]):
            a1latout = np.arange(lat0, lat0+dlatlon*ny+0.01, dlatlon)+0.5*dlatlon
            a1lonout = np.arange(lon0, lon0+dlatlon*nx+0.01, dlatlon)+0.5*dlatlon
            np.save(mapDir +'/lat.npy', a1latout)
            np.save(mapDir +'/lon.npy', a1lonout)
