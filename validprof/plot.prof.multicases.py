from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, os, sys
import numpy as np
import calendar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

iYM = [2017,1]
eYM = [2017,1]
lYM = util.ret_lYM(iYM,eYM)
lat0   = -60
lon0   = -180
dlatlon= 1.0
ny,nx  = 120,360

#**** Functions *********************
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

#**** Month loop ********************
for Year,Mon in lYM:
    eDay = calendar.monthrange(Year,Mon)[1]
    #eDay = 2
    iDTime = datetime(Year,Mon,1)
    eDTime = datetime(Year,Mon,eDay)
    dDTime = timedelta(days=1)
    lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
    #lDTime = lDTime[:1]
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
    
            a1precrad =np.load(srcDir + '/precrad.%06d.npy'%(oid))
            a1precpmw=np.load(srcDir + '/precpmw.%06d.npy'%(oid))
    
            a2profrad=np.load(srcDir + '/profrad.%06d.npy'%(oid))
            a2profpmw=np.load(srcDir + '/profpmw.%06d.npy'%(oid))

            ny,nx = 3, 2  # figure dimension
            w     = 0.8/nx
            h     = 0.8/ny
            lbbox = [[[-50,-30],[-45,-25]],
                     [[-20,-30],[-15,-25]],
                     [[-2.5,-30],[2.5, -25]]]

            fig = plt.figure(figsize=(15,10))
            axcount = 0
            for ibbox, bbox in enumerate(lbbox):
                [[lllat,lllon],[urlat,urlon]] = bbox
                a1flaglat = ma.masked_inside(a1lat, lllat,urlat).mask
                a1flaglon = ma.masked_inside(a1lon, lllat,urlon).mask
                a1flagrad = ma.masked_greater(a1precrad,0).mask
                a1flagprof= ma.masked_not_equal(a2profrad,-9999.9).mask.any(1) * ma.masked_not_equal(a2profpmw,-9999.9).mask.any(1)
                a1flag    = a1flaglat * a1flaglon * a1flagrad * a1flagprof

                print 'flag.sum()=',a1flag.sum()
                print 'A'
                if a1flag is np.bool_(False):
                    continue 

                if a1flag.sum()==0:
                    continue

                print 'B'
                alattmp = a1lat[a1flag]
                alontmp = a1lon[a1flag]

                a1precradtmp = a1precrad[a1flag]
                a2profradtmp = a2profrad[a1flag,:]
                a2profpmwtmp = a2profpmw[a1flag,:]
                #a2mask1 = ma.masked_less_equal(a2profradtmp, -9999.9).mask
                #a2mask2 = ma.masked_less_equal(a2profpmwtmp, -9999.9).mask
                #a2mask  = a2mask1 + a2mask2
                #a2profradtmp = ma.masked_where(a2mask, a2profradtmp)
                #a2profpmwtmp = ma.masked_where(a2mask, a2profpmwtmp)

                a2profradtmp = ma.masked_less(a2profradtmp,0)
                a2profpmwtmp = ma.masked_less(a2profpmwtmp,0)

                #if a2mask is np.bool_(False):
                #    continue

                #print 'C'
                #a1s = taylor_index(a2profradtmp, a2profpmwtmp, miss=-9999.9)
                a1profradtmp = a2profradtmp[0,:]
                a1profpmwtmp = a2profpmwtmp[0,:]
                #-- figure ---
                x = ibbox/ny
                y = ny*x - ibbox
                ax = fig.add_axes([0.05+x,0.05+y, w,h])
                
                a1y   = np.arange(36)
                ax.plot(a1profradtmp, a1y, '-',color='k') 
                ax.plot(a1profpmwtmp, a1y, '--',color='r')
                plt.xlim([0,1.0])
                stitle = 'lat=%s:%s lon=%s:%s'%(lllat,urlat, lllon, urlon)
                axcount += 1 
            #-- save figure --
            figDir = '/home/utsumi/temp/validprof/caseprof'
            util.mk_dir(figDir)
            figPath= figDir + '/prof.%02d%02d.%06d.png'%(Mon,Day,oid)

            if axcount >0:
                plt.savefig(figPath)
                print figPath
                print ''
            else:
                plt.clf()
                print 'No plots',oid

            #if oid == 16194:
            #if oid == 16164:
            #if oid == 16156:  # with fig
            if oid == 16163:
                print 'a1profradtmp',a1profradtmp
                print 'a1profpmwtmp',a1profpmwtmp
                sys.exit()
