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
iYM  = [2014,12]
eYM  = [2015,5]
lYM  = util.ret_lYM(iYM,eYM)
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
#lrettype = ['gprof']
lrettype = ['epc']
#lrettype = ['epc','gprof']
lvar = ['profpmw','profrad']
#lvar = ['profrad']
#lvar = ['profpmw','profrad','top-profpmw']
#lvar = ['top-profpmw']
lstype= ['all','sea','land','veg','snow','coast']
lptype= ['all','conv','stra']
#lprrange=[[0.5,999],[1,3],[8,12]]
lprrange=[[0.5,999]]
lprrange= map(tuple, lprrange)
lph     = ['L','H','A']
#lph     = ['L','H']
#lph     = ['A']
lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360

for Year,Mon in lYM:
    for rettype in lrettype:
        for var in lvar:
            if (rettype=='gprof')and(var=='top-profpmw'):
                continue

            if var in ['profpmw','profrad','top-profpmw']:
                nz = 25
            else:
                nz = 1
       
            eDay   = calendar.monthrange(Year,Mon)[1]
            iDTime = datetime(Year,Mon,1)
            eDTime = datetime(Year,Mon,eDay)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
        
            lkey = [(stype,ptype,ph,prrange)
                            for stype in lstype
                            for ptype in lptype
                            for ph    in lph
                            for prrange in lprrange]

            #** Initialize **********
            d3ss  = {}
            d3sum = {}
            d2num = {}
            d3numprof= {}
            for key in lkey:
                d3ss [key] = np.zeros([ny,nx,nz],float32)
                d3sum[key] = np.zeros([ny,nx,nz],float32)
                d2num[key] = np.zeros([ny,nx],int32)
                d3numprof[key]=np.zeros([ny,nx,nz],int32)


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
        
                    a2var  = np.load(srcDir + '/%s.%06d.npy'%(var,oid))[:,:nz] # Stored in bottom to top order.
                    a1lat  = np.load(srcDir + '/Latitude.%06d.npy'%(oid)).astype(float64)
                    a1lon  = np.load(srcDir + '/Longitude.%06d.npy'%(oid)).astype(float64)
                    a1prec = np.load(srcDir + '/precrad.%06d.npy'%(oid))
        
        
                    a1stype= np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))
                    a1ptype= np.load(srcDir + '/typePreciprad.%06d.npy'%(oid)) 
                    #a1ph   = np.load(srcDir + '/stoprad.%06d.npy'%(oid)) 
                    a1ph   = np.load(srcDir + '/heightStormToprad.%06d.npy'%(oid)) 
                    a1freez= np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid)) 
                    a1dprx = np.load(srcDir + '/dprx.%06d.npy'%(oid)) 

                    #print a2var.shape[0], a1lat.shape, a1lon.shape, a1prec.shape, a1stype.shape, a1ptype.shape, a1ph.shape, a1freez.shape, a1dprx.shape

                    #-- Screen invalid (nan) ----- 
                    a2var  = ma.masked_invalid(a2var).filled(miss_out) 
                    #-- Projection over grid map --
                    a1y = np.floor((a1lat - lat0)/dlatlon).astype(int32)
                    a1x = np.floor((a1lon - lon0)/dlatlon).astype(int32)
       
                    a1x = ma.masked_equal(a1x,-1).filled(nx-1)
                    a1x = ma.masked_equal(a1x,nx).filled(0) 
                
                    #---------- 
                    if rettype =='gprof':
                        a1flagQ = ma.masked_equal(np.load(srcDir + '/qualityFlag.%06d.npy'%(oid)) ,0).mask   # Good quality (flag=0)
       
                        a1flagS = ma.masked_equal(np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid)) ,2)
                        a1flagS = ma.masked_inside(a1flagS,8,11)
                        a1flagS = ma.masked_equal(a1flagS,14)
                        a1flagS = ~(a1flagS.mask)
                        a1flag  = a1flagQ * a1flagS 
        
                    else:
                        a1flag = True

                    #-- DPR angle bin --
                    a1flagdprx = ma.masked_inside(a1dprx, 14, 34).mask # nadir=24
                    a1flag  = a1flag * a1flagdprx

                    #---------- 
                    a1flagy = ma.masked_inside(a1y,0,ny-1).mask
                    a1flag  = a1flag * a1flagy
                    #----------
                    #print oid,a1flag.shape, a1y.shape, a2var.shape, a1lat.shape, a1prec.shape, a1stype.shape, a1ptype.shape, a1ph.shape, a1freez.shape

            
                    a1y = a1y[a1flag]
                    a1x = a1x[a1flag]
                    a2var=a2var[a1flag,:]
                    a1lat  = a1lat[a1flag]
                    a1lon  = a1lon[a1flag]
                    a1prec = a1prec[a1flag]
                    a1stype= a1stype[a1flag]
                    a1ptype= a1ptype[a1flag]
                    a1ph   = a1ph[a1flag]
                    a1freez= a1freez[a1flag]
 
                    a2bit= ma.masked_greater_equal(a2var,0).mask.astype(int32) 
                    a2var= ma.masked_less(a2var,0).filled(0.0)

                    if a2bit.max()==0:
                        continue

                    #-- Precip type ----------
                    a1ptype = ma.masked_less(a1ptype,0).astype('int16')
                    a1stra  = (a1ptype %10).astype('int16')
                    a1conv  = (a1ptype%100-a1stra).astype('int16')/10
                    a1other = (a1ptype/100).astype('int16')
                    a1all   = (a1stra + a1conv + a1other).astype(float32) 
                    a1stra  = a1stra / a1all
                    a1conv  = a1conv / a1all 

                    #-- Surface type flag -------------------
                    a1flagsea  = ma.masked_equal(a1stype,1).mask
                    a1flagveg  = ma.masked_inside(a1stype,3,7).mask
                    a1flagsnow = ma.masked_inside(a1stype,8,11).mask
                    a1flagcoast= ma.masked_equal(a1stype,13).mask
                    a1flagland = a1flagveg + a1flagsnow
                    #-- Precipitation type flag -------------
                    a1flagconv = ma.masked_greater(a1conv,0.6).mask
                    a1flagstra = ma.masked_greater(a1stra,0.6).mask
                    for key in lkey:
                        stype,ptype,ph,prrange = key
                        #-- Surface type ---
                        if stype == 'sea':
                            a1flagstype = a1flagsea
                        elif stype=='veg':
                            a1flagstype = a1flagveg
                        elif stype=='snow':
                            a1flagstype = a1flagsnow
                        elif stype=='coast':
                            a1flagstype = a1flagcoast
                        elif stype=='land':
                            a1flagstype = a1flagland
                        elif stype=='all':
                            a1flagstype = np.array([True]*len(a1stype))
                        else:
                            print 'check stype',stype
                            sys.exit()

                        if type(a1flagstype) is np.bool_:
                            a1flagstype = np.array([a1flagstype]*len(a1stype))


                        #-- Precip type ----
                        if ptype == 'conv':
                            a1flagptype = a1flagconv
                        elif ptype=='stra':
                            a1flagptype = a1flagstra
                        elif ptype=='all':
                            a1flagptype = np.array([True]*len(a1ptype))
                        else:
                            print 'stype',stype
                            sys.exit()

                        if type(a1flagptype) is np.bool_:
                            a1flagptype = np.array([a1flagptype]*len(a1ptype))

                        #-- Precipitation height ------------
                        a1fpd = ma.masked_less(a1ph,0) - ma.masked_less(a1freez,0)  # freezing precipitation depth
                        if ph =='L':
                            a1flagph = ma.masked_less(a1fpd, -0.5).mask  # 500m lower than FL
                        elif ph=='H':
                            a1flagph = ma.masked_greater(a1fpd, 0.5).mask # 500m higher than FL 
                        elif ph=='A':
                            a1flagph = np.array([True]*len(a1fpd))
                        
                        else:
                            print 'check ph',ph
                            sys.exit()

                        #-- Precipitation range -------------
                        thpr0,thpr1 = prrange
                        a1flagp = ma.masked_inside(a1prec, thpr0, thpr1).mask
                        
                        #-- Screen --------------------------
                        a1flag = a1flagstype * a1flagptype * a1flagph * a1flagp

                        if a1flag.sum()==0:
                            continue 

                        a2varTmp = a2var[a1flag,:]
                        a2bitTmp = a2bit[a1flag,:]
                        a1yTmp   = a1y[a1flag]
                        a1xTmp   = a1x[a1flag]


                        a3sum = d3sum[key]
                        a3ss  = d3ss [key]
                        a2num = d2num[key]
                        a3numprof = d3numprof[key]

                        for i in range(a2varTmp.shape[0]):
                            y = a1yTmp[i]
                            x = a1xTmp[i]
            
                            a3sum[y,x,:] = a3sum[y,x,:] + a2varTmp[i]
                            a3ss [y,x,:] = a3ss[y,x,:]  + a2varTmp[i]**2
                            a2num[y,x] = a2num[y,x] + 1
                            a3numprof[y,x,:] = a3numprof[y,x,:] + a2bitTmp[i]

                        d3sum[key] = a3sum
                        d3ss [key] = a3ss
                        d2num[key] = a2num
                        d3numprof[key] = a3numprof
 
            #-- Save ---
            for key in lkey:
                stype,ptype, ph, prrange = key
                thpr0,thpr1 = prrange
                a3sum = d3sum[key]
                a3ss  = d3ss [key]
                a2num = d2num[key]
                a3numprof = d3numprof[key]

                if rettype =='epc':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
                elif rettype =='gprof':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
                else:
                    print 'check rettype (in outDir)',rettype
                    sys.exit()
        
                util.mk_dir(outDir)
                stamp  = 's-%s.p-%s.ph-%s.pr-%.1f-%.1f.%04d%02d'%(stype,ptype,ph,thpr0,thpr1,Year,Mon)
                sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                numprofPath= outDir  + '/%s.numprof.%s.sp.one.npy' %(var,stamp)
                sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
                np.save(sumPath, a3sum)
                np.save(sum2Path, a3ss)
                np.save(numPath, a2num)
                np.save(numprofPath, a3numprof)
                print sumPath
