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
iYM  = [2014,12]
eYM  = [2015,2]
lYM  = util.ret_lYM(iYM,eYM)

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
#lrettype = ['epc']
lrettype = ['gprof']
#lvar= ['stoprad','top-stoppmw']
#lvar= ['stoprad','top-stoppmw']
#lvar = ['precrad','precpmw','zeroDegAltituderad']
lvar = ['precpmw']
#lvar = ['convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['precrad','precpmw','stoprad','top-stoppmw','convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['convfreqpmw','stratfreqpmw']
#lvar = ['zeroDegAltituderad']

lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]


#lstype= ['all','sea','land','veg','snow','coast']
lstype= ['all','sea','veg','snow','coast']
#lstype= ['coast']
lptype= ['all','conv','stra']
lph   = ['H','L','A']
#lprrange=[[0.5,999],[1,3],[8,12]]
lprrange=[[0.5,999]]
lprrange= map(tuple, lprrange)
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

            eDay   = calendar.monthrange(Year,Mon)[1]
            iDTime = datetime(Year,Mon,1)
            eDTime = datetime(Year,Mon,eDay)
            lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
            #lDTime = lDTime[:1]  # test

            lkey = [(stype,ptype,ph,prrange)
                    for stype   in lstype
                    for ptype   in lptype
                    for ph      in lph
                    for prrange in lprrange]

            #** Initialize **********
            d2ss  = {}
            d2sum = {}
            d2num = {}
            for key in lkey:
                d2ss [key] = np.zeros([ny,nx],float32)
                d2sum[key] = np.zeros([ny,nx],float32)
                d2num[key] = np.zeros([ny,nx],int32)
         

            for DTime in lDTime:
                print rettype,var,DTime
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


                    a1stype= np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))
                    a1ptype= np.load(srcDir + '/typePreciprad.%06d.npy'%(oid))
                    #a1ph   = np.load(srcDir + '/stoprad.%06d.npy'%(oid))
                    a1ph   = np.load(srcDir + '/stop-profrad.%06d.npy'%(oid))
                    a1freez= np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid))

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

                    #---------- 
                    a1flagy = ma.masked_inside(a1y,0,ny-1).mask
                    a1flag  = a1flagy
                    #----------
            
                    a1y = a1y[a1flag]
                    a1x = a1x[a1flag]
                    avar= avar[a1flag]
                    a1lat=a1lat[a1flag]
                    a1lon=a1lon[a1flag]

                    a1prec = a1prec[a1flag]
                    a1stype= a1stype[a1flag]
                    a1ptype= a1ptype[a1flag]
                    a1ph   = a1ph   [a1flag]
                    a1freez= a1freez[a1flag]
 
                    avar= ma.masked_less(avar,0).filled(0.0)

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
                        if ph =='L':
                            a1flagph = a1ph < a1freez - 500
                        elif ph=='H':
                            a1flagph = a1ph > a1freez + 500
                        elif ph=='A':
                            a1flagph = np.array([True]*len(a1ph))

                        else:
                            print 'check ph',ph
                            sys.exit()

                        a1flagph = a1flagph * ma.masked_not_equal(a1ph, -9999).mask
                        a1flagph = a1flagph * ma.masked_not_equal(a1freez, -9999).mask


                        #-- Precipitation range -------------
                        if a1prec.shape[0]==0:
                            continue

                        thpr0,thpr1 = prrange
                        a1flagp = ma.masked_inside(a1prec, thpr0, thpr1).mask

                        #-- Screen --------------------------
                        a1flag = a1flagstype * a1flagptype * a1flagph * a1flagp

                        if a1flag.sum()==0:
                            continue

                        avarTmp  = avar[a1flag]
                        a1yTmp   = a1y[a1flag]
                        a1xTmp   = a1x[a1flag]


                        a2sum = d2sum[key]
                        a2ss  = d2ss [key]
                        a2num = d2num[key]

                        for i in range(avarTmp.shape[0]):
                            y = a1yTmp[i]
                            x = a1xTmp[i]

                            a2sum[y,x] = a2sum[y,x] + avarTmp[i]
                            a2ss [y,x] = a2ss[y,x]  + avarTmp[i]**2
                            a2num[y,x] = a2num[y,x] + 1

                        d2sum[key] = a2sum
                        d2ss [key] = a2ss
                        d2num[key] = a2num
       
            #-- Save ---
            for key in lkey:
                stype,ptype,ph,prrange = key
                thpr0,thpr1 = prrange
                a2sum = d2sum[key]
                a2ss  = d2ss [key]
                a2num = d2num[key]

                if rettype=='epc':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
                elif rettype=='gprof':
                    outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
                else:
                    print 'check var (in outDir)',var
                    sys.exit()
        
                stamp  = 's-%s.p-%s.ph-%s.pr-%.1f-%.1f.%04d%02d'%(stype,ptype,ph,thpr0,thpr1,Year,Mon)
                sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)

                np.save(sumPath, a2sum)
                np.save(sum2Path, a2ss)
                np.save(numPath, a2num)
                print sumPath
