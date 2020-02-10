import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
#calcflag = False
calcflag = True
iYM  = [2014,7]
eYM  = [2014,7]
lYM  = util.ret_lYM(iYM,eYM)
lskipdates = [[2014,6,4],[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
#lrettype = ['gprof']
lrettype = ['epc']
#lrettype = ['epc','gprof']
lvar = ['profrad','profpmw']
#lvar = ['profrad']
#lvar = ['profpmw']
#lvar = ['profpmw','profrad','top-profpmw']
#lvar = ['top-profpmw']
#lstype= ['all','sea','land','veg','snow','coast']
lstype= ['veg']
#lptype= ['all','conv','stra']
lptype= ['conv']
#lprrange=[[0.5,999],[1,3],[8,12]]
lprrange=[[0.5,999]]
lprrange= map(tuple, lprrange)
#lph     = ['L','H','A']
#lph     = ['L','H']
#lph     = ['L']
lph     = ['A']
lat0 = -60.
lon0 = -180.
dlatlon=1.0
ny,nx = 120,360

#region = 'TRO'
region = 'MIDN'

dBBox = {
         'TRO':  [[-15,-180],[15,180]]
        ,'SUBN': [[15,-180],[30,180]]
        ,'MIDN': [[35,-180],[50,180]]
        ,'AMZ':  [[-5,-70],[5,-53]]
        ,'CUS':  [[35,-105],[45,-95]]
        ,'EUS':  [[30,-90],[40,-80]]
        ,'TIB':  [[30,85],[35,95]]
        ,'NETP':  [[0,-120],[10,-110]]
        ,'SETP':  [[-10,-120],[0,-110]]
        ,'NTA' : [[0,-35],[10,-25]]
        ,'STA' : [[-10,-35],[0,-25]]
        ,'WTP':  [[0,140],[10,150]]
        ,'ETI':  [[-5,85],[5,95]]
        ,'WMP':  [[30,143],[40,158]]
        ,'WMA':  [[30,-74],[40,-59]]
        ,'TAF':  [[0, 15],[15,30]]
        ,'NEA':  [[40,120],[50,130]]
        ,'SEC':  [[22,105],[37,115]]
        ,'NIN':  [[20,75],[25,85]]
        }

    
#********************************
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


            for key in lkey:

                stype,ptype,ph,prrange = key

                #** Initialize **********
                a2prof = None

                #------------------------

                for DTime in lDTime:

                    if calcflag==False: continue

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
                        #a1ph   = np.load(srcDir + '/heightStormToprad.%06d.npy'%(oid))
                        a1ph   = np.load(srcDir + '/stop-profrad.%06d.npy'%(oid))
                        a1freez= np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid)) 
                        a1dprx = np.load(srcDir + '/dprx.%06d.npy'%(oid)) 


                        #-- Screen invalid (nan) ----- 
                        a2var  = ma.masked_invalid(a2var).filled(miss_out) 
                        #-- Projection over grid map --
                        #-- DPR angle bin --
                        a1flagdprx = ma.masked_inside(a1dprx, 14, 34).mask # nadir=24
                        a1flag  = a1flagdprx

                        #---------- 
                
                        a2var=a2var[a1flag,:]
                        a1lat  = a1lat[a1flag]
                        a1lon  = a1lon[a1flag]
                        a1prec = a1prec[a1flag]
                        a1stype= a1stype[a1flag]
                        a1ptype= a1ptype[a1flag]
                        a1ph   = a1ph[a1flag]
                        a1freez= a1freez[a1flag]
 
                        a2bit= ma.masked_greater_equal(a2var,0).mask.astype(int32) 
                        #a2var= ma.masked_less(a2var,0).filled(0.0)

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
                        thpr0,thpr1 = prrange
                        a1flagp = ma.masked_inside(a1prec, thpr0, thpr1).mask
                        

                        #-- Region --------------------------
                        [[lat0,lon0],[lat1,lon1]] = dBBox[region]

                        a1flaglat = ma.masked_inside(a1lat, lat0, lat1).mask
                        a1flaglon = ma.masked_inside(a1lon, lon0, lon1).mask

                        a1flagregion= a1flaglat * a1flaglon

 
                        #-- Screen --------------------------
                        a1flag = a1flagstype * a1flagptype * a1flagph * a1flagp * a1flagregion

                        if a1flag.sum()==0:
                            continue 

                        a2varTmp = a2var[a1flag,:]

                        if a2prof is None:
                            a2prof = a2varTmp
                        else:
                            a2prof = np.concatenate([a2prof, a2varTmp],axis=0)

                        print oid,a2prof.shape

                #**********************************************
                # Save
                #**********************************************
                stype,ptype, ph, prrange = key
                thpr0,thpr1 = prrange

                outDir = '/home/utsumi/temp/ret/tempdat'
 
                util.mk_dir(outDir)
                stamp  = '%s-s-%s.p-%s.ph-%s.pr-%.1f-%.1f.%s.%04d%02d'%(var, stype,ptype,ph,thpr0,thpr1,region,Year,Mon)
                stackPath    = outDir  + '/%s.stack.%s.npy' %(var,stamp)

                if calcflag == True:
                    np.save(stackPath, a2prof)
                    print stackPath

                #**********************************************
                # Load
                #**********************************************
                a2prof = np.load(stackPath) 


                a2prof = ma.masked_invalid(a2prof)
                a2prof = ma.masked_less(a2prof, 0)
                a1prof = a2prof.mean(axis=0)
                a1sum = a2prof.sum(axis=0)
                a1num = a2prof.count(axis=0)

                print a2prof
                print ''
                print a2prof.min(),a2prof.max()
                #**********************************************
                # mean profile
                #**********************************************
                ymin, ymax = 0, 12

                fig = plt.figure(figsize=(2.5,3.2))
                ax  = fig.add_axes([0.25,0.15,0.65,0.7])

                a1y = 0.25 + np.arange(nz) * 0.5 # [km]

                ax.plot( a1prof, a1y, '-', c='k', linewidth=2, label='CMB')

                ax.set_ylim([ymin,ymax])
                #ax.set_xlim([0,None])

                stitle = '%s %s %s %s %s'%(var, region, stype, ptype, ph)
                stitle = stitle + '\n'+'frompair %.1f-%.1fmm/h'%(thpr0,thpr1)
                plt.title(stitle, fontsize=9)
                #plt.legend()


                plt.xlabel('precipitation water (g/m3)')
                plt.ylabel('hight (above sea level) (km)')
                figPath = outDir + '/prof.drct.%s.%s.png'%(stamp,region)
                plt.savefig(figPath)
                print figPath


                #**********************************************
                # Sum profile
                #**********************************************
                ymin, ymax = 0, 12

                fig = plt.figure(figsize=(2.5,3.2))
                ax  = fig.add_axes([0.25,0.15,0.65,0.7])

                a1y = 0.25 + np.arange(nz) * 0.5 # [km]

                ax.plot( a1sum, a1y, '-', c='k', linewidth=2, label='CMB')

                ax.set_ylim([ymin,ymax])
                #ax.set_xlim([0,None])

                stitle = 'Sum %s %s %s %s'%(region, stype, ptype, ph)
                stitle = stitle + '\n'+'%.1f-%.1fmm/h'%(thpr0,thpr1)
                plt.title(stitle, fontsize=9)
                #plt.legend()


                plt.xlabel('Sum (g/m3)')
                plt.ylabel('hight (above sea level) (km)')
                figPath = outDir + '/sum.drct.%s.%s.png'%(stamp,region)
                plt.savefig(figPath)
                print figPath


                #**********************************************
                # Sum profile
                #**********************************************
                ymin, ymax = 0, 12

                fig = plt.figure(figsize=(2.5,3.2))
                ax  = fig.add_axes([0.25,0.15,0.65,0.7])

                a1y = 0.25 + np.arange(nz) * 0.5 # [km]

                ax.plot( a1num, a1y, '-', c='k', linewidth=2, label='CMB')

                ax.set_ylim([ymin,ymax])
                #ax.set_xlim([0,None])

                stitle = 'Num %s %s %s %s'%(region, stype, ptype, ph)
                stitle = stitle + '\n'+'%.1f-%.1fmm/h'%(thpr0,thpr1)
                plt.title(stitle, fontsize=9)
                #plt.legend()


                plt.xlabel('Num ')
                plt.ylabel('hight (above sea level) (km)')
                figPath = outDir + '/num.drct.%s.%s.png'%(stamp,region)
                plt.savefig(figPath)
                print figPath



