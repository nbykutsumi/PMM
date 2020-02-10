import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import numpy as np
import myfunc.util as util
import glob
from datetime import datetime, timedelta
import sys, os
from numpy import *

#calcflag= True
calcflag= False
iDTime = datetime(2014,7,1)
eDTime = datetime(2014,7,31)
#eDTime = datetime(2014,7,2)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

thpr = 0.5
#lstype = ['sea','land']
lstype = ['land']
lptype = ['conv']
#lregion = ['TRO']
lregion = ['MIDN']

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

nz = 88
ix, ex = 14, 34  # center xpy=24
#ix, ex = 23,25  # center xpy=24

for region in lregion:
    [[lat0,lon0],[lat1,lon1]] = dBBox[region]
    for stype in lstype:
        for ptype in lptype:


            a2wat = None
            for DTime in lDTime:
                if calcflag==False: continue
    
                Year,Mon,Day = DTime.timetuple()[:3]
                srcDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
                ssearch = srcDir + '/2B.GPM.DPRGMI.CORRA2018.20140701-S030429-E043702.001921.V06A.HDF5'
                ssearch = srcDir + '/2B.GPM.DPRGMI.*.V06A.HDF5'
                lsrcPath= np.sort(glob.glob(ssearch))
            
                for srcPath in lsrcPath:
                    with h5py.File(srcPath,'r') as h:
                        a2elev = h['NS/Input/surfaceElevation'][:,ix:ex+1]
                        a2surf = h['NS/Input/surfaceType'][:,ix:ex+1]
                        a2ptype= h['NS/Input/precipitationType'][:,ix:ex+1]
                        a2lat  = h['NS/Latitude'][:,ix:ex+1]
                        a2lon  = h['NS/Longitude'][:,ix:ex+1]
                        a3wat  = h['NS/precipTotWaterCont'][:,ix:ex+1]
                        a2prec = h['NS/surfPrecipTotRate'][:,ix:ex+1]
            
                    #**************************
                    # Screening
                    #--------------------------
    
                    # Precipitation rate
                    a2flagprec= ma.masked_greater_equal(a2prec, thpr).mask
    
                    # Region
                    a2flaglat = ma.masked_inside(a2lat,lat0,lat1).mask
                    a2flaglon = ma.masked_inside(a2lon,lon0,lon1).mask
                    a2flagregion = a2flaglat * a2flaglon        
    
                    # Surface type
                    if stype =='sea':
                        a2flagsurf = ma.masked_inside(a2surf,0,99).mask
                    elif stype=='land':
                        a2flagsurf = ma.masked_inside(a2surf,100,199).mask
                    else:
                        print 'check stype',stype
                        sys.exit()
    
                    # Precipitaiton type
                    a2ptype = (a2ptype/10000000.).astype('int16')
                    if ptype=='stra': 
                        a2flagptype = ma.masked_equal(a2ptype,1).mask
                    elif ptype=='conv': 
                        a2flagptype = ma.masked_equal(a2ptype,2).mask
                    else:
                        print 'check ptype',ptype
                        sys.exit()

                    #-- Flag --------
                    a2flag = a2flagprec * a2flagregion * a2flagsurf * a2flagptype
    
                    a2watTmp= a3wat[a2flag]
    
                    if a2wat is None:
                        a2wat = a2watTmp
                    else:
                        a2wat = np.concatenate([a2wat, a2watTmp], axis=0)
   
                    print DTime,a2wat.shape, a2flagprec.sum(), a2flagregion.sum(), a2flagsurf.sum(), a2flagptype.sum()
            #--- Save --- 
            outDir = '/home/utsumi/temp/ret/tempdat'
            outPath= outDir + '/prof.hdf.x-%02d-%02d.%s.%s.%s.npy'%(ix,ex,stype,ptype,region)
            if calcflag==True:
                np.save(outPath, a2wat)
                print outPath
                print a2wat.shape


            ymin,ymax = [0,12]
            a2wat = np.load(outPath)
            a1wat = ma.masked_less(a2wat,0).mean(axis=0)
            a1num = ma.masked_less(a2wat,0).count(axis=0)
            a1wat = a1wat[::-1]  # bottom to top
            a1num = a1num[::-1] 

            #--- Figure (profile) ---
            fig = plt.figure(figsize=(2.5,3.2))
            ax  = fig.add_axes([0.25,0.15,0.65,0.7])

            #a1y = 0.25 + np.arange(nz) * 0.5 # [km]
            a1y = 0.125 + np.arange(nz) * 0.25 # [km]

            ax.plot( a1wat, a1y, '-', c='k', linewidth=2, label='CMB')

            ax.set_ylim([ymin,ymax])
            #ax.set_xlim([0,None])
            ax.set_xlim([0,1.2])

            stitle = 'WatCont %s %s %s %s'%('rad', region, stype, ptype)
            stitle = stitle + '\n'+'from HDF x=%d-%d'%(ix,ex)
            plt.title(stitle, fontsize=9)
            outPath= outDir + '/prof.hdf.x-%02d-%02d.%s.%s.%s.png'%(ix,ex,stype,ptype,region)
            plt.savefig(outPath)
            print outPath
            plt.clf()
            #plt.legend()

            
            #--- Figure (profile) ---
            fig = plt.figure(figsize=(2.5,3.2))
            ax  = fig.add_axes([0.25,0.15,0.65,0.7])

            #a1y = 0.25 + np.arange(nz) * 0.5 # [km]
            a1y = 0.125 + np.arange(nz) * 0.25 # [km]

            ax.plot( a1num, a1y, '-', c='k', linewidth=2, label='CMB')

            ax.set_ylim([ymin,ymax])
            ax.set_xlim([0,None])

            stitle = 'Num %s %s %s %s'%('rad', region, stype, ptype)
            stitle = stitle + '\n'+'from HDF'
            plt.title(stitle, fontsize=9)
            outPath= outDir + '/num.hdf.x-%02d-%02d.%s.%s.%s.png'%(ix,ex,stype,ptype,region)
            plt.savefig(outPath)
            print outPath
            plt.clf()
            #plt.legend()

            





