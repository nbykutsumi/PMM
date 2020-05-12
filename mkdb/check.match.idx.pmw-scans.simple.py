from numpy import *
import myfunc.util as util
import myfunc.IO.GPM.l1_gmi as l1_gmi
from f_match_fov import *
import sys, os, glob
from datetime import datetime, timedelta
import numpy as np
import h5py

'''  Simple matchiup '''

iDTime = datetime(2018,1,1)
eDTime = datetime(2018,1,1)
#iDTime = datetime(2014,9,12)
#eDTime = datetime(2015,11,30)


dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)


#lDTime = [DTime for DTime in lDTime if not ((datetime(2017,9,26)<=DTime)&(datetime(2017,9,29))]
lDTime = [DTime for DTime in lDTime if not (datetime(2017,9,26)<=DTime)&(DTime<=datetime(2017,9,29))]

mainscan  = 'high'
mainscan  = 'low'
amsr2     = ["GCOMW1","AMSR2"]
ssmis_f16 = ["F16","SSMIS"]
ssmis_f17 = ["F17","SSMIS"]
ssmis_f18 = ["F18","SSMIS"]
atms_npp  = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]
mhs_metopa= ["METOPA","MHS"]
mhs_metopb= ["METOPB","MHS"]
mhs_noaa18= ["NOAA18","MHS"]
mhs_noaa19= ["NOAA19","MHS"]

lspec = [ssmis_f16]
#lspec = [amsr2]
version = '05'

dnx = {'AMSR2':[243,486], 'SSMIS':[90,180]}

#****************************
def calc_dist(lat1, lon1, lat2, lon2):
    RADEARTH = 6371
    DTR = 0.017453

    #RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2))

    return RADEARTH*np.arccos(np.cos(DTR*lon1-DTR*lon2)*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))

def pickyx_low2high_amsr2(ny_high):
    nxhi  = 486
    axpick = ((np.arange(nxhi) -1)/2).astype('int16')
    axpick[0] = 0    # [0,0,0,1,1,2,2,3,3,...]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_low2high_ssmis(ny_high):
    nxhi = 180
    axpick = (np.arange(nxhi) /2).astype('int16')  # [0,0,1,1,2,2,3,3,...,89,89]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_amsr2(ny_high):
    nxlw  = 243
    axpick = (np.arange(nxlw)*2+2).astype('int16')
    axpick[-1] = nxlw*2-1   # [2,4,6,...482,484,485]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_ssmis(ny_high):
    nxlw = 90
    axpick = (np.arange(nxlw)*2+1).astype('int16') # [1,3,5,...,177,179]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

#****************************
for (sate,sensor) in lspec:
    baseDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/1C/V%s'%(sate, sensor, version)

    #obaseDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.%s.%s.V%s/%s.ABp%03d-%03d.smpl'%(ver)
    #obaseDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.%s.%s.V%s/%s.ABp%03d-%03d.smpl'%(ver)

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]   
        print DTime 

        srcDir = baseDir +'/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch    = srcDir + '/1C.%s.%s.*.HDF5'%(sate,sensor)
        lsrcPath   = sort(glob.glob(ssearch))

        if len(lsrcPath)==0:
            print 'No 1C file',Year,Mon,Day
            print ssearch
            continue

        for srcPath in lsrcPath:
            print srcPath
            oid       = int(srcPath.split('.')[-3])
            if sensor=='AMSR2':
                if oid != 29925: continue  # AMSR2
            elif sensor=='SSMIS':
                if oid != 73300: continue  # F16 SSMIS

            with h5py.File(srcPath,'r') as h:
                if sensor=='AMSR2':
                    a2lat1 = h['S1/Latitude'][:]
                    a2lon1 = h['S1/Longitude'][:]
                    a2lat2 = h['S2/Latitude'][:]                 
                    a2lon2 = h['S2/Longitude'][:]
                    a2lat3 = h['S3/Latitude'][:]                  
                    a2lon3 = h['S3/Longitude'][:]
                    a2lat4 = h['S4/Latitude'][:]                    
                    a2lon4 = h['S4/Longitude'][:]
                    a2lat5 = h['S5/Latitude'][:]                    
                    a2lon5 = h['S5/Longitude'][:]

                elif sensor=='SSMIS':
                    a2lat1 = h['S1/Latitude'][:]
                    a2lon1 = h['S1/Longitude'][:]
                    a2lat2 = h['S2/Latitude'][:]                 
                    a2lon2 = h['S2/Longitude'][:]
                    a2lat3 = h['S3/Latitude'][:]                  
                    a2lon3 = h['S3/Longitude'][:]
                    a2lat4 = h['S4/Latitude'][:]                    
                    a2lon4 = h['S4/Longitude'][:]
                    print '---------in'
            print sensor    

            if mainscan=='high':
                if sensor=='AMSR2':
                    nxlw = dnx[sensor][0]
                    nyhi,nxhi = a2lon5.shape
                    a2ypick, a2xpick = pickyx_low2high_amsr2(nyhi)

                    o2lon1 = a2lon1[a2ypick, a2xpick]
                    o2lon2 = a2lon2[a2ypick, a2xpick]
                    o2lon3 = a2lon3[a2ypick, a2xpick]
                    o2lon4 = a2lon4[a2ypick, a2xpick]
                    o2lat1 = a2lat1[a2ypick, a2xpick]
                    o2lat2 = a2lat2[a2ypick, a2xpick]
                    o2lat3 = a2lat3[a2ypick, a2xpick]
                    o2lat4 = a2lat4[a2ypick, a2xpick]

                    o2lat5 = a2lat5
                    o2lon5 = a2lon5

                    a2dist1 = calc_dist( o2lat1, o2lon1, o2lat5, o2lon5)
                    a2dist2 = calc_dist( o2lat2, o2lon2, o2lat5, o2lon5)
                    a2dist3 = calc_dist( o2lat3, o2lon3, o2lat5, o2lon5)
                    a2dist4 = calc_dist( o2lat4, o2lon4, o2lat5, o2lon5)

                elif sensor =='SSMIS':
                    nxlw = dnx[sensor][0]
                    nyhi,nxhi = a2lon4.shape
                    a2ypick, a2xpick = pickyx_low2high_ssmis(nyhi)

                    o2lon1 = a2lon1[a2ypick, a2xpick]
                    o2lon2 = a2lon2[a2ypick, a2xpick]
                    o2lat1 = a2lat1[a2ypick, a2xpick]
                    o2lat2 = a2lat2[a2ypick, a2xpick]

                    o2lat3 = a2lat3
                    o2lon3 = a2lon3
                    o2lat4 = a2lat4
                    o2lon4 = a2lon4


                    a2dist1 = calc_dist( o2lat1, o2lon1, o2lat4, o2lon4)
                    a2dist2 = calc_dist( o2lat2, o2lon2, o2lat4, o2lon4)
                    a2dist3 = calc_dist( o2lat3, o2lon3, o2lat4, o2lon4)

            elif mainscan =='low':
                if sensor=='AMSR2':
                    nxlw = dnx[sensor][0]
                    nyhi,nxhi = a2lon5.shape
                    a2ypick, a2xpick = pickyx_high2low_amsr2(nyhi)

                    o2lon1 = a2lon1
                    o2lon2 = a2lon2
                    o2lon3 = a2lon3
                    o2lon4 = a2lon4
                    o2lat1 = a2lat1
                    o2lat2 = a2lat2
                    o2lat3 = a2lat3
                    o2lat4 = a2lat4

                    o2lat5 = a2lat5[a2ypick, a2xpick]
                    o2lon5 = a2lon5[a2ypick, a2xpick]

                    a2dist2 = calc_dist( o2lat2, o2lon2, o2lat1, o2lon1)
                    a2dist3 = calc_dist( o2lat3, o2lon3, o2lat1, o2lon1)
                    a2dist4 = calc_dist( o2lat4, o2lon4, o2lat1, o2lon1)
                    a2dist5 = calc_dist( o2lat5, o2lon5, o2lat1, o2lon1)

                elif sensor =='SSMIS':
                    nxlw = dnx[sensor][0]
                    nyhi,nxhi = a2lon4.shape
                    a2ypick, a2xpick = pickyx_high2low_ssmis(nyhi)

                    o2lon1 = a2lon1
                    o2lon2 = a2lon2
                    o2lat1 = a2lat1
                    o2lat2 = a2lat2

                    o2lat3 = a2lat3[a2ypick, a2xpick]
                    o2lon3 = a2lon3[a2ypick, a2xpick]
                    o2lat4 = a2lat4[a2ypick, a2xpick]
                    o2lon4 = a2lon4[a2ypick, a2xpick]

                    a2dist2 = calc_dist( o2lat2, o2lon2, o2lat1, o2lon1)
                    a2dist3 = calc_dist( o2lat3, o2lon3, o2lat1, o2lon1)
                    a2dist4 = calc_dist( o2lat4, o2lon4, o2lat1, o2lon1)


            ##---- test ----------
            nx0 = nxlw
            #x1  = int(0.9*nx0/2)-1
            x1  = nx0-1
            y1  = int(0.8*0.25*nyhi)
            #y1  = nyhi-1
            print 'y1,x1=',y1,x1
            #print 5,o2lat5[y0,x0], o2lon5[y0,x0]
            print 1,o2lat1[y1,x1], o2lon1[y1,x1]
            print 2,o2lat2[y1,x1], o2lon2[y1,x1], 'dist=',a2dist2[y1,x1]
            #print 3,o2lat3[y1,x1], o2lon3[y1,x1], 'dist=',a2dist3[y1,x1]
            print 4,o2lat4[y1,x1], o2lon4[y1,x1], 'dist=',a2dist4[y1,x1]
            #print 5,o2lat5[y1,x1], o2lon5[y1,x1], 'dist=',a2dist5[y1,x1]

            sys.exit()



            ##---- test ----------
            #nx0 = nxlw
            ##x1  = int(0.9*nx0/2)*2 -2
            #x1  = 4
            ##y1  = int(0.8*0.25*nyhi)
            #y1  = nyhi-1
            #print 'y1,x1=',y1,x1
            ##print 5,o2lat5[y0,x0], o2lon5[y0,x0]
            #print 4,o2lat4[y1,x1], o2lon4[y1,x1]
            #print 1,o2lat1[y1,x1], o2lon1[y1,x1], 'dist=',a2dist1[y1,x1]

            #print 2,o2lat2[y1,x1], o2lon2[y1,x1], 'dist=',a2dist2[y1,x1]

            ##print 3,o2lat3[y1,x1], o2lon3[y1,x1], 'dist=',a2dist3[y1,x1]
            ##print 4,o2lat4[y1,x1], o2lon4[y1,x1], 'dist=',a2dist4[y1,x1]

            #sys.exit()


