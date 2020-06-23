from numpy import *
import numpy as np
import h5py
import myfunc.util as util
import glob
import gzip
from collections import deque
from datetime import datetime, timedelta
import os, sys, socket

iDTime = datetime(2018,1,24)
eDTime = datetime(2018,1,24)

#----------------------------
gpr_amsr2      = ["GCOMW1","AMSR2","2A-CLIM","gprof","V05"]
gpr_ssmis_f16  = ["F16","SSMIS","2A-CLIM","gprof","V05"]
gpr_atms_noaa20= ["NOAA20","ATMS","2A-CLIM","gprof","V05"]
gpr_mhs_metopa = ["METOPA","MHS","2A-CLIM","gprof","V05"]

gmi       = ["GPM","GMI","1C","1C","V05"]
amsr2     = ["GCOMW1","AMSR2","1C","1C","V05"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05"]
ssmis_f18 = ["F18","SSMIS","1C","1C","V05"]
atms_npp  = ["NPP","ATMS","1C","1C","V05"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05"]

mhs_metopa= ["METOPA","MHS","1C","1C","V05"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05"]
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05"]

#lspec = [amsr2, ssmis_f16, ssmis_f17, ssmis_f18, atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
lspec = [ssmis_f16]

dnscan = {'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}



#*****************************************
# Start sate,sensor loop
#-----------------------------------------
for spec in lspec:
    print 'spec=',spec
    print ''
    sate      = spec[0]
    sensor    = spec[1]
    prdName   = spec[2]
    prj       = spec[3]
    ver       = spec[4]

    #*****************************************
    myhost = socket.gethostname()
    if myhost =="well":
        workDir   = '/home/utsumi/mnt/lab_work'
        tankDir   = '/home/utsumi/mnt/lab_tank'
        tbbaseDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
    else:
        workDir   = '/work'
        tankDir   = '/tank'
        tbbaseDir = '/work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
    #-----------------------------------------

    iY,iM = iDTime.timetuple()[:2]
    eY,eM = eDTime.timetuple()[:2]
    lYM   = util.ret_lYM([iY,iM],[eY,eM])

    ltbPathAll = []
    liescanAll = []
    for (Year,Mon) in lYM:
        listPath = tankDir + '/utsumi/PMM/US/obtlist/overpass.%s.%s.%04d.%02d.csv'%(sensor,sate,Year,Mon)
        f=open(listPath,'r'); lines=f.readlines(); f.close()
        for line in lines:
            _,_,Day,oid,iscan,escan = map(int,line.strip().split(','))

            #if (sensor=='AMSR2')&(oid <= 30345): continue # test
            if (iDTime<=datetime(Year,Mon,Day))&(datetime(Year,Mon,Day)<=eDTime):
                ltbPathTmp = sorted(glob.glob(tbbaseDir + '/%04d/%02d/%02d/*.%06d.????.HDF5'%(Year,Mon,Day,oid)))
                ltbPathAll = ltbPathAll + ltbPathTmp
                liescanAll.append([Year,Mon,Day,oid,iscan,escan])

    #*****************
    # Start oid loop
    #*****************
    for (tbPath, iescan) in zip(ltbPathAll, liescanAll):
        Year,Mon,Day,oid,iscan,escan = iescan
        DTime = datetime(Year,Mon,Day)
        #if DTime != datetime(2015,1,9): continue  # test

        #*****************
        # Read L1C
        #*****************
        with h5py.File(tbPath,'r') as h:
            a2latpmw = h['S1/Latitude'][iscan:escan+1]
            a2lonpmw = h['S1/Longitude'][iscan:escan+1]
            a1yyyy = h['S1/ScanTime/Year'][iscan:escan+1]
            a1mm   = h['S1/ScanTime/Month'][iscan:escan+1]
            a1dd   = h['S1/ScanTime/DayOfMonth'][iscan:escan+1]
            a1hh   = h['S1/ScanTime/Hour'][iscan:escan+1]
            a1mn   = h['S1/ScanTime/Minute'][iscan:escan+1]
            a1ss   = h['S1/ScanTime/Second'][iscan:escan+1]

        a1dtime = [datetime(yyyy,mm,dd,hh,mn,ss) for (yyyy,mm,dd,hh,mn,ss) in zip(a1yyyy,a1mm,a1dd,a1hh,a1mn,a1ss)]

        nypmw,nxpmw = a2latpmw.shape

        a1latpmw = a2latpmw.flatten()
        a1lonpmw = a2lonpmw.flatten()

        #*****************
        # Make MRMS list
        #*****************
        mrmsDir = workDir + '/hk02/PMM/MRMS/level2/%s/%04d/%02d'%(sate,Year,Mon)
        ssearch = mrmsDir + '/PRECIPRATE.GC.????????.??????.%05d.dat.gz'%(oid)
        lmrmsPath = sort(glob.glob(ssearch))
        if len(lmrmsPath)==0:
            print 'No MRMS',Year,Mon,Day,oid

        #*****************
        # Initialize output    
        #*****************
        a2sum = zeros([nypmw,nxpmw],float64)
        a2num = zeros([nypmw,nxpmw],int32)

        #*****************
        # Read MRMS
        #*****************
        ''' first row = northern end '''
        ''' Header
        ncols 7000
        nrows 3500
        xllcenter -129.995000
        yllcenter 20.005000
        cellsize 0.010000
        '''

        for mrmsPath in lmrmsPath:    
            hhmnss = os.path.basename(mrmsPath).split('.')[3]
            slabel = '.'.join(os.path.basename(mrmsPath).split('.')[2:5])
            rqiPath = mrmsDir + '/RQI.%s.asc.gz'%(slabel) 
            maskPath= mrmsDir + '/MASK.%s.asc.gz'%(slabel)

            a2mrms = []
            a2rqi  = []
            a2mask = []

            with gzip.open(mrmsPath, 'rt') as f:
                lines = f.readlines()
                for line in lines[6:]:
                    a2mrms.append(map(float,line.split()))
            a2mrms = np.flipud(np.array(a2mrms))

            sys.exit()


            #*****************
            # Project MRMS over GMI track 
            #*****************
            for i in range(len(a1lon)):
                lat = a1lat[i]
                lon = a1lon[i]
                dat = a1dat[i]
                #lat = 46.291878  # test
                #lon = -129.8957  # test

                RADEARTH = 6371
                DTR      = 0.017453
                adist = RADEARTH*np.arccos(np.cos(DTR*a1lonpmw-DTR*lon)*np.cos(DTR*a1latpmw)*np.cos(DTR*lat) + np.sin(DTR*a1latpmw)*np.sin(DTR*lat))

                j = adist.argmin()
                #a1out[j] = dat
                a1sum[j] = a1sum[j] + dat
                a1num[j] = a1num[j] + 1

        a1out = (ma.masked_where(a1num==0, a1sum)/a1num).filled(-9999.)
        ##*****************
        ## Save
        ##*****************
        #a2out = a1out.reshape(nypmw,nxpmw).astype(float32)
        #outDir= '/work/hk01/PMM/MRMS/match-GMI-orbit'
        #util.mk_dir(outDir)
        #outPath= outDir + '/GMI.MRMS.130W_55W_20N_55N.%04d%02d%02d.%06d.%05d-%05d.npy'%(Year,Mon,Day,oid,iy,ey)
        #np.save(outPath,a2out)
        #print outPath
   
     
