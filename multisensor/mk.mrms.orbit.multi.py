# %%
import matplotlib.pyplot as plt
#%matplotlib inline
from numpy import *
import numpy as np
import h5py
import myfunc.util as util
import glob
import gzip
from collections import deque
from datetime import datetime, timedelta
import os, sys, socket
from bisect import bisect_left

#iDTime = datetime(2018,3,1)
#eDTime = datetime(2018,12,31)

iDTime = datetime(2018,1,24)
eDTime = datetime(2018,12,31)

#iDTime = datetime(2018,1,23)
#eDTime = datetime(2018,2,28)




#----------------------------
gmi       = ["GPM","GMI","1C","1C","V05"]
amsr2     = ["GCOMW1","AMSR2","1C","1C","V05"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05"]      # 37V is missing after April 2016
ssmis_f18 = ["F18","SSMIS","1C","1C","V05"]
atms_npp  = ["NPP","ATMS","1C","1C","V05"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05"]   # MRMS is not available

mhs_metopa= ["METOPA","MHS","1C","1C","V05"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05"]   # Not available at arthurhou.pps after 2018/10/21
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05"]

#lspec = [amsr2, ssmis_f16, ssmis_f17, ssmis_f18, atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19, gmi]
lspec = [mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19, gmi]
#lspec = [gmi]
#lspec = [ssmis_f16, ssmis_f18, atms_npp, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19, gmi]
#lspec = [mhs_noaa18]
#lspec = [gmi]

#scantype='epc'
scantype='gpr'

ny_mrms = 3500
nx_mrms = 7000

#*****************************************
def calc_km(lat1,lon1,lat2,lon2):
    RADEARTH = 6371
    DTR = 0.017453
    km= RADEARTH*np.arccos(np.cos(DTR*(lon1-lon2))*np.cos(DTR*lat1)*np.cos(DTR*lat2) + np.sin(DTR*lat1)*np.sin(DTR*lat2))
    return km

#*****************************************
# Start sate,sensor loop
#-----------------------------------------
for spec in lspec:
    print('spec=',spec)
    print('')
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
        mrmsbaseDir='/media/disk2/data/PMM/MRMS'
    else:
        workDir   = '/work'
        tankDir   = '/tank'
        tbbaseDir = '/work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
        mrmsbaseDir='/work/hk02/PMM/MRMS'

    #*****************************************
    # Make path list
    #-----------------------------------------
    iY,iM = iDTime.timetuple()[:2]
    eY,eM = eDTime.timetuple()[:2]
    lYM   = util.ret_lYM([iY,iM],[eY,eM])

    ltbPathAll = []
    liescanAll = []
    for (Year,Mon) in lYM:

        if (sate=='NOAA18')&(datetime(Year,Mon,1,0)>=datetime(2018,10,21)):continue

        listPath = tankDir + '/utsumi/PMM/US/obtlist/overpass.%s.%s.%04d.%02d.csv'%(sensor,sate,Year,Mon)
        f=open(listPath,'r'); lines=f.readlines(); f.close()
        for line in lines:
            _,_,Day,oid,iscan,escan = list(map(int,line.strip().split(',')))
            #print(line)
            #if ('sate'=='AMSR2')&(datetime(2018,8,31) < datetime(Year,Mon,Day)): continue  # No NOAA18 data after 2018/10/21

            #if (sate=='GCOMW1'): continue  # test
            #if (sate=='NPP')&(Mon<3): continue  # test
            #if (sate=='GPM')&(Mon<3): continue  # test
            #if (sate=='METOPA')&(Mon<3): continue  # test
            #if (sate=='METOPB')&(Mon<3): continue  # test
            #if (sate=='NOAA18')&(Mon<3): continue  # No NOAA18 data after 2018/10/21
            #if (sate=='NOAA19')&(Mon<3): continue  # 
            #if (sate=='F16'): continue  # test
            #if (sate=='F18')&(Mon<3): continue  # test
            #if (sate=='GPM')&(oid !=22370): continue  # test
            if (iDTime<=datetime(Year,Mon,Day))&(datetime(Year,Mon,Day)<=eDTime):
                ltbPathTmp = sorted(glob.glob(tbbaseDir + '/%04d/%02d/%02d/*.%06d.????.HDF5'%(Year,Mon,Day,oid)))
                ltbPathAll = ltbPathAll + ltbPathTmp
                liescanAll.append([Year,Mon,Day,oid,iscan,escan])

                #print(tbbaseDir + '/%04d/%02d/%02d/*.%06d.????.HDF5'%(Year,Mon,Day,oid))
    #*****************
    # Start oid loop
    #*****************
    for (tbPath, iescan) in zip(ltbPathAll, liescanAll):
        Year,Mon,Day,oid,iscan,escan = iescan
        DTime = datetime(Year,Mon,Day)

        #--- test ----------
        #if (sate=='GCOMW1')&(DTime <datetime(2018,7,21)):continue  #
        #if (sate=='F18')&(DTime <datetime(2018,7,13)):continue  # SSMIS F18
        #if (sate=='NPP')&(DTime <datetime(2018,6,10)):continue  # SSMIS F18
        #if (sate=='NPP')&(oid <34296):continue  # SSMIS F18
        if (sate=='NOAA18')&(DTime >datetime(2018,10,21)):continue  # No NOAA18 data after 2018/10/21
        #if (sate=='NOAA19')&(DTime <datetime(2018,6,18)):continue  #
        if (sate=='METOPA')&(oid<62004):continue
        #if (sate=='METOPB')&(DTime <datetime(2018,11,23)):continue
        #if (sate=='GPM')&(DTime <datetime(2018,2,13)):continue

        #if oid != 65343: continue  # test

        #*****************
        # Read L1C
        #*****************
        if scantype=='epc':
            mainscan = 1
        elif scantype=='gpr':
            mainscan={'GMI':1, 'AMSR2':5, 'SSMIS':1, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]

        scanname = 'S%d'%(mainscan)
        with h5py.File(tbPath,'r') as h:
            a2latpmw = h[scanname+'/Latitude'][iscan:escan+1].astype('float64')
            a2lonpmw = h[scanname+'/Longitude'][iscan:escan+1].astype('float64')
            a1yyyy = h[scanname+'/ScanTime/Year'][iscan:escan+1]
            a1mm   = h[scanname+'/ScanTime/Month'][iscan:escan+1]
            a1dd   = h[scanname+'/ScanTime/DayOfMonth'][iscan:escan+1]
            a1hh   = h[scanname+'/ScanTime/Hour'][iscan:escan+1]
            a1mn   = h[scanname+'/ScanTime/Minute'][iscan:escan+1]
            a1ss   = h[scanname+'/ScanTime/Second'][iscan:escan+1]
            a2inc  = h[scanname+'/incidenceAngle'][iscan:escan+1, :, 0]


        a1dtime = [datetime(yyyy,mm,dd,hh,mn,ss) for (yyyy,mm,dd,hh,mn,ss) in zip(a1yyyy,a1mm,a1dd,a1hh,a1mn,a1ss)]

        nypmw,nxpmw = a2latpmw.shape

        a1latpmw = a2latpmw.flatten()
        a1lonpmw = a2lonpmw.flatten()


        #*****************
        # PMW effective area
        #*****************
        dpixres       = {'GMI':15, 'AMSR2':15, 'SSMIS':25, 'ATMS':-9999, 'MHS':-9999}
        dpixres_nadir = {'GMI':-9999, 'AMSR2':-9999, 'SSMIS':-9999, 'ATMS':30, 'MHS':15}

        dcross_scan   = {'GMI':0, 'AMSR2':0, 'SSMIS':0, 'ATMS':1, 'MHS':1} 

        pixres       = dpixres[sensor]
        pixres_nadir = dpixres_nadir[sensor]
        a2inc = ma.masked_outside(a2inc,-90,90)
        if dcross_scan[sensor]==1:
            DTR = np.pi/180.
            a2pixres = pixres_nadir / np.cos(DTR*a2inc)
            a2pixara = np.square(a2pixres)
        else:
            a2pixara = np.square(np.full([nypmw,nxpmw], pixres))

        a1pixara = ma.masked_where(a2inc[0]<-90, a2pixara[0])   # Representative for the entire granule
        a1pixara = ma.masked_where(a2inc[0]>90, a2pixara[0])   # Representative for the entire granule

        #*****************
        # Determine the number of MRMS pixels
        #*****************
        ara_mrms = 1  # km2
        a1nrad = (np.sqrt(a1pixara/ara_mrms) /2).astype('int32')  

        #*****************
        # Determine the corresponding MRMS location
        #*****************
        lat0_mrms = 20
        lon0_mrms = -130
        lat1_mrms = 55
        lon1_mrms = -60
        dlatlon_mrms = 0.01

        a2latpmw = ma.masked_outside(a2latpmw, lat0_mrms, lat1_mrms)
        a2lonpmw = ma.masked_outside(a2lonpmw, lon0_mrms, lon1_mrms)
        a2masklatlon = a2latpmw.mask + a2lonpmw.mask        
        a2latpmw = ma.masked_where(a2masklatlon, a2latpmw)
        a2lonpmw = ma.masked_where(a2masklatlon, a2lonpmw)

        a2y_mrms = ((a2latpmw - lat0_mrms)/dlatlon_mrms).astype('int32')
        a2x_mrms = ((a2lonpmw - lon0_mrms)/dlatlon_mrms).astype('int32')
         
        #*****************
        # Make MRMS list
        #*****************
        mrmsDir = mrmsbaseDir + '/level2/%s/%04d/%02d'%(sate,Year,Mon)
        ssearch = mrmsDir + '/PRECIPRATE.GC.????????.??????.%05d.dat.gz'%(oid)
        lmrmsPath = sort(glob.glob(ssearch))
        if len(lmrmsPath)==0:
            print('No MRMS',Year,Mon,Day,oid)

        #*****************
        # Initialize output    
        #*****************
        a2ave = np.full([nypmw,nxpmw], -9999., dtype='float32')
        a2num = np.full([nypmw,nxpmw], -9999., dtype='float32')
        a2qlt = np.full([nypmw,nxpmw], -9999., dtype='float32')

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

        a2ave_mrms = np.full([nypmw,nxpmw], -9999., 'float32')
        a2goodfrac = np.full([nypmw,nxpmw], -9999., 'float32')

        a2num_mrms = np.full([nypmw,nxpmw], -9999., 'int16')

        lptype = [0,1,2,3,4,6,7,10,91,96]
        d2vfc_type = {}   # volume fraction
        d2num_type = {}
        for ptype in lptype:
            d2vfc_type[ptype] = np.full([nypmw,nxpmw], -9999., 'float32')
            d2num_type[ptype] = np.full([nypmw,nxpmw], -9999., 'int16')

        for mrmsPath in lmrmsPath:    
            slabel = '.'.join(os.path.basename(mrmsPath).split('.')[2:5])
            rqiPath = mrmsDir + '/RQI.%s.asc.gz'%(slabel) 
            typePath= mrmsDir + '/MASK.%s.asc.gz'%(slabel)

            yyyymmdd=os.path.basename(mrmsPath).split('.')[2]
            hhmnss = os.path.basename(mrmsPath).split('.')[3]
            dtimegv = datetime.strptime(yyyymmdd+hhmnss, '%Y%m%d%H%M%S')

            # MRMS file contains the backward 2-min average
            # Confirmed with Pierre.
            idtimegv= dtimegv - timedelta(seconds=60*2)
            edtimegv= dtimegv

            if (idtimegv > a1dtime[-1]): continue
            if (edtimegv < a1dtime[0]): continue
            print(slabel, dtimegv)
            print(mrmsPath)

            try:
                with gzip.open(mrmsPath, 'rt') as f:
                    a2mrms=np.array(f.read().split()[12:], 'float64').reshape(ny_mrms,nx_mrms)
    
                with gzip.open(rqiPath, 'rt') as f:
                    a2rqi =np.array(f.read().split()[12:], 'float32').reshape(ny_mrms,nx_mrms)

                with gzip.open(typePath, 'rt') as f:
                    a2type =np.array(f.read().split()[12:], 'float32').reshape(ny_mrms,nx_mrms)

            except:
                print('')
                print('Error')
                print('Skip')
                print('')
                continue


            #iy=561
            #ix=719

            #plt.figure()
            #plt.imshow(ma.masked_less_equal(a2mrms[iy-10:iy+10,ix-10:ix+10],0))
            #plt.colorbar()
            #plt.show()

            #plt.figure()
            #plt.imshow(ma.masked_less_equal(a2type[iy-10:iy+10,ix-10:ix+10],0))
            #plt.colorbar()
            #plt.show()

            #if a2mrms.max()>0:
            #    a = ma.masked_where(a2type<1,a2mrms)
            #    a = ma.masked_where(a2rqi<1, a)
            #    print a.min(), a.max()
            #    sys.exit()
            ##--- test -------
            ##a1latmrms= np.arange(20+0.005, 55-0.005+0.001, 0.01)
            ##a1lonmrms= np.arange(-130+0.005, -60-0.005+0.001, 0.01)
            ##a2lonmrms, a2latmrms = np.meshgrid(a1lonmrms, a1latmrms)
            #a1lonmrms = np.linspace(0,10,7000)
            #a1latmrms = np.linspace(0,10,3500)
            #a2lonmrms, a2latmrms = np.meshgrid(a1lonmrms, a1latmrms)
            #a2mrms    = a2lonmrms
            #a2rqi = ma.ones([ny_mrms,nx_mrms])

            ##----------------

            a2mrms = np.flipud(np.array(a2mrms))
            a2rqi  = np.flipud(np.array(a2rqi ))
            a2type = np.flipud(np.array(a2type))

            #*****************
            # Padding outer area of MRMS
            #*****************
            dpad = a1nrad[0] + 1  # take the largest rad

            a2mrms_pad = np.full([ny_mrms+dpad*2,nx_mrms+dpad*2],-9999., dtype='float64')
            a2mrms_pad[dpad:dpad+ny_mrms, dpad:dpad+nx_mrms] = a2mrms 

            #*****************
            # Partial extraction of PMW based on time
            #*****************
            iypart = bisect_left(a1dtime, idtimegv) 
            eypart = bisect_left(a1dtime, edtimegv) 
            nypart = eypart - iypart
            if nypart ==0: continue

            a2y_mrms_tmp = a2y_mrms[iypart:eypart] 
            a2x_mrms_tmp = a2x_mrms[iypart:eypart] 

            nypmw_tmp = eypart - iypart+1
            nxpmw_tmp = nxpmw

            #print 'iypart,eypart=',iypart,eypart
            #*****************
            # Process each angle bin
            #*****************
            for ix in range(nxpmw):
                a1x_mrms = a2x_mrms_tmp[:,ix]  # Keeps masks. Should be handled later.
                a1y_mrms = a2y_mrms_tmp[:,ix]


                nrad = a1nrad[ix]

                a2offsetx, a2offsety = np.meshgrid(list(range(-nrad,nrad+1)), list(range(-nrad,nrad+1)))
                a3collectx = np.array([a2offsetx + x for x in a1x_mrms], 'int32')
                a3collecty = np.array([a2offsety + y for y in a1y_mrms], 'int32')
                #a3collectx = np.apply_along_axis(add, 1, a1x_mrms[:,None], a2offsetx).astype('int32')
                #a3collecty = np.apply_along_axis(add, 1, a1y_mrms[:,None], a2offsety).astype('int32')

                #- If center MRMS pixel is out of boundary then replace ---
                a1xy_mrms_mask = a1x_mrms.mask + a1y_mrms.mask
                a3collectx[a1xy_mrms_mask,:,:] = -9999
                a3collecty[a1xy_mrms_mask,:,:] = -9999
                #----------------------------------------------------------
                a1collectx = a3collectx.flatten()
                a1collecty = a3collecty.flatten()

                #- If MRMS pixel location is out of boundary, then mask and replace --
                a3collectxy_mask = ma.masked_outside(a3collectx, 0, nx_mrms-1).mask + ma.masked_outside(a3collecty, 0, ny_mrms-1).mask   # masks are kept with flatten
                a1collectxy_mask = ma.masked_outside(a1collectx, 0, nx_mrms-1).mask + ma.masked_outside(a1collecty, 0, ny_mrms-1).mask   # masks are kept with flatten

                # temporarily replace with zero
                a1collectx[a1collectxy_mask] = 0 
                a1collecty[a1collectxy_mask] = 0

                #- Extract MRMS data ----------------
                a1collect_prec = a2mrms[a1collecty, a1collectx]
                a1collect_rqi  = a2rqi [a1collecty, a1collectx]
                a1collect_type = a2type[a1collecty, a1collectx]

                a3collect_prec = a1collect_prec.reshape(-1,nrad*2+1, nrad*2+1) 
                a3collect_rqi  = a1collect_rqi .reshape(-1,nrad*2+1, nrad*2+1) 
                a3collect_type = a1collect_type.reshape(-1,nrad*2+1, nrad*2+1) 

                #- Mask of invalid MRMS data -----------
                a3collect_mask = ma.masked_less(a3collect_rqi, 0.8).mask + ma.masked_less(a3collect_prec, 0).mask 

                #- Add mask of invalid pixel location --
                a3collect_mask = a3collect_mask + a3collectxy_mask

                if type(a3collect_mask) == np.bool_:
                    a3collect_mask = np.full(a3collect_prec.shape, a3collect_mask)

                a1goodfrac = ((nrad*2+1)**2 - a3collect_mask.sum(axis=(1,2))) / float( (nrad*2+1)**2)

                a3collect_prec = ma.masked_where(a3collect_mask, a3collect_prec)
                a3collect_type = ma.masked_where(a3collect_mask, a3collect_type)

                #a2ave_mrms[iypart:eypart,ix] = ma.masked_where(a3collect_mask, a3collect_prec).mean(axis=(1,2)).filled(-9999.)
                a2ave_mrms[iypart:eypart,ix] = a3collect_prec.mean(axis=(1,2)).filled(-9999.)
                a2goodfrac[iypart:eypart,ix] = a1goodfrac

                #-- precipitation type info ---
                a1sum_mrms = a3collect_prec.sum(axis=(1,2))

                a2num_mrms[iypart:eypart,ix] = (2*nrad+1)*(2*nrad+1)

                for ptype in lptype: 
                    a3collect_type_mask = ma.masked_not_equal(a3collect_type, ptype).mask

                    if a3collect_type_mask is np.True_:
                        d2num_type[ptype][iypart:eypart,ix] = 0
                        d2vfc_type[ptype][iypart:eypart,ix] = 0.

                    elif a3collect_type_mask is np.False_:
                        d2num_type[ptype][iypart:eypart,ix] = (2*nrad+1)*(2*nrad+1)
                        d2vfc_type[ptype][iypart:eypart,ix] = 1.0

                    else:
                        d2num_type[ptype][iypart:eypart,ix] = (~a3collect_type_mask).sum(axis=(1,2))

                        if ptype !=0:
                            d2vfc_type[ptype][iypart:eypart,ix] = ma.masked_invalid( ma.masked_where(a3collect_type_mask, a3collect_prec).sum(axis=(1,2)) / a1sum_mrms ).filled(-9999.)


                ##-- test ----
                #for ytmp in range(iypart,eypart+1):
                #    xtmp  = ix
                #    ptype = 1
                #    if d2num_type[ptype][ytmp,xtmp] >0:
                #        print 'num_type'
                #        print d2num_type[ptype][ytmp,xtmp]
                #        print 'collect_type'
                #        print a3collect_type[ytmp-iypart]
                #        print ''
                #        print 'vfc_type'
                #        print d2vfc_type[ptype][ytmp,xtmp]
                #        print 'mrms_prec'
                #        print a3collect_prec[ytmp-iypart]
                #        print ''
                #        print 'mrms'
                #        print a2ave_mrms[ytmp,ix]

                #        #sys.exit()
                #        #"print ''
                #        #"ymrms,xmrms= a1y_mrms[ytmp-iypart], a1x_mrms[ytmp-iypart]
                #        #"print a2mrms[ymrms-nrad:ymrms+nrad+1,xmrms-nrad:xmrms+nrad+1]
                #        #"print ''
                #        #"print ma.masked_where(a3collect_mask, a3collect_prec).mean(axis=(1,2)).filled(-9999.)[ytmp-iypart]
                #        #"print a3collect_prec[ytmp-iypart].shape
                #        #"print 'nrad=',nrad
                #        #"print 'pixres_nadir',pixres_nadir
                #        #"print 'pixres=',pixres
                #        #"print 'inc=',a2inc[0,ix]
                #        #"print 'pixara=', a1pixara[ix]




                ##-- test ----
                #ytmp, xtmp = 181, 5
                #if (ytmp >=iypart)&(ytmp < eypart):
                #    if xtmp == ix:
                #        print a2ave_mrms[ytmp,xtmp]
                #        print a3collect_prec[ytmp-iypart]
                #        print ''
                #        ymrms,xmrms= a1y_mrms[ytmp-iypart], a1x_mrms[ytmp-iypart]
                #        print a2mrms[ymrms-nrad:ymrms+nrad+1,xmrms-nrad:xmrms+nrad+1]
                #        print ''
                #        print ma.masked_where(a3collect_mask, a3collect_prec).mean(axis=(1,2)).filled(-9999.)[ytmp-iypart]
                #        print a3collect_prec[ytmp-iypart].shape
                #        print 'nrad=',nrad
                #        print 'pixres_nadir',pixres_nadir
                #        print 'pixres=',pixres
                #        print 'inc=',a2inc[0,ix]
                #        print 'pixara=', a1pixara[ix]

        a2goodfrac = ma.masked_where(a2ave_mrms<0, a2goodfrac).filled(-9999.)
        #*****************
        # Save
        #*****************
        outDir= tankDir + '/utsumi/PMM/MRMS/level2-pixel-match/%s.%s.%s/%04d/%02d/%02d'%(scantype, sensor, sate, Year, Mon, Day)
        util.mk_dir(outDir)
        ave_mrms_path= outDir + '/mrms.%06d.%05d-%05d.npy'%(oid,iscan,escan)
        goodfrac_path= outDir + '/goodfrac.%06d.%05d-%05d.npy'%(oid,iscan,escan)
        np.save(ave_mrms_path, a2ave_mrms.astype('float32'))
        np.save(goodfrac_path, a2goodfrac.astype('float32'))

        num_mrms_path = outDir + '/num-all.%06d.%05d-%05d.npy'%(oid,iscan,escan)
        np.save(num_mrms_path, a2num_mrms.astype('int16'))

        for ptype in lptype:
            ave_type_path= outDir + '/pr-%02d.%06d.%05d-%05d.npy'%(ptype, oid,iscan,escan)
            num_type_path= outDir + '/num-%02d.%06d.%05d-%05d.npy'%(ptype, oid,iscan,escan)

            np.save(num_type_path, d2num_type[ptype].astype('int16'))
            
            if ptype !=0:
                np.save(ave_type_path, d2vfc_type[ptype].astype('float32'))


        print(ave_mrms_path)
   
     


# %%
