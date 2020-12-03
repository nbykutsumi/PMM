# %%
import matplotlib
matplotlib.use("Agg")
% matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
from datetime import datetime, timedelta

calcflag  = True

amsr2     = ["GCOMW1","AMSR2"]
ssmis_f16 = ["F16","SSMIS"]
ssmis_f17 = ["F17","SSMIS"]   # 37V is missing since Apr 2016.
ssmis_f18 = ["F18","SSMIS"]
atms_npp  = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]   # MRMS does not have NOAA20 for any year.

mhs_metopa= ["METOPA","MHS"]
mhs_metopb= ["METOPB","MHS"]
mhs_noaa18= ["NOAA18","MHS"]  # Not available at arthurhou.pps after 2018/10/21
mhs_noaa19= ["NOAA19","MHS"]

#lspec = [amsr2, ssmis_f16, ssmis_f18, atms_npp, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
mainscan_amsr2 = 1   # for AMSR2, SSMIS
mainscan_ssmis = 1   # for AMSR2, SSMIS
lspec = [amsr2]
lrettype= ['epc']
#lrettype= ['gpr']

iym = [2018,1]
eym = [2018,1]

#prmin = 0.1
#prmin = 0.01
prmin = 0.0
#** Constants ******

lsurftype = ['all','ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'all': 'All surface' 
            ,'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

#**************************************
# Functions
#--------------------------------------
def read_orbitlist(Year,Mon):
    listPath = listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath); lines=f.readlines(); f.close()
    lout=[]
    for line in lines:
        line = map(int,line.strip().split(','))
        lout.append(line)
    return lout

def pickyx_low2high_amsr2(nyhi):
    nxhi  = 486
    axpick = ((np.arange(nxhi) -1)/2).astype('int16')
    axpick[0] = 0    # [0,0,0,1,1,2,2,3,3,...]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_low2high_ssmis(nyhi):
    nxhi = 180
    axpick = (np.arange(nxhi) /2).astype('int16')  # [0,0,1,1,2,2,3,3,...,89,89]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_amsr2(nyhi):
    nxlw  = 243
    axpick = (np.arange(nxlw)*2+2).astype('int16')
    axpick[-1] = nxlw*2-1   # [2,4,6,...482,484,485]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_ssmis(nyhi):
    nxlw = 90
    axpick = (np.arange(nxlw)*2+1).astype('int16') # [1,3,5,...,177,179]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick


#--------------------------------------

dvnummax = {}
for (rettype,sate,sensor) in [(rettype,sate,sensor)
                for rettype  in lrettype
                for (sate,sensor) in lspec]:

    print rettype,sate, sensor

    myhost = socket.gethostname()
    if myhost == 'shui':
        tankdir = '/tank'
        workdir = '/work'
    elif myhost == 'well':
        tankdir = '/home/utsumi/mnt/lab_tank'
        workdir = '/home/utsumi/mnt/lab_work'

    if rettype=='epc':
        retbasedir = tankdir + '/utsumi/PMM/retepc/us.v01'

    elif rettype=='gpr':
        sys.exit()
    else: 
        print 'check rettype', sys.exit()

    #------------------------------
    # make orbit list
    #------------------------------
    lYM   = util.ret_lYM(iym,eym)
    lorbit = []
    for (Year,Mon) in lYM:
        listPath=tankdir + '/utsumi/PMM/US/obtlist/overpass.%s.%s.%04d.%02d.csv'%(sensor,sate,Year,Mon)
        f=open(listPath,'r'); lines=f.readlines(); f.close()
        for line in lines:
            _,_,Day,oid,iscan,escan = map(int,line.strip().split(','))

            #if (idtime<=datetime(Year,Mon,Day))&(datetime(Year,Mon,Day)<=edtime):
            if (datetime(Year,Mon,Day) < datetime(2018,1,23)): continue
            #if (datetime(Year,Mon,Day) < datetime(2018,1,31)): continue
            lorbit.append(map(int,line.strip().split(',')))


        a1epc = []
        a1gpr = []
        a1mrms= []
        a1lat = []
        a1lon = []
        a1stype = []
        a1good  = []
        a1gprqc = []
        for year,mon,day,oid,iscan,escan in lorbit:
            print year,mon,day,oid

            #----------------------------
            # epc
            #----------------------------
            if rettype =='epc':
                #retdir = retbasedir + '/us.v01/ATMS.NPP/2018/01/25'
                retdir = retbasedir + '/%s.%s/%04d/%02d/%02d'%(sensor,sate,year,mon,day)
                prpath = retdir + '/nsurfNScmb.%06d.y%04d-%04d.nrec10000.npy'%(oid,iscan,escan)
                latpath= retdir + '/lat.%06d.y%04d-%04d.nrec10000.npy'%(oid,iscan,escan)
                lonpath= retdir + '/lon.%06d.y%04d-%04d.nrec10000.npy'%(oid,iscan,escan)
                if not os.path.exists(prpath):
                    print 'no file', prpath
                    continue
                

                a2epc = np.load(prpath)
                a2lat = np.load(latpath)
                a2lon = np.load(lonpath)

            #----------------------------
            # gprof
            #----------------------------
            #gprbasedir = workdir + '/hk02/PMM/NASA/GPM.GMI/2A/V05'
            gprbasedir = '/media/disk2/data/PMM/NASA/%s.%s/2A-CLIM/V05'%(sate,sensor)
            gprdir = gprbasedir + '/%04d/%02d/%02d'%(year,Mon,Day)
            ssearch  = gprdir + '/2A*.HDF5'
            lgprpath = sort(glob.glob(ssearch))

            if len(lgprpath)==0:
                print 'No GPROF file'
                print ssearch
                sys.exit()

            gprpath = lgprpath[0]
            with h5py.File(gprpath, 'r') as h:
                a2gpr   = h['S1/surfacePrecipitation'][:]
                a2gprqc = h['S1/qualityFlag'][:]
                a2stype = h['S1/surfaceTypeIndex'][:]
                a2glat  = h['S1/Latitude'][:]
                a2glon  = h['S1/Longitude'][:]



            nyobt = a2glat.shape[0]
            if sensor=='AMSR2':
                if mainscan_amsr2==1:
                    a2ypick, a2xpick = pickyx_high2low_amsr2(nyobt)

                elif mainscan_amsr2==5:
                    a2ypick, a2xpick = pickyx_low2high_amsr2(nyobt)

                else:
                    print 'check mainscan', sensor, mainscan

            elif sensor=='SSMIS':
                if mainscan_ssmis==1:
                    a2ypick, a2xpick = pickyx_high2low_ssmis(nyobt)

                elif mainscan_ssmis==4:
                    a2ypick, a2xpick = pickyx_low2high_ssmis(nyobt)

                else:
                    print 'check mainscan', sensor, mainscan

            if sensor in ['AMSR2','SSMIS']:
                a2gpr   = a2gpr  [a2ypick,a2xpick]
                a2gprqc = a2gprqc[a2ypick,a2xpick]
                a2stype = a2stype[a2ypick,a2xpick]
                a2glat  = a2glat [a2ypick,a2xpick]
                a2glon  = a2glon [a2ypick,a2xpick]

            #--- extract iscan to escan --------
            a2gpr   = a2gpr  [iscan:escan+1]
            a2gprqc = a2gprqc[iscan:escan+1]
            a2stype = a2stype[iscan:escan+1]
            a2glat  = a2glat [iscan:escan+1]
            a2glon  = a2glon [iscan:escan+1]

            #-- test: check distance ---
            #RADEARTH = 6371
            #DTR      = 0.017453
            #a2dist = RADEARTH*np.arccos(np.cos(DTR*a2lon-DTR*a2glon)*np.cos(DTR*a2lat)*np.cos(DTR*a2glat) + np.sin(DTR*a2lat)*np.sin(DTR*a2glat))

            #----------------------------
            # mrms
            #----------------------------
            mrmsdir ='/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match/%s.%s.%s/%04d/%02d/%02d'%(rettype,sensor,sate,year,mon,day)

            mrmspath = mrmsdir + '/mrms.%06d.%05d-%05d.npy'%(oid,iscan,escan)
            gfracpath= mrmsdir + '/goodfrac.%06d.%05d-%05d.npy'%(oid,iscan,escan)
            if not os.path.exists(mrmspath):
                print 'no file', prpath
                continue

            a2mrms = np.load(mrmspath) 
            a2good = np.load(gfracpath)

            a2maskgood = ma.masked_less(a2good, 0.8).mask

            a1epc  .extend( ma.masked_where(a2maskgood, a2epc).compressed().astype('float32'))
            a1gpr  .extend( ma.masked_where(a2maskgood, a2gpr).compressed().astype('float32'))
            a1mrms .extend( ma.masked_where(a2maskgood, a2mrms).compressed().astype('float32'))
            a1lat  .extend( ma.masked_where(a2maskgood, a2lat).compressed().astype('float32'))
            a1lon  .extend( ma.masked_where(a2maskgood, a2lon).compressed().astype('float32'))
            a1stype.extend( ma.masked_where(a2maskgood, a2stype).compressed().astype('int32'))
            a1gprqc.extend( ma.masked_where(a2maskgood, a2gprqc).compressed().astype('int8'))
            a1good.extend( ma.masked_where(a2maskgood, a2good).compressed().astype('float32'))

        #--- Save ----------
        outdir = tankdir + '/utsumi/multi/us.pair/%04d/%02d'%(year,mon)
        util.mk_dir(outdir)

        np.save( outdir + '/nsurf-epc.npy', np.array(a1epc))
        np.save( outdir + '/nsurf-gpr.npy', np.array(a1gpr))
        np.save( outdir + '/gprqc.npy', np.array(a1gprqc))
        np.save( outdir + '/nsurf-mrms.npy', np.array(a1mrms))
        np.save( outdir + '/good-mrms.npy', np.array(a1good))
        np.save( outdir + '/surftype.npy', np.array(a1stype))
        np.save( outdir + '/lat.npy', np.array(a1lat))
        np.save( outdir + '/lon.npy', np.array(a1lon))

        print outdir

        a1epc = ma.masked_less(np.array(a1epc),0).filled(0)
        a1gpr  = ma.masked_less(np.array(a1gpr),0).filled(0)
        a1mrms = ma.masked_less(np.array(a1mrms),0).filled(0)

        plt.figure()
        plt.hist(np.array(a1gprqc))
        #plt.figure()
        #plt.plot(a1mrms, a1epc, 'o')
        #plt.show()

        #plt.figure()
        #plt.plot(a1mrms, a1gpr, 'o')
        #plt.show()
# %%
