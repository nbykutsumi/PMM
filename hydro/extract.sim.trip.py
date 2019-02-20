from numpy import *
import numpy as np
import myfunc.util as util
from datetime import datetime, timedelta

iYear  = 1998
eYear  = 2009
lYear  = range(iYear,eYear+1)

#prcp   = 'JRA55'
lprcp   = ['GPCC']

#----------------------------------
def read_a3flw(prcp, Year):
    iDTime = datetime(Year,1,1)
    eDTime = datetime(Year,12,31)
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
    srcDir = '/data2/PMM/agcm5.6GW/out/PMM.Prcp_%s/river/0.4/riv_out_'%(prcp)

    nday   = (eDTime - iDTime).days + 1
    a3out  = zeros([nday,180,360],float32)
    for i,DTime in enumerate(lDTime):
        Mon,Day = DTime.timetuple()[1:3]
        srcPath = srcDir + '/PMM_____%04d%02d%02d.one'%(Year,Mon,Day)
        a2dat   = fromfile(srcPath,'float32').reshape(180,360).byteswap()
        a3out[i]= a2dat
    return a3out 
#----------------------------------
for prcp in lprcp:
    #-- Read GRDC location info for TRIP ---
    #loclistDir  = '/work/data1/hjkim/cf.mizu'
    loclistDir  = '/work/hk01/utsumi/PMM/hydro/list'
    loclistPath = loclistDir + '/siteinfo_1.0_utm_v02.csv'
    f=open(loclistPath,'r'); lines=f.readlines(); f.close()
    
    lstnid = []
    lxy    = []
    for line in lines[1:]:
        line  = line.split(',')
        stnid = int(line[0])
        lontrip=float(line[5])
        lattrip=float(line[6])
        if lontrip >=180:
            lontrip = lontrip - 360

        x     = int(lontrip - (-180)/1.0)
        y     = int((90-lattrip)/1.0)
        lstnid.append(stnid)
        lxy   .append([x,y])
    
    #-- Read simulation data --
    for Year in lYear:
        #flwPath    = simDir + '/outflw%04d.bin'%(Year)
        #a3flw      = fromfile(flwPath,float32).reshape(-1,720,1440)
        a3flw      = read_a3flw(prcp, Year)
    
        for istn, stnid in enumerate(lstnid):
            x,y = lxy[istn]
            a1flw = a3flw[:,y,x]
            print istn,len(a1flw)
    
            outDir = '/work/hk01/utsumi/PMM/hydro/outflw.trip/%s/%04d'%(prcp,Year)
            util.mk_dir(outDir)
            outPath= outDir + '/flwout.%07d.npy'%(stnid)
            np.save(outPath, a1flw)
            print outPath
    
        

