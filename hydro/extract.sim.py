from numpy import *
import numpy as np
import myfunc.util as util


iYear  = 1998
eYear  = 1999
lYear  = range(iYear,eYear+1)

#prcp   = 'JRA55'
lprcp   = ['GPCC']

for prcp in lprcp:
    #simDir = '/work/hk01/hjkim/CaMa-Flood/CaMa-Flood_v3.6.2_20140909/out/global_15min.soc18'
    simDir = '/data2/PMM/out/PMM.Prcp_%s/river'%(prcp)
    #loclistDir  = '/work/hk01/hjkim/CaMa-Flood/CaMa-Flood_v3.6.2_20140909/map/global_15min'
    #loclistPath = loclistDir + '/grdc_loc.txt'
    loclistPath = '/work/hk01/utsumi/PMM/hydro/list/grdc_loc_rev_1998_2009.txt'
    
    #upareaPath = loclistDir + '/uparea.bin'
    #a2uparea   = fromfile(upareaPath,float32).reshape(720,1440)
    
    #-- Read GRDC location info ---
    f=open(loclistPath,'r'); lines=f.readlines(); f.close()
    
    lstnid = []
    lxy    = []
    luparea= []
    for line in lines:
        line  = line.split()
        stnid = int(line[0])
        x     = int(line[6])  # 0,1,2 ..
        y     = int(line[7])  # 0,1,2 ..
        uparea= float(line[9])
        lstnid.append(stnid)
        lxy   .append([x,y])
        luparea.append(uparea)
    
    #-- Read simulation data --
    for Year in lYear:
        flwPath    = simDir + '/outflw%04d.bin'%(Year)
        a3flw      = fromfile(flwPath,float32).reshape(-1,720,1440)
    
    
        for istn, stnid in enumerate(lstnid):
            x,y = lxy[istn]
            a1flw = a3flw[:,y,x]
            print istn,len(a1flw)
    
            outDir = '/work/hk01/utsumi/PMM/hydro/outflw/%s/%04d'%(prcp,Year)
            util.mk_dir(outDir)
            outPath= outDir + '/flwout.%07d.npy'%(stnid)
            np.save(outPath, a1flw)
            print outPath
    
        

