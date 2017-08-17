import os, sys, socket
import myfunc.util      as util
import calendar
from   datetime         import datetime, timedelta
from   numpy            import *
from   myfunc.IO        import GPM
from   hcell_fsub import *

hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    orootDir   = "/media/disk2/share"

orootDir = os.path.join(orootDir, "HCELL/PDF/Type")


iYM     = [1997,12]
eYM     = [2014,1]
#eYM     = [1998,1]
lYM     = util.ret_lYM(iYM,eYM)

#prj     = "TRMM.PR"
#prdLv   = "L2A23"
#prdVer  = "07"
#gpm     = GPyM.GPM(prj, prdLv, prdVer)

version = "07"
gpm     = GPM.L2A23(version=version)
Var     = "stormH"
lat_first = -37.0
lon_first = 0.0
dlat      = 0.5  # deg
dlon      = 2.0  # deg
ny        = 148
nx        = 180

Bin       = arange(1000,12000+1,1000)
lrtype    = ["strat","conv","other"]
dlower    = {"strat":100, "conv":200, "other":300}


for Year,Mon in lYM:

    aNoRain= zeros([ny,nx], int32)
    dCount = {rtype:zeros([len(Bin)+1,ny,nx], int32) for rtype in lrtype}

    lPath = gpm.list_granule(Year,Mon)
    for srcPath in lPath:
    #for srcPath in lPath[0:1]:
        print srcPath

        stormH  = gpm.load_var_granule(srcPath,Var)
        Lat     = gpm.load_var_granule(srcPath,"Latitude")
        Lon     = gpm.load_var_granule(srcPath,"Longitude")
        rType   = gpm.load_var_granule(srcPath,"rainType")

        #--------------------------
        # Mask Missing value
        # -8888 No rain
        # -1111 Rain is not present with high level conficence
        # -9999 Missing value
        #--------------------------
        Htmp = stormH


        aMsk    = ma.masked_equal(stormH,-9999).mask
        stormH  = ma.masked_where(aMsk,stormH).compressed()
        Lat     = ma.masked_where(aMsk, Lat  ).compressed()
        Lon     = ma.masked_where(aMsk, Lon  ).compressed()
        rType   = ma.masked_where(aMsk, rType).compressed()

        #--------------------------
        # No rain or zero height
        aMsk     = ma.masked_greater(stormH,0).mask
        stormHtmp   = ma.masked_where(aMsk, stormH).compressed()
        Lattmp   = ma.masked_where(aMsk, Lat   ).compressed()
        Lontmp   = ma.masked_where(aMsk, Lon   ).compressed()
        rTypetmp = ma.masked_where(aMsk, rType ).compressed()

        if ((len(stormHtmp)==0)or(len(rTypetmp)==0)):
            aTmp = zeros([nx,ny],int32)
        else:
            aTmp     = hcell_fsub.obt2grid_count(Lontmp, Lattmp, lon_first, lat_first, dlon, dlat, nx, ny)

        aNoRain = aNoRain + aTmp.T

        #print "*"*50
        #print "zero=",len(stormHtmp)
        #print "aNoRain",aNoRain.sum()


        #--------------------------
        # Rain events
        aMsk   = ma.masked_less_equal(stormH,0).mask
        stormH = ma.masked_where(aMsk, stormH).compressed()
        Lat    = ma.masked_where(aMsk, Lat   ).compressed()
        Lon    = ma.masked_where(aMsk, Lon   ).compressed()
        rType  = ma.masked_where(aMsk, rType ).compressed()

        #--------------------------
        # Rain-type
        
        tmpsum0 = 0
        tmpsum1 = 0
 
        for rtype in lrtype:
            lower = dlower[rtype]

            # screen by rain type
            aMsk     = ma.masked_outside(rType, lower, lower+99).mask
            stormHtmp= ma.masked_where(aMsk, stormH).compressed()
            Lattmp   = ma.masked_where(aMsk, Lat).compressed()
            Lontmp   = ma.masked_where(aMsk, Lon).compressed()
            rTypetmp = ma.masked_where(aMsk, rType).compressed()


            if ((len(stormHtmp)==0)or(len(rTypetmp)==0)):
                aTmp = zeros([nx,ny,len(Bin)+1],int32)
            else:    
                aTmp = hcell_fsub.obt2grid_hist(stormHtmp, Lontmp, Lattmp, Bin, lon_first, lat_first, dlon, dlat, nx, ny)
          
            for iz in range(len(Bin)+1):
                dCount[rtype][iz] = dCount[rtype][iz] + aTmp[:,:,iz].T
    
            #print " "*50
            #print "non-zero",rtype,len(stormHtmp)
            #print "dCount",dCount[rtype].sum()
            tmpsum0 = tmpsum0 + len(stormHtmp)
            tmpsum1 = tmpsum1 + dCount[rtype].sum()
 
    #print " "*50
    #print tmpsum0, tmpsum1
    #------------------------------
    # Write
    #------------------------------
    outDir = os.path.join(orootDir, "%04d"%(Year))
    util.mk_dir(outDir)

    # Write : No rain
    outPath= os.path.join(outDir, "num.non.%02d.%03dx%03d"%(Mon, ny, nx))
    aNoRain.tofile(outPath)

    # Write : with rain
    for rtype in lrtype:        
        outPath= os.path.join(outDir, "num.%s.%02d.%02dx%03dx%03d"%(rtype,Mon, len(Bin)+1,ny, nx))
        dCount[rtype].tofile(outPath)
        print outPath
    

