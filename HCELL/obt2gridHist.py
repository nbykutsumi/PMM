import os, sys, socket
import myfunc.util      as util
import calendar
from   datetime         import datetime, timedelta
from   numpy            import *
#from   myfunc.IO        import GPyM
from   myfunc.IO        import GPM
from   hcell_fsub import *

hostname = socket.gethostname()
if  hostname =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    orootDir   = "/media/disk2/share"

orootDir = os.path.join(orootDir, "HCELL/PDF/ALL")


iYM     = [1999,1]
eYM     = [2003,12]
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

for Year,Mon in lYM:

    aCount= zeros([len(Bin)+2,ny,nx],int32)
    lPath = gpm.list_granule(Year,Mon)
    for iPath,srcPath in enumerate(lPath):
    #for iPath,srcPath in enumerate(lPath[153:156]):
        print iPath,srcPath

        stormH  = gpm.load_var_granule(srcPath,Var)
        Lat     = gpm.load_var_granule(srcPath,"Latitude")
        Lon     = gpm.load_var_granule(srcPath,"Longitude")
        #rainType= gpm.load_var_granule(srcPath,"rainType").flatten()

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

        #--------------------------
        # No rain or zero height
        aMsk     = ma.masked_greater(stormH,0).mask
        stormHtmp   = ma.masked_where(aMsk, stormH).compressed()
        Lattmp   = ma.masked_where(aMsk, Lat   ).compressed()
        Lontmp   = ma.masked_where(aMsk, Lon   ).compressed()

        if len(stormHtmp)==0:
            aTmp = zeros([nx,ny],int32)
        else:
            aTmp     = hcell_fsub.obt2grid_count(Lontmp, Lattmp, lon_first, lat_first, dlon, dlat, nx, ny)

        aCount[0] = aCount[0] + aTmp.T


#        print "zero=",len(stormHtmp)
#        print "aCount[0]",aCount[0].sum()
        #--------------------------
        # Rain events
        aMsk     = ma.masked_less_equal(stormH,0).mask
        stormHtmp= ma.masked_where(aMsk, stormH).compressed()
        Lattmp   = ma.masked_where(aMsk, Lat   ).compressed()
        Lontmp   = ma.masked_where(aMsk, Lon   ).compressed()

        if len(stormHtmp)==0:
            aTmp = zeros([nx,ny,len(Bin)+1],int32)
        else:
            aTmp  = hcell_fsub.obt2grid_hist(stormHtmp, Lontmp, Lattmp, Bin, lon_first, lat_first, dlon, dlat, nx, ny)
      
        for iz in range(len(Bin)+1):
            aCount[iz+1] = aCount[iz+1] + aTmp[:,:,iz].T

#        print "non-zero=",len(stormHtmp)
#        print "aCount[1:]",aCount[1:].sum()
#
#        tmp  = ma.masked_greater_equal(stormHtmp,12000)
#        print "0< <12000",len(tmp.compressed())
#        print "aCount[1:13]",aCount[1:13].sum()


    #------------------------------
    # Write
    #------------------------------
    outDir = os.path.join(orootDir, "%04d"%(Year))
    outPath= os.path.join(outDir, "num.%02d.%02dx%03dx%03d"%(Mon,len(Bin)+2, ny, nx))
    util.mk_dir(outDir)
    aCount.tofile(outPath)
    print outPath


