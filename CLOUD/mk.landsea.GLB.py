from numpy import *
import socket, sys, os
import PIL.Image as Image
import myfunc.regrid.Regrid as Regrid
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util          as util

#------
Image.MAX_IMAGE_PIXELS = 933120000  # To avoid PIL.Image.DecompressionBombError
#------


if   socket.gethostname()=="well":
  baseDir = "/media/disk2/share"
  srcDir   = os.path.join(baseDir, "data/GLCC2/data")
elif socket.gethostname()=="mizu":
  baseDir = "/home/utsumi/mnt/wellshare"
  srcDir   = os.path.join(baseDir, "data/GLCC2/data")
elif socket.gethostname()=='shui':
  baseDir = '/tank/utsumi'
  srcDir   = os.path.join(baseDir, "GLCC2/data")
  print 'baseDir=',baseDir
  print 'srcDir=',srcDir
else:
  print 'hostname=',socket.gethostname()
  sys.exit()

#srcDir   = baseDir + "/data/GLCC2/data"
#srcDir   = "/home/utsumi/bin/PMM/CLOUD"
srcPath  = srcDir  + "/gblulcgeo20.tif"
#srcDir   = baseDir + "/data/GLCC2/seto"
#srcPath  = srcDir  + "/gusgs2_0ll.img"

NY     = 21600
NX     = 43200
nhead  = 173172 # Bytes
#nhead  = 0
res    = 30.0   # arc seconds

LatUp  = arange(-59.95,59.95+0.01,0.1)
LonUp  = arange(0.05, 359.95+0.01,0.1)

#LatUp  = arange(-89.95,89.95+0.01,0.1)
#LonUp  = arange(0.05, 359.95+0.01,0.1)

#LatUp  = arange(-37+0.25,37-0.25+0.01,0.5)
#LonUp  = arange(0.0+1.0, 360-1.0+0.01,2.0)

#BBox   = [[-89.95, 0.5],[89.95, 359.95]]  # degree
#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # degree
BBox    = [[-60, 0],[60, 360]]  # degree


llLAT_sec = -323985. -15.0   # arc seconds, lower boundary of the grid box
llLON_sec = -647985. -15.0   # arc seconds, left boundary of the grid box
urLAT_sec =  323985. +15.0   # arc seconds, lower boundary of the grid box
urLON_sec =  647985. +15.0   # arc seconds, left boundary of the grid box

lllat_sec = BBox[0][0] *3600.
lllon_sec = BBox[0][1] *3600.
urlat_sec = BBox[1][0] *3600.
urlon_sec = BBox[1][1] *3600.

iy   = int(floor((urLAT_sec - urlat_sec)/res))
ey   = int(floor((urLAT_sec - lllat_sec)/res))

ix   = int(floor((lllon_sec - llLON_sec)/res))
ex   = int(floor((urlon_sec - llLON_sec)/res))

IM = Image.open(srcPath)

if ((0<=ix)&(ex<NX)):
  im = IM.crop( (ix,iy,ex+1,ey+1) )
  a2fin = asarray(im)
elif ((ix<0)&(ex<NX)):
  im1 = IM.crop( (NX+ix,iy,NX,ey+1) )
  im2 = IM.crop( (0    ,iy,ex+1,ey+1) )
  a2fin = c_[asarray(im1), asarray(im2)]
elif ((0<=ix)&(NX<=ex)):
  im1 = IM.crop( (ix,iy,NX     ,ey+1) )
  im2 = IM.crop( (0 ,iy,ex-NX+1,ey+1) )
  a2fin = c_[asarray(im1), asarray(im2)]
else:
  print "check BBox"
  print BBox
  sys.exit()

a2fin = flipud(a2fin)

a2fin = ma.masked_not_equal(a2fin, 16).filled(1.0)  # land
a2fin = ma.masked_equal(a2fin, 16).filled(0.0)  # land
a2fin = a2fin.astype(float32)

# Upscale
cl     = CLOUDTYPE.CloudWNP()
us     = Regrid.UpScale()

LatOrg = arange(BBox[0][0], BBox[1][0]+1./120/2.0, 1./120.)
LonOrg = arange(BBox[0][1], BBox[1][1]+1./120/2.0, 1./120.)

us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)


a2up   = us.upscale(a2fin, pergrid=False, miss_in=-9999., miss_out=-9999.)
#outDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MASK"
#outDir = "/home/utsumi/mnt/wellshare/PMM/WNP.261x265/MASK"
#outDir = "/home/utsumi/mnt/wellshare/data/const"
outDir = '/tank/utsumi/validprof/const'
util.mk_dir(outDir)
outPath= outDir + "/landfrac.37SN.sa.%dx%d"%(len(LatUp),len(LonUp) )
a2up.tofile(outPath)
print outPath

