from numpy import *
import socket, sys
import Image
import myfunc.regrid.Regrid as Regrid
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

if   socket.gethostname()=="well":
  baseDir = "/media/disk2/share"
elif socket.gethostname()=="mizu":
  baseDir = "/home/utsumi/mnt/wellshare"

#srcDir   = baseDir + "/data/GLCC2/data"
srcDir   = "/home/utsumi/bin/PMM/CLOUD"
srcPath  = srcDir  + "/gblulcgeo20.tif"
#srcDir   = baseDir + "/data/GLCC2/seto"
#srcPath  = srcDir  + "/gusgs2_0ll.img"

NY     = 21600
NX     = 43200
nhead  = 173172 # Bytes
#nhead  = 0
res    = 30.0   # arc seconds

llLAT_sec = -323985. -15.0   # arc seconds, lower boundary of the grid box
llLON_sec = -647985. -15.0   # arc seconds, left boundary of the grid box
urLAT_sec =  323985. +15.0   # arc seconds, lower boundary of the grid box
urLON_sec =  647985. +15.0   # arc seconds, left boundary of the grid box

BBox   = [[-89.95, 0.5],[89.95, 359.95]]  # degree
#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # degree
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

#LatUp  = cl.Lat
#LonUp  = cl.Lon

LatUp  = arange(-89.95,89.95+0.01,0.1)
LonUp  = arange(0.05, 359.95+0.01,0.1)

us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)


a2up   = us.upscale(a2fin, pergrid=False, miss_in=-9999., miss_out=-9999.)
#outDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MASK"
outDir = "/home/utsumi/mnt/wellshare/PMM/WNP.261x265/MASK"
outPath= outDir + "/landfrac.%dx%d"%(len(LatUp),len(LonUp) )
a2up.tofile(outPath)
print outPath

