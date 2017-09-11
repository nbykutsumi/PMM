from   numpy import *
from   myfunc.fig import Fig
import myfunc.util as util
import os, sys, socket
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')

#lYear    = [1998]
lYear    = [1998,2005,2012]

hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

rootDir = os.path.join(rootDir, "HCELL/SUMNUM/Type")

ny  = 148
nx  = 180
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)
BBox= [[-37,0],[37,360]]
#-------------------------------
# Function
#---------
def ret_sumnum_single(sumnum,Year,Mon,rtype):
    if rtype !="all":
        srcDir = os.path.join(rootDir, "%04d"%(Year))
        srcPath= os.path.join(srcDir, "%s.%s.%02d.%03dx%03d"%(sumnum, rtype,Mon, ny, nx))
        a2in   = fromfile(srcPath, int32).reshape(ny,nx)
    else:
        a2strat = ret_sumnum_single(sumnum,Year,Mon,"strat") 
        a2conv  = ret_sumnum_single(sumnum,Year,Mon,"conv") 
        a2other = ret_sumnum_single(sumnum,Year,Mon,"other")
        return a2strat+a2conv+a2other 

    return a2in

def ret_sumnum(sumnum,iYM,eYM,rtype):
    lYM  = util.ret_lYM(iYM,eYM)
    a2dat= zeros([ny,nx],int32)
    for Year,Mon in lYM:
       a2dat = a2dat + ret_sumnum_single(sumnum,Year,Mon,rtype)
    return a2dat

def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))


def mkClimMon(a2in):
    # Seasonal climatology
    a2clim    = empty(a2in.shape)
    for iM in range(12):
        a1climMon= a2in[:,iM::12].mean(axis=1).reshape(ny,-1)
        a2clim[:,iM::12] = a1climMon
    return a2clim

def mkClim(a2in):
    # Simple climatology
    ny,nx     = a2in.shape
    a2clim    = empty(a2in.shape)
    a1clim    = a2in.mean(axis=1)
    for ix in range(nx):
        a2clim[:,ix] = a1clim
    return a2clim


def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))


#------------------------------

figsize = (7,2)
for Year in lYear:
    iYM = [Year,1]
    eYM = [Year,12]
    #------------------------------
    # Rain event
    a2num = ret_sumnum("num",iYM, eYM, "all")
    
    cmap    = "jet"
    stitle  = "counts rain %d"%(Year)
    figDir  = "/tank/utsumi/PMM/HCELL/pict"
    figname = os.path.join(figDir, "map.count.rain.%d.png"%(Year))
    cbarname= os.path.join(figDir, "cbar.count.png")
    
    bnd     = range(0,4000,200)
    Fig.DrawMapSimple(a2in=a2num, a1lat=Lat, a1lon=Lon, BBox=BBox, figname=figname, bnd=bnd, cbarname=cbarname, stitle=stitle, cmap=cmap, figsize=figsize)
    
    
    # all observations
    a2rain = ret_sumnum("num",iYM, eYM, "all")
    a2no   = ret_sumnum("num",iYM, eYM, "non")
    a2num  = a2rain + a2no
    
    cmap    = "jet"
    stitle  = "counts obs %d"%(Year)
    figDir  = "/tank/utsumi/PMM/HCELL/pict"
    figname = os.path.join(figDir, "map.count.obs.%d.png"%(Year))
    cbarname= os.path.join(figDir, "cbar.count.obs.png")
    
    bnd     = range(0,100000,10000)
    Fig.DrawMapSimple(a2in=a2num, a1lat=Lat, a1lon=Lon, BBox=BBox, figname=figname, bnd=bnd, cbarname=cbarname, stitle=stitle, cmap=cmap, figsize=figsize)
