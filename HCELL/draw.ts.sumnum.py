from numpy import *
import myfunc.util as util
import os, sys, socket
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')


iYM = [1997,12]
eYM = [2013,11]
#eYM = [2006,11]
dMon= 1
lYM = util.ret_lYM(iYM, eYM)[::dMon]

hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

rootDir = os.path.join(rootDir, "HCELL/SUMNUM/Type")

#lrtype = ["strat","conv","other","all"]
lrtype = ["all"
ny  = 148
nx  = 180
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)
llndsea = ["sea","lnd"]
#-------------------------------
# Function
#---------
def load_landfrac():
    srcPath = os.path.join(rootDir,"data/const/landfrac.37SN.148x180")
    return fromfile(srcPath,float32).reshape(148,180)


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


def DrawTimeLat(a2in, figPath, stitle, cmap="YlGnBu", vmin=None, vmax=None,season=False):
    if season ==False:
        fig = plt.figure(figsize=(6,3))
    else:
        fig = plt.figure(figsize=(4,2.5))
    ax  = fig.add_axes([0.15,0.15,0.80,0.64])

    ny,nx = a2in.shape
    x   = range(nx)
    y   = Lat
    X,Y = meshgrid(x,y)

    im  = ax.pcolormesh(X,Y,a2in,cmap=cmap,vmin=vmin,vmax=vmax)

    # Colorbar
    cb  = plt.colorbar(im, orientation="vertical")
    cb.ax.tick_params(labelsize=8)

    # X-tick label
    if season ==False:
        plt.xticks(x[1::12])
        lxlabel = range(iYM[0]+1,eYM[0]+1)
    else:
        plt.xticks(x[::2])
        lxlabel = range(iYM[0]+1,eYM[0]+1,2)

    ax.xaxis.set_ticklabels(lxlabel, fontsize=8)

    # Y-tick label
    yticks  = arange(-35,35+0.1,5)
    plt.yticks(yticks)
    ax.yaxis.set_ticklabels(yticks, fontsize=8)

    # Title
    plt.title(stitle,fontsize=9)

    # Save
    plt.savefig(figPath)
    plt.close()
    print figPath



def PlotTime(a2in, lname, figPath, legPath, stitle):
    fig = plt.figure(figsize=(3,3))
    ax  = fig.add_axes([0.15,0.15,0.80,0.64])

    nline,nx = a2in.shape
    x   = range(nx)

    [ax.plot(x,a2in[iline]) for iline in range(nline)]

    plt.xticks(x[::2])
    lxlabel = range(iYM[0]+1,eYM[0]+1,2)
    ax.xaxis.set_ticklabels(lxlabel, fontsize=8)

    # Title
    plt.title(stitle,fontsize=9)

    # Save
    plt.savefig(figPath)
    print figPath

    # Legend file
    if legPath == None:
        return

    lines   = ax.get_lines()
    figleg  = plt.figure(figsize=(2,2))
    figleg.legend(lines, lname)
    figleg.savefig(legPath)   
    plt.close()

def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))

#-------------------------------

d2dat = {rtype,lndsea:zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}
d2num = {rtype,lndsea:zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}
d2sum = {rtype,lndsea:zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}

for it,(iYear,iMon) in enumerate(lYM):
    eYear,eMon = util.shift_YM(iYear,iMon,dMon-1)
    print iYear,iMon
    print eYear,eMon
    for rtype in lrtype:
        for lndsea in llndsea:
            a2num = ret_sumnum("num",[iYear,iMon],[eYear,eMon], rtype)
            a2sum = ret_sumnum("sum",[iYear,iMon],[eYear,eMon], rtype)
            a1num= a2num.mean(axis=1) 
            a1sum= a2sum.mean(axis=1) 
               
            a1dat= ma.masked_invalid(a1sum.astype(float32)/a1num).filled(0) 
    
            d2sum[rtype,lndsea][:,it] = a1sum
            d2num[rtype,lndsea][:,it] = a1num
            d2dat[rtype,lndsea][:,it] = a1dat


#-- monthly ----
for rtype in lrtype:
    figDir  = "/tank/utsumi/PMM/HCELL/pict"
#    for lndsea in llndsea:
#        # Storm Height
#        a2out   = d2dat[rtype]
#        figPath = os.path.join(figDir, "TimeLat.sumnum.stormH.%s.png"%(rtype))
#        stitle  = "stormH (%s) 2A23"%(rtype)
#        DrawTimeLat(a2out, figPath, stitle) 
#   
#        # Normalized
#        a2clim  = mkClimMon(d2dat[rtype])
#        a2out   = (d2dat[rtype]-a2clim) / a2clim
#   
#        figPath = os.path.join(figDir, "TimeLat.sumnum.stormH.Norm.%s.png"%(rtype))
#        stitle  = "Normalized stormH (%s) 2A23"%(rtype)
#        cmap    = "RdBu_r"
#        DrawTimeLat(a2out, figPath, stitle, cmap=cmap, vmin=-0.1, vmax=0.1)


#-- Seasonal ----
lseason = ["DJF","MAM","JJA","SON"]
dik     = {"DJF":0,"MAM":3,"JJA":6,"SON":9}
if iYM[1] != 12:
    print "iYM=",iYM
    print "the first month must be Dec."
    sys.exit()
if eYM[1] != 11:
    print "eYM=",eYM
    print "the last month must be Nov."
    sys.exit()


for rtype in ["all"]:
    for lndsea in lndsea:
        figDir  = "/tank/utsumi/PMM/HCELL/pict"
    
        for season in lseason:
            ik   = dik[season] 
    
            # Storm Height
            a2sum =   d2sum[rtype,lndsea][:,ik::12]\
                     +d2sum[rtype,lndsea][:,ik+1::12]\
                     +d2sum[rtype,lndsea][:,ik+2::12]
    
            a2num =   d2num[rtype,lndsea][:,ik::12]\
                     +d2num[rtype,lndsea][:,ik+1::12]\
                     +d2num[rtype,lndsea][:,ik+2::12]
    
            a2dat = ma.masked_invalid(a2sum/a2num)
    
            a2out = a2dat
            figPath = os.path.join(figDir, "TimeLat.sumnum.stormH.%s.%s.%s.png"%(season,rtype),lndsea)
            stitle  = "stormH (%s) @%s %s 2A23"%(rtype,season,lndsea)
            DrawTimeLat(a2out, figPath, stitle, vmin=1000., season=True)
    
    
            # Normalized 
            a2sum =   d2sum[rtype,lndsea][:,ik::12]\
                     +d2sum[rtype,lndsea][:,ik+1::12]\
                     +d2sum[rtype,lndsea][:,ik+2::12]
    
            a2num =   d2num[rtype,lndsea][:,ik::12]\
                     +d2num[rtype,lndsea][:,ik+1::12]\
                     +d2num[rtype,lndsea][:,ik+2::12]
    
            a2dat = ma.masked_invalid(a2sum/a2num)
            a2clim= mkClim(a2dat)
    
            a2out = ma.masked_invalid((a2dat - a2clim)/a2clim)
    
            figPath = os.path.join(figDir, "TimeLat.sumnum.stormH.Norm.%s.%s.%s.png"%(season,rtype,lndsea))
            stitle  = "stormH (%s) @%s %s 2A23"%(rtype,season,lndsea)
            cmap    = "RdBu_r"
            DrawTimeLat(a2out, figPath, stitle, cmap=cmap, vmin=-0.1, vmax=0.1, season=True)
    #
    
    #    # Line plot
    #    for season in lseason:
    #        ik    = dik[season] 
    #        a2sum =   d2sum[rtype,lndsea][:,ik::12]\
    #                 +d2sum[rtype,lndsea][:,ik+1::12]\
    #                 +d2sum[rtype,lndsea][:,ik+2::12]
    #
    #        a2num =   d2num[rtype,lndsea][:,ik::12]\
    #                 +d2num[rtype,lndsea][:,ik+1::12]\
    #                 +d2num[rtype,lndsea][:,ik+2::12]
    #
    #        a2dat = ma.masked_invalid(a2sum/a2num)
    #
    #        lielat= [[5,10],[10,15],[15,20],[25,30],[30,35]] 
    #        nx    = a2dat.shape[1]
    #        a2out = empty([len(lielat),nx])
    #        for iline,(ilat,elat) in enumerate(lielat):
    #            iy    = ret_y(ilat)
    #            ey    = ret_y(elat)
    #            a1out = a2dat[iy:ey].mean(axis=0)
    #            a2out[iline] = a1out
    #
    #        figPath = os.path.join(figDir, "plot.sumnum.stormH.%s.%s.png"%(season,rtype))
    #        stitle  = "stormH (%s) @%s %s 2A23"%(rtype, season, lndsea)
    #
    #        legPath = os.path.join(figDir, "leg.plot.stormH.png")
    #        lname   = ["%.1f-%.1f"%(ilat,elat) for (ilat,elat) in lielat]
    #        PlotTime(a2out, lname, figPath, legPath, stitle)
    #

