from numpy import *
import myfunc.util as util
import os, sys, socket
import matplotlib.pyplot as plt
import matplotlib
from   matplotlib.ticker import AutoMinorLocator
from   scipy import stats

#matplotlib.use('Agg')


iYear = 1998
eYear = 2013
hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

baseDir = os.path.join(rootDir, "HCELL/SUMNUM/Type")

#lrtype = ["strat","conv","other","all"]
lrtype = ["all"]
ny  = 148
nx  = 180
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)
llndsea = ["sea","lnd"]
#llndsea = ["sea"]
#-------------------------------
# Function
#---------
def load_landfrac():
    srcPath = os.path.join(rootDir,"data/const/landfrac.37SN.148x180")
    return fromfile(srcPath,float32).reshape(148,180)


def ret_sumnum_single(sumnum,Year,Mon,rtype):
    if rtype !="all":
        srcDir = os.path.join(baseDir, "%04d"%(Year))
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


def DrawLatH(dH,stitle, figPath, legPath):

    lYear   = dH.keys()
    figplot = plt.figure(figsize=(5,3))
    axplot  = figplot.add_axes([0.15,0.15,0.8,0.7])

    # color
    cmap  = plt.cm.jet
    nYear = len(lYear)
    lcmap = [cmap(int(i)) for i in linspace(0,cmap.N,nYear)]

    x     = Lat
    for i,Year in enumerate(lYear):
        axplot.plot(x[1:-1],dH[Year][1:-1],"-",color=lcmap[i], linewidth=1)

    # X-tick label
    LatTick = arange(-35,35+0.01,5).astype(int)
    plt.xticks(LatTick)
    axplot.xaxis.set_ticklabels(LatTick)

    ## X-tick minor
    #axplot.xaxis.set_minor_locator(AutoMinorLocator())

    # Y-lim
    if lndsea=="lnd":
        axplot.set_ylim((3000,10000))
    elif lndsea=="sea":
        axplot.set_ylim((3000,6000))

    # Title
    plt.title(stitle, fontsize=8)


    # Legend file
    if legPath == None:
        return

    legPath = legPath
    lines   = axplot.get_lines()
    lname   = lYear
    figleg  = plt.figure(figsize=(1,5))
    figleg.legend(lines, lname)
    figleg.savefig(legPath)   
    plt.close()

    # Save
    plt.savefig(figPath)
    plt.close()
    print figPath





def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))

#-------------------------------
a2landfrac= load_landfrac()
dlndsea={"lnd" :ma.masked_not_equal(a2landfrac,1).mask
         ,"sea" :ma.masked_not_equal(a2landfrac,0).mask
         }


figDir  = "/tank/utsumi/PMM/HCELL/pict"
#lseason = ["MAM"]
lseason = ["DJF","MAM","JJA","SON"]
ddMon   = {"DJF":-1,"MAM":2,"JJA":5,"SON":8}
for season in lseason:

    dHlnd  = {}
    dHsea  = {}

    dMon= ddMon[season]
    for Year in range(iYear,eYear+1):
        iY,iM = util.shift_YM(Year,1,dMon)
        eY,eM = util.shift_YM(Year,1,dMon+2)
        a2num = ret_sumnum("num",[iY,iM],[eY,eM],"all")
        a2sum = ret_sumnum("sum",[iY,iM],[eY,eM],"all")

        a2h   = ma.masked_invalid(a2sum/a2num)
        dHlnd[Year] = ma.masked_where(dlndsea["lnd"], a2h).mean(axis=1)
        dHsea[Year] = ma.masked_where(dlndsea["sea"], a2h).mean(axis=1)

    sys.exit()
    # Figure
    for lndsea in llndsea:
        if lndsea=="sea":
            dDat = dHsea
        elif lndsea=="lnd":
            dDat = dHlnd

        stitle    = "Mean Height %s %s"%(lndsea,season)
        figPath   = os.path.join(figDir,"plot.LatH.%s.%s.png"%(lndsea,season))
        legPath   = os.path.join(figDir,"legend.plot.LatH.png")
    
        DrawLatH(dDat,stitle,figPath,legPath)
