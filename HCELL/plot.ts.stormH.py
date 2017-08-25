from numpy import *
import myfunc.util as util
import os, sys, socket
import matplotlib.pyplot as plt
import matplotlib
from   scipy import stats

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

baseDir = os.path.join(rootDir, "HCELL/SUMNUM/Type")

#lrtype = ["strat","conv","other","all"]
lrtype = ["all"]
ny  = 148
nx  = 180
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)
llndsea = ["sea","lnd"]
lhemi   = ["N","S"]
#lhemi   = ["S"]

lielatN= [[0,5],[5,10],[10,15],[15,20],[20,25],[25,30],[30,35]] 
lielatS= [[-x[1],-x[0]] for x in lielatN]
dlielat = {"N":lielatN, "S":lielatS}
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


def DrawTimeLat(a2in, figPath, stitle, cmap="viridis", vmin=None, vmax=None,season=False):
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


def PlotTimeFit(a2in, la2fit, la1xfit, lslp, lpv, lname, figPath, legPath, stitle):
    fig = plt.figure(figsize=(3,4.5))
    ax  = fig.add_axes([0.16,0.08,0.80,0.85])

    nline,nx = a2in.shape
    x   = range(iYM[0]+1,eYM[0]+1)

    [ax.plot(x,a2in[iline]) for iline in range(nline)]

    plt.xticks(x[::2])
    lxlabel = range(iYM[0]+1,eYM[0]+1,2)
    ax.xaxis.set_ticklabels(lxlabel, fontsize=8)

    # Title
    plt.title(stitle,fontsize=9)


    # Legend file
    if legPath == None:
        return

    lines   = ax.get_lines()
    figleg  = plt.figure(figsize=(2,2))
    figleg.legend(lines[::-1], lname[::-1])
    figleg.savefig(legPath)   
    plt.close()

    # plot fit line
    colors = [line.get_color() for line in lines]
    #colors = ["r","b","g","k","yellow"]
    dlinestyle = {0:"--",1:":"}
    for iperiod,a2fit in enumerate(la2fit):
        for iline,a1fit in enumerate(a2fit):
            x = la1xfit[iperiod]
            ax.plot(x,a1fit,color=colors[iline],linewidth=1,linestyle=dlinestyle[iperiod])
#    

    # set ylim
    if   lndsea=="sea":
        ymin = 3200
        ymax = 6600
        textpos = 5500
    elif lndsea=="lnd":
        ymin = 3600
        ymax = 7200
        textpos = 3700


    ax.set_ylim((ymin,ymax))

    # text
    for iperiod,a1slp in enumerate(la1slp):
        for iline,slp in enumerate(a1slp):
            pv = la1pv[iperiod][iline]
            plt.text(1999+8*iperiod,textpos+120*iline,"%.1f (%.2f)"%(slp,pv), color=colors[iline])

    # Save main plot
    plt.savefig(figPath)
    print figPath



def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))

#-------------------------------
a2landfrac= load_landfrac()
da2lndsea={"lnd":ma.masked_not_equal(a2landfrac,1).mask
         ,"sea" :ma.masked_not_equal(a2landfrac,0).mask
         }


d2dat = {(rtype,lndsea):zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}
d2num = {(rtype,lndsea):zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}
d2sum = {(rtype,lndsea):zeros([ny,len(lYM)],int32) for rtype in lrtype for lndsea in llndsea}

for it,(iYear,iMon) in enumerate(lYM):
    eYear,eMon = util.shift_YM(iYear,iMon,dMon-1)
    print iYear,iMon
    print eYear,eMon
    for rtype in lrtype:
        for lndsea in llndsea:
            a2lndsea = da2lndsea[lndsea]
            a2num = ret_sumnum("num",[iYear,iMon],[eYear,eMon], rtype)
            a2sum = ret_sumnum("sum",[iYear,iMon],[eYear,eMon], rtype)

            # mask lndsea
            a2num = ma.masked_where(a2lndsea, a2num)
            a2sum = ma.masked_where(a2lndsea, a2sum)
            a1num= a2num.mean(axis=1) 
            a1sum= a2sum.mean(axis=1) 
               
            a1dat= ma.masked_invalid(a1sum.astype(float32)/a1num).filled(0) 
    
            d2sum[rtype,lndsea][:,it] = a1sum
            d2num[rtype,lndsea][:,it] = a1num
            d2dat[rtype,lndsea][:,it] = a1dat


figDir  = "/tank/utsumi/PMM/HCELL/pict"
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

for hemi in lhemi:
    lielat = dlielat[hemi]
    print lielat
    for rtype in ["all"]:
        # Line plot
        for lndsea in llndsea:
            for season in lseason:
                ik    = dik[season] 
                a2sum =   d2sum[rtype,lndsea][:,ik::12]\
                         +d2sum[rtype,lndsea][:,ik+1::12]\
                         +d2sum[rtype,lndsea][:,ik+2::12]
        
                a2num =   d2num[rtype,lndsea][:,ik::12]\
                         +d2num[rtype,lndsea][:,ik+1::12]\
                         +d2num[rtype,lndsea][:,ik+2::12]
        
                a2dat = ma.masked_invalid(a2sum/a2num)
        
                nx    = a2dat.shape[1]
                a2out = empty([len(lielat),nx])
    
                a2fit0= empty([len(lielat),16])
                a2fit1= empty([len(lielat),7])
    
                a1slp0= empty([len(lielat)])
                a1slp1= empty([len(lielat)])
    
                a1pv0 = empty([len(lielat)])
                a1pv1 = empty([len(lielat)])
    
                for iline,(ilat,elat) in enumerate(lielat):
                    iy    = ret_y(ilat)
                    ey    = ret_y(elat)
                    print "iy,ey=",iy,ey
                    a1out = a2dat[iy:ey].mean(axis=0)
                    a2out[iline] = a1out
    
                    # Linear fit for all years
                    x     = arange(1998,2013+1)
                    y     = a1out # 1998-2013
    
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                    print lndsea,ilat,elat,slope, p_value
                    a2fit0[iline] = x*slope + intercept
                    a1slp0[iline] = slope
                    a1pv0 [iline] = p_value
    
                    # Linear fit for 2002-2008
                    x     = arange(2002,2008+1)
                    y     = a1out[4:4+7] # 2002-2008
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                    print lndsea,ilat,elat,slope, p_value
                    a2fit1[iline] = x*slope + intercept
                    a1slp1[iline] = slope
                    a1pv1 [iline] = p_value
    
                la2fit = [a2fit0, a2fit1]
                la1xfit= [arange(1998,2013+1),arange(2002,2008+1)]
                la1slp = [a1slp0, a1slp1]
                la1pv  = [a1pv0,  a1pv1]
    
        
                figPath = os.path.join(figDir, "plot.ts.stormH.%s.%s.%s.%s.png"%(hemi,season,rtype,lndsea))
                stitle  = "stormH %s-hemi (%s) @%s %s 2A23"%(hemi, rtype, season, lndsea)
        
                legPath = os.path.join(figDir, "leg.plot.stormH.%s.png"%(hemi))
                lname   = ["%.1f-%.1f"%(ilat,elat) for (ilat,elat) in lielat]
    #            PlotTime(a2out, lname, figPath, legPath, stitle)
                PlotTimeFit(a2out, [a2fit0,a2fit1], [arange(1998,2013+1),arange(2002,2008+1)], la1slp, la1pv, lname, figPath, legPath, stitle)
    
    
