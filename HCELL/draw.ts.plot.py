import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
import os, sys, socket
import matplotlib.pyplot as plt
from   scipy import stats
import hcell_func
import myfunc.grids as grids


iYM = [1997,12]
eYM = [2013,11]
#eYM = [2006,11]
dMon= 1
lYM = util.ret_lYM(iYM, eYM)[::dMon]
#var = "num"
var = "hgt"
hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

baseDir = os.path.join(rootDir, "HCELL/PDF/Type")
figDir  = "/tank/utsumi/PMM/HCELL/pict"

lregion = ["NAT","NAF","ASI","NPA","NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]
#lregion = ["NAT"]

#lrtype = ["strat","conv","other","all"]
lrtype = ["all"]
ny  = 148
nx  = 180
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)
llndsea = ["sea","lnd"]
#llndsea = ["sea"]

BinBnd = arange(0,12000+1,1000)
nz     = len(BinBnd)
liH    = [5,8,10]
#liH    = [10]

dhemi   = {}
for region in ["NAT","NAF","ASI","NPA","NAM"]:
    dhemi[region] = "N"
for region in ["SAT","SAF","SIN","OCE","SWP","SEP","SAM"]:
    dhemi[region] = "S"

#-------------------------------
# Function
#---------
def ret_a3sumnum_single(sumnum,Year,Mon,rtype):
    srcDir = os.path.join(baseDir, "%04d"%(Year))
    if rtype !="all":
        srcPath= os.path.join(srcDir, "%s.%s.%02d.%02dx%03dx%03d"%(sumnum,rtype,Mon, nz, ny, nx))
        a3in   = fromfile(srcPath, int32).reshape(nz,ny,nx)
    else:
        a3in   = ret_a3sumnum_single(sumnum,Year,Mon,"strat")\
                +ret_a3sumnum_single(sumnum,Year,Mon,"conv") \
                +ret_a3sumnum_single(sumnum,Year,Mon,"other")
    return a3in


def ret_a3sumnum(sumnum,iYM,eYM,rtype):
    lYM  = util.ret_lYM(iYM,eYM)
    a3dat= zeros([nz,ny,nx],int32)
    for Year,Mon in lYM:
       a3dat = a3dat + ret_a3sumnum_single(sumnum,Year,Mon,rtype)
    return a3dat


def load_landfrac():
    srcPath = os.path.join(rootDir,"data/const/landfrac.37SN.148x180")
    return fromfile(srcPath,float32).reshape(148,180)


def mk_regionMask(region,landsea):
    lndfracPath   = os.path.join(rootDir,"data/const/landfrac.37SN.148x180")
    a2lndfrac= fromfile(lndfracPath, float32).reshape(148,180)
    BBox = hcell_func.ret_regionBBox(region)
    a2region = grids.mk_mask_BBox(Lat,Lon,BBox,miss=-9999)
    if landsea == "lnd":
        a2region = ma.masked_where(a2lndfrac <0.9, a2region).filled(-9999.)
    elif landsea == "sea":
        a2region = ma.masked_where(a2lndfrac >0.1, a2region).filled(-9999.)
    elif landsea == "all":
        pass
    else:
        print "check landsea",landsea
        print "landsea must be","lnd","sea","all"
        sys.exit()
    return a2region


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
    figleg.legend(lines, lname)
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
    if var=="num":
        ax.set_ylim(top=a2in.max()*1.3)
    elif var=="hgt":
        print BinBnd[iH]
        ax.set_ylim(bottom=BinBnd[iH])
        ax.set_ylim(top=BinBnd[iH]+ (a2in.max() - BinBnd[iH])*1.3)


    # y-axis
    ax.tick_params(axis="y", which="major", pad=0)
    # text
    for iperiod,a1slp in enumerate(la1slp):
        for iline,slp in enumerate(a1slp):
            pv = la1pv[iperiod][iline]
            #plt.text(1998+8*iperiod,6000+100*iline,"%.1f (%.2f)"%(slp,pv), color=colors[iline])
            plt.text(0.1+0.4*iperiod,0.75+0.04*iline,"%.1f (%.2f)"%(slp,pv), color=colors[iline], transform=ax.transAxes)

    # Save main plot
    plt.savefig(figPath)
    print figPath



def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))

#*******************************
# Load Data
#-------------------------------

lkey = [(rtype,lndsea,region,iH)
        for rtype   in lrtype
        for lndsea  in llndsea
        for region  in lregion
        for iH      in liH
       ]
d2dat = {key:zeros([ny,len(lYM)],int32) for key in lkey}
d2num = {key:zeros([ny,len(lYM)],int32) for key in lkey}
d2sum = {key:zeros([ny,len(lYM)],int32) for key in lkey}

for it,(iYear,iMon) in enumerate(lYM):
    eYear,eMon = util.shift_YM(iYear,iMon,dMon-1)
    print iYear,iMon
    for rtype in lrtype:
        if var=="num":
            a3var = ret_a3sumnum("num",[iYear,iMon],[eYear,eMon], rtype)
        elif var=="hgt":
            a3num = ret_a3sumnum("num",[iYear,iMon],[eYear,eMon], rtype)
            a3sum = ret_a3sumnum("sum",[iYear,iMon],[eYear,eMon], rtype)
        else:
            print "check var",var
            sys.exit()

        for iH in liH:
            if var=="num":
                a2var= a3var[iH:,:,:].sum(axis=0)
            elif var=="hgt":
                a2num= a3num[iH:,:,:].sum(axis=0).astype(float32)
                a2sum= a3sum[iH:,:,:].sum(axis=0).astype(float32)
                a2var= ma.masked_invalid(a2sum/a2num)     # [m]
                #a2var= ma.masked_invalid(a2sum/a2num) / 1000.    # [km]
            else:
                print "check var",var
                sys.exit()

            lKey =[[lndsea,region] for lndsea in llndsea for region in lregion]
            for [lndsea,region] in lKey:
                a2regionMsk = mk_regionMask(region,lndsea)
                a1var       = ma.masked_where(a2regionMsk==miss, a2var).mean(axis=1).filled(0.0)
                d2dat[rtype,lndsea,region,iH][:,it] = a1var



#-- Seasonal ----
#lseason = ["DJF","MAM","JJA","SON"]
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
    # Line plot
    lKey =[[lndsea,region] for lndsea in llndsea for region in lregion]
    for [lndsea,region] in lKey:
        for iH in liH:    
            for season in lseason:
                ik    = dik[season] 
                if var=="num":
                    a2dat =  d2dat[rtype,lndsea,region,iH][:,ik::12]\
                            +d2dat[rtype,lndsea,region,iH][:,ik+1::12]\
                            +d2dat[rtype,lndsea,region,iH][:,ik+2::12]

                if var=="hgt":
                     a3dat = array([\
                             d2dat[rtype,lndsea,region,iH][:,ik::12]\
                            ,d2dat[rtype,lndsea,region,iH][:,ik+1::12]\
                            ,d2dat[rtype,lndsea,region,iH][:,ik+2::12]
                            ])
                     a2dat = ma.masked_equal(a3dat,0.0).mean(axis=0)
        
                lielat= [[5,10],[10,15],[15,20],[20,25],[25,30],[30,35]] 
                nx    = a2dat.shape[1]
                a2out = empty([len(lielat),nx])
    
                a2fit0= empty([len(lielat),16])
                a2fit1= empty([len(lielat),7])
    
                a1slp0= empty([len(lielat)])
                a1slp1= empty([len(lielat)])
    
                a1pv0 = empty([len(lielat)])
                a1pv1 = empty([len(lielat)])
    
                for iline,(ilat,elat) in enumerate(lielat):
                    if   dhemi[region]=="N":
                        iy    = ret_y(ilat)
                        ey    = ret_y(elat)
                    elif dhemi[region]=="S":
                        iy    = ret_y(-elat)
                        ey    = ret_y(-ilat)
                    else:
                        print "check hemisphere",region,dhemi[region]
                        sys.exit()

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
    
        
                figPath = os.path.join(figDir, "plot.TS.stormH.%s.%s.%s.%s.%s.%02dkm.png"%(var,season,rtype,lndsea,region,int(BinBnd[iH]/1000)))
                stitle  = "stormH (%s) @%s %s %s >%02dkm"%(rtype, season, lndsea,region,int(BinBnd[iH]/1000)) 
                legPath = os.path.join(figDir, "leg.plot.stormH.png")
                lname   = ["%.1f-%.1f"%(ilat,elat) for (ilat,elat) in lielat]
    #            PlotTime(a2out, lname, figPath, legPath, stitle)
                PlotTimeFit(a2out, [a2fit0,a2fit1], [arange(1998,2013+1),arange(2002,2008+1)], la1slp, la1pv, lname, figPath, legPath, stitle)
    
    
