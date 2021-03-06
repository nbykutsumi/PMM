import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
import myfunc.grids as grids
import os, sys, socket
import matplotlib.pyplot as plt
import hcell_func



iYM = [1997,12]
eYM = [2013,11]
nYear=eYM[0]-(iYM[0]+1)+1
lYM = util.ret_lYM(iYM, eYM)



hostname = socket.gethostname()
if  hostname =="mizu":
    rootDir   = "/home/utsumi/mnt/wellshare"
elif hostname=="well":
    rootDir   = "/media/disk2/share"

baseDir = os.path.join(rootDir, "HCELL/PDF/Type")

Bin = arange(500,12500+1,1000)
BinBnd = arange(0,12000+1,1000)

#lrtype = ["all","strat","conv","other"]
lrtype = ["all"]
ny  = 148
nx  = 180
nz  = len(Bin)
miss= -9999.
Lat = arange(-37+0.25,37-0.25+0.01,0.5)
Lon = arange(0+1.0,360-1.0+0.01,2.0)



lielat= [(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35)]  # degree
#lielat= [(20,25)]  # degree
#-------------------------------
# Function
#---------
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

def ret_num_single(Year,Mon,rtype):
    srcDir = os.path.join(baseDir, "%04d"%(Year))
    if rtype !="all":
        srcPath= os.path.join(srcDir, "num.%s.%02d.%02dx%03dx%03d"%(rtype,Mon, nz, ny, nx))
        a3in   = fromfile(srcPath, int32).reshape(nz,ny,nx)
    else:
        a3in   = ret_num_single(Year,Mon,"strat")\
                +ret_num_single(Year,Mon,"conv") \
                +ret_num_single(Year,Mon,"other")
    return a3in

def ret_num(iYM,eYM,rtype):
    lYM  = util.ret_lYM(iYM,eYM)
    a3dat= zeros([nz,ny,nx],int32)
    for Year,Mon in lYM:
       a3dat = a3dat + ret_num_single(Year,Mon,rtype)
    return a3dat

def ret_y(lat):
    lat_first = -37.0
    dlat      = 0.5
    return int(floor((lat - lat_first)/dlat))



def DrawTimeH(a2in, figPath, stitle, cmap=None, vmin=None, vmax=None):
    fig = plt.figure(figsize=(3.3,1.5))
    ax  = fig.add_axes([0.09,0.13,0.90,0.7])

    ny,nx = a2in.shape
    x   = range(nx)
    y   = Bin
    X,Y = meshgrid(x,y)

    im  = ax.pcolormesh(X,Y,a2in,cmap=cmap,vmin=vmin,vmax=vmax)

    # Colorbar
    cb  = plt.colorbar(im, orientation="vertical")
    cb.ax.tick_params(labelsize=8)

    # X-tick label
    plt.xticks(x[::2])
    lxlabel = range(iYM[0]+1,eYM[0]+1)[::2]
    ax.xaxis.set_ticklabels(lxlabel, fontsize=8)

    # Y-tick label
    yticks  = BinBnd[::2]
    plt.yticks(yticks)
    ax.yaxis.set_ticklabels(yticks/1000, fontsize=8)

    # lines at 2001 & 2009
    ax.plot([2002-1998,2002-1998],[Bin[0],Bin[-1]],"-",color="k")
    ax.plot([2009-1998,2009-1998],[Bin[0],Bin[-1]],"-",color="k")


    # Title
    plt.title(stitle,fontsize=9)

    # Save
    plt.savefig(figPath)
    plt.close()
    print figPath

#-------------------------------
llandsea= ["sea","lnd"]
lregion = ["NAT","NAF","ASI","NPA","NAM","SAT","SAF","SIN","OCE","SWP","SEP","SAM"]
#lregion = ["SAF","SIN","OCE","SWP","SEP","SAM"]

lseason = ["DJF","MAM","JJA","SON"]
dik     = {"DJF":0,"MAM":3,"JJA":6,"SON":9}

dhemi   = {}
for region in ["NAT","NAF","ASI","NPA","NAM"]:
    dhemi[region] = "N"
for region in ["SAT","SAF","SIN","OCE","SWP","SEP","SAM"]:
    dhemi[region] = "S"

lkey = [[season,rtype,landsea] 
                       for season in lseason
                       for rtype  in lrtype
                       for landsea in llandsea
                        ]
#-------------------
# Load data
for [season,rtype,landsea] in lkey:
    ik  = dik[season]
    d2hist  = {(region,ielat):zeros([nz,nYear],int32) 
                for region  in lregion
                for ielat   in lielat}

    d2region = {region: ma.masked_equal(mk_regionMask(region, landsea),miss) for region in lregion}
 
    for it,(iYear,iMon) in enumerate(lYM[ik::12]):
        eYear,eMon = util.shift_YM(iYear,iMon,2)
        a3in = ret_num([iYear,iMon],[eYear,eMon], rtype)
        for region in lregion:
            a2regionMask = d2region[region]

            for ilat,elat in lielat:
                ielat = (ilat,elat)
                if dhemi[region]=="N":
                    iy,ey = ret_y(ilat), ret_y(elat)
                elif dhemi[region]=="S":
                    iy,ey = ret_y(-elat), ret_y(-ilat)
                    #print ielat,"iy,ey=",iy,ey
                for iz in range(nz):
                    d2hist[region,ielat][iz,it] = ma.masked_where(a2regionMask.mask, a3in[iz])[iy:ey+1,:].mean()
                    #print d2hist[region,ielat][iz,it]
#-------------------   
    
    # Counts
    for region in lregion:
        for ielat in lielat:
            ilat,elat = ielat
            stitle = "stormH %s %s %s %s Lat:%.1f-%.1f deg"%(landsea, region, rtype, season, ilat, elat)
            figDir = "/tank/utsumi/PMM/HCELL/pict"
            figPath= os.path.join(figDir,"TS.Count.%s.%s.%s.%s.Lat.%.1f.%.1f.png"%(landsea,region,rtype,season,ilat,elat))

            #print d2hist[(region,ielat)]   
            DrawTimeH(d2hist[(region,ielat)], figPath, stitle)


    # Normalized Counts
    for region in lregion:
        for ielat in lielat:
            ilat,elat = ielat
    
            a2norm    = empty([nz,nYear])
            a1clim    = d2hist[(region,ielat)][:,2002-1998:2008-1998+1].mean(axis=1)
            for iz in range(nz):
                a2norm[iz] = (d2hist[(region,ielat)][iz]-a1clim[iz])/a1clim[iz]
    
            a2norm = ma.masked_invalid(a2norm) 
            stitle = "Norm.H %s %s %s %s Lat:%.1f-%.1f deg"%(landsea,region, rtype, season, ilat, elat)
            figDir = "/tank/utsumi/PMM/HCELL/pict"
            figPath= os.path.join(figDir,"TS.Count.Norm.%s.%s.%s.%s.Lat.%.1f.%.1f.png"%(landsea,region,rtype,season,ilat,elat))
            cmap   = "RdBu_r" 
            DrawTimeH(a2norm, figPath, stitle, cmap,vmin=-1.0, vmax=1.0)

    
    """
    # PDF
    for ielat in lielat:
        ilat,elat = ielat


        a2pdf    = empty([nz,nYear])
        a1sum    = d2hist[ielat].sum(axis=0).astype(float32)
        for it in range(nYear):
            a2pdf[:,it] = d2hist[ielat][:,it]/a1sum[it]
        a2out    = a2pdf

        stitle = "stormH %s %s %s Lat:%.1f-%.1f deg"%(landsea, rtype, season, ilat, elat)
        figDir = "/tank/utsumi/PMM/HCELL/pict"
        figPath= os.path.join(figDir,"TS.PDF.%s.%s.%s.Lat.%.1f.%.1f.png"%(landsea,rtype,season,ilat,elat))
    
        DrawTimeH(a2out, figPath, stitle)
    """

    """
    # Anomary of PDF
    for ielat in lielat:
        ilat,elat = ielat


        a2pdf    = empty([nz,nYear])
        a1sum    = d2hist[ielat].sum(axis=0).astype(float32)
        for it in range(nYear):
            a2pdf[:,it] = d2hist[ielat][:,it]/a1sum[it]

        a2norm   = empty([nz,nYear])
        a1clim   = a2pdf.mean(axis=1)
        for iz in range(nz):
            a2norm[iz] = (a2pdf[iz]-a1clim[iz]) / a1clim[iz]
        a2norm = ma.masked_invalid(a2norm)

        stitle = "H RelF %s %s %s Lat:%.1f-%.1f deg"%(landsea, rtype, season, ilat, elat)
        figDir = "/tank/utsumi/PMM/HCELL/pict"
        figPath= os.path.join(figDir,"TS.PDF.Norm.%s.%s.%s.Lat.%.1f.%.1f.png"%(landsea,rtype,season,ilat,elat)) 
        cmap   = "RdBu_r" 
        DrawTimeH(a2norm, figPath, stitle, cmap,vmin=-1.0, vmax=1.0)

    """

