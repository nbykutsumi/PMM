import matplotlib
matplotlib.use("Agg")
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import myfunc.util as util
import GPMGV
import sys

def load_gtopo(lat,lon):
    orogDir  = "/data1/hjkim/GTOPO30"

    ullat = int( (lat - (-60))/50. )*50. + 50 -60.
    ullon = int( (lon - (-180))/40.)*40. -180.

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<ullon)&(ullon<0):
        WE = "W" 
    else:
        WE = "E" 

    #orogPath = orogDir + "/E060N40.DEM"
    orogPath = orogDir + "/%s%03d%s%02d.DEM"%(WE, abs(ullon), SN, abs(ullat))
    dmwPath  = orogDir + "/%s%03d%s%02d.DMW"%(WE, abs(ullon), SN, abs(ullat))   
    print dmwPath
    """
    BYTEORDER      M
    LAYOUT       BIL
    NROWS         6000
    NCOLS         4800
    NBANDS        1
    NBITS         16
    BANDROWBYTES         9600
    TOTALROWBYTES        9600
    BANDGAPBYTES         0
    NODATA        -9999
    ULXMAP        60.00416666666667
    ULYMAP        39.99583333333333
    XDIM          0.00833333333333
    YDIM          0.00833333333333
    """
    
    ny, nx = 6000, 4800
    a = flipud(fromfile(orogPath, "int16").byteswap().reshape(ny,nx))

    # load DMW file
    f=open(dmwPath, "r"); lines=f.readlines(); f.close()
    lonmin = float(lines[4])
    latmax = float(lines[5])
    lonmax = lonmin + 0.00833333333333*(nx-1)
    latmin = latmax - 0.00833333333333*(ny-1)

    dlat = 0.00833333333333
    dlon = 0.00833333333333
    para = {}
    para["ny"] = ny
    para["nx"] = nx
    para["lllat"]=latmin
    para["lllon"]=lonmin
    para["urlat"]=latmax
    para["urlon"]=lonmax
    para["dlat" ]=dlat
    para["dlon" ]=dlon

    para["Lat"]  = arange(latmin,latmax+0.5*dlat, dlat)
    para["Lon"]  = arange(lonmin,lonmax+0.5*dlon, dlon)

    return a,para
 


gv = GPMGV.GPMGV()
#gv.load_sitelist()
gv.load_sitelist_reclassified()

ldomain = gv.domains
#for domain in ldomain[5:]:
for domain in ldomain:

    #if domain not in ["KWAJALEIN-KWA","KWAJALEIN-RMI"]: continue
    #if domain not in ["VIRGINIA-NSWD-N"]: continue
    #if domain not in ["VIRGINIA-NASA-SE","VIRGINIA-NASA-NE"]: continue
    dgNames = gv.dgNames
    lkey = [ (domain, gName) for gName in gv.dgNames[domain]]
    
    llatlon = []
    for key in lkey:
        latlon = gv.dlatlon[key]
        llatlon.append(latlon)
    
    alatlon = array(llatlon)
    latmax   = alatlon[:,0].max()
    latmin   = alatlon[:,0].min()
    lonmax   = alatlon[:,1].max()
    lonmin   = alatlon[:,1].min()
    
    dlat = latmax - latmin
    dlon = lonmax - lonmin


    wlatlon = max(dlat, dlon)
    if wlatlon ==0:
        wlatlon = 0.5
    urlat = latmax + wlatlon*0.2
    lllat = latmin - wlatlon*0.2
    urlon = lonmax + wlatlon*0.2
    lllon = lonmin - wlatlon*0.2
    
    
    figmap = plt.figure(figsize=(4,4))
    axmap  = figmap.add_axes([0.2,0.2,0.7,0.7])
    M      = Basemap(resolution="i", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax= axmap)

    #- load DEM

    if domain in ["CALIFORNIA-ERK"]:
        clat   = (urlat+lllat)*0.5
        clon   = (urlon+lllon)*0.5
        a2dem0, para0 = load_gtopo(clat,clon)
        Lat0 = para0["Lat"]
        Lon0 = para0["Lon"]

        a2dem1, para1 = load_gtopo(lllat,clon)
        Lat1 = para1["Lat"]
        Lon1 = para1["Lon"]

        a2dem = concatenate([a2dem1,a2dem0],axis=0)
        Lat   = concatenate([Lat1, Lat0])
        Lon   = Lon0
   
    else:
        clat   = (urlat+lllat)*0.5
        clon   = (urlon+lllon)*0.5
        a2dem, para = load_gtopo(clat,clon)
        Lat = para["Lat"]
        Lon = para["Lon"]

    #- draw elevation
    X,Y = meshgrid(Lon, Lat)

    a2dem = ma.masked_equal(a2dem,-9999 )
    m = M.pcolormesh(X,Y,a2dem, vmin=0, vmax=6000, cmap="terrain")

    #- colorbar
    M.colorbar(m)

    #- plot sites
    Lat, Lon = zip(*llatlon)
    M.plot(Lon,Lat, "o", color="r")
     
    #- coastlines
    M.drawcoastlines()
    dtick = int(wlatlon*1000/5) /1000.
    #- parallels & meridians
    if dtick ==0:
        parallels = arange(lllat,urlat,wlatlon/4.)
        meridians = arange(lllon,urlon,wlatlon/4.)
    else:
        parallels = arange(-90,90,dtick)
        meridians = arange(-180,180,dtick)
    M.drawparallels(parallels,labels=[1,0,0,0])
    M.drawmeridians(meridians, labels=[0,0,0,1], rotation=60)
   
    #- title
    stitle  = "%s"%(domain)
    plt.title(stitle)

    figDir  = "/work/a01/utsumi/GPMGV/fig"
    figPath = figDir + "/sites.recl.%s.png"%(domain)
    plt.savefig(figPath)
    print figPath
