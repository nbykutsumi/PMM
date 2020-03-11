import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket

myhost = socket.gethostname()
if myhost == 'shui':
    workbaseDir= '/work'
    tankbaseDur= '/tank'
    listDir    = '/work/hk01/utsumi/PMM/US/obtlist'
    srcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()


expr = 'glb.wprof.org'

argv = sys.argv
if len(argv)==2:
    largv = argv[1].split()
    dargv = {}
    for argv in largv:
        key,param = argv.split('=')
        dargv[key] = param

    print dargv
    oid    = int(dargv['oid'])
    Year   = int(dargv['Year'])
    Mon    = int(dargv['Mon'])
    Day    = int(dargv['Day'])
    iy     = int(dargv['iy'])
    ey     = int(dargv['ey'])
    clat   = float(dargv['clat'])
    clon   = float(dargv['clon'])
    DB_MAXREC = int(dargv['DB_MAXREC'])
    dlatlon= int(dargv['dlatlon'])
    xpos   = int(dargv['xpos'])

    if clon >=180:
        clon = -180 + (clon-180)


elif len(argv)==1:

    ## Africa case
    #oid = 2421
    #iy, ey = 2029, 2089
    #clat    = 14 # Africa case
    #clon    = 2  # -180 - +180
    
    # SE.US case, oid=003556, 2014/10/14
    oid = 3556
    Year,Mon,Day = 2014,10,14
    iy, ey = -9999,-9999
    #iy, ey = 987, 1047
    #iy,ey = 917,1117
    #clat    = 34    # SE.US case. oid = 003556
    clat    = 31    # SE.US case. oid = 003556
    clon    = -86   # 2014/10/14  05:42:03 UTC
    DB_MAXREC = 10000
    DB_MINREC = 1000
    expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    #expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    xpos    = 100  # x-position for cross section


    #
    ## SW.Japan typhoon case, oid=019015, 2017/07/03
    #oid = 19015
    #Year,Mon,Day = 2017,7,3
    #iy, ey = 1793, 1973
    #clat    = 33    # SE.US case. oid = 003556
    #clon    = 130   # 2014/10/14  05:42:03 UTC
    #DB_MAXREC = 20000
    
    
    ## QJRMS case, oid=012149, 2016/4/18
    #oid = 12149
    #iy, ey = 1038, 1098
    #clat    = 32. # QJRMS case, oid=012149
    #clon    = -94 # -180 - +180
    #DB_MAXREC = 20000
    #Year,Mon,Day = 2016,4,18

    dlatlon = 5
    #dscan   = 30
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
    [[lllat,lllon],[urlat,urlon]] = BBox

miss = -9999.
#--------------------------------
#*********************************************
def ret_domain_cy(a2lat, a2lon, clat, clon, dlatlon):
    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,nxTmp/2]
    a1lon = a2lon[:,nxTmp/2]

    idx_latmax = np.argmax(a1lat)
    a1lat0 = a1lat[:idx_latmax+1]
    a1lat1 = a1lat[idx_latmax+1:]
    a1lon0 = a1lon[:idx_latmax+1]
    a1lon1 = a1lon[idx_latmax+1:]

    #-- search first half: ascending --
    found_domain = 0
    idx_c  = bisect_left(a1lat0, clat)
    latTmp = a1lat0[idx_c]
    lonTmp = a1lon0[idx_c]
    if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
        found_domain = 1
    else:
        #-- search second half: descending --
        idx_c  = bisect_left(a1lat1[::-1], clat)
        idx_c  = len(a1lat) - idx_c -1
        latTmp = a1lat[idx_c]
        lonTmp = a1lon[idx_c]

        if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
            found_domain =1

    if found_domain==1:
        return idx_c

    else:
        print 'No matching scans in the target domain are found.'
        print 'Exit'
        sys.exit()

#--------------------------------
def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr,nzdpr= a3in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,miss).mean(axis=0)
    a2datTmp[a1dprmask,:] = miss
    return a2datTmp

#--------------------------------

def average_2ranges_3d(a3in,miss=None,dtype=float32, fill=True):
    '''
    a3in: (ny, nx, nz) --> a2out: (ny, nx, nz/2)
    nz should be an even number
    '''
    ny,nx,nz = a3in.shape
    a4out = empty([2,ny,nx,nz/2], dtype)
    a4out[0] = a3in[:,:,::2]
    a4out[1] = a3in[:,:,1::2]
    if fill==True:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).filled(miss).astype(dtype)
    else:
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).astype(dtype)
    return a3out


#--------------------------------
def load_GMI_Tb13(Year,Mon,Day,oid):
    #-- Corresponding S2 yx ----------
    #srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    srcDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    a1xgmi2 = ma.masked_less(np.load(srcDir + '/Xpy.1.%06d.npy'%(oid)).flatten(), 0)
    a1ygmi2 = ma.masked_less(np.load(srcDir + '/Ypy.1.%06d.npy'%(oid)).flatten(), 0)

    a1flag  = a1xgmi2.mask + a1ygmi2.mask
    
    a1xgmi2 = a1xgmi2.filled(0) 
    a1ygmi2 = a1ygmi2.filled(0) 
    #-- Read HDF -----
    gmiDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
    #ssearch= gmiDir + '/1C.GPM.GMI.XCAL2016-C.20141014-S081337-E094609.003558.V05A.HDF5'
    ssearch= gmiDir + '/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(oid)
    gmiPath= glob.glob(ssearch)[0]
    with h5py.File(gmiPath, 'r') as h:
        a3tbNS1    = h['/S1/Tc'][:]
        a3tbNS2Org = h['/S2/Tc'][:]

    nytmp, nxtmp, nztmp = a3tbNS1.shape 
    a1tmp  = a3tbNS2Org[a1ygmi2, a1xgmi2,:]
    a1tmp[a1flag] = -9999. 
    a3tbNS2 = a1tmp.reshape(nytmp, nxtmp, -1)
     
    return np.concatenate([a3tbNS1, a3tbNS2], axis=2)


#*****************
#- Read My data ----
#srcDir = '/home/utsumi/temp/out/my'
#stamp = '%06d.y%04d-%04d.nrec%d'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
srcDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)
#a2topzmMS  = np.load(srcDir + '/top-zmMS.%s.npy'%(stamp))[:, xpos, :]
a2topzmNS  = np.load(srcDir + '/top-zmNS.%s.npy'%(stamp))[:, xpos, :]
#a2topprNS  = np.load(srcDir + '/top-prprofNS.%s.npy'%(stamp))[:, xpos, :]
#a2topprNScmb= np.load(srcDir + '/top-prprofNScmb.%s.npy'%(stamp))[:, xpos, :]
a2toptbNS  = np.load(srcDir + '/top-tbNS.%s.npy'%(stamp))[:, xpos, :][:,::-1]  # Flip

#a2prNS   = np.load(srcDir + '/prprofNS.%s.npy'%(stamp))[:, xpos, :]
#a2prNScmb= np.load(srcDir + '/prprofNScmb.%s.npy'%(stamp))[:,xpos,:]

a2latMy  = np.load(srcDir + '/lat.%s.npy'%(stamp))
a2lonMy  = np.load(srcDir + '/lon.%s.npy'%(stamp))
a1latMy  = a2latMy[:,xpos]
a1lonMy  = a2latMy[:,xpos]

#--- Find iyTmp & eyTmp for drawing if iy==ey==-9999 ---
if ((iy<0) or (ey<0)):
    cy = ret_domain_cy(a2latMy, a2lonMy, clat, clon, dlatlon)
    #dscan = 90
    dscan = 50

    print a2latMy.shape[0],cy,dscan
    iyTmp = max(0, cy-dscan)
    eyTmp = min(a2latMy.shape[0], cy+dscan)
    iyMy  = iyTmp
    eyMy  = eyTmp

else:
    iyTmp = iy
    eyTmp = ey
    iyMy = 0
    eyMy = ey-iy



#***********************************
# Extract domain from my data
#***********************************
a2topzmNS = a2topzmNS[iyMy:eyMy+1,:]
#a2topzmMS = a2topzmMS[iyMy:eyMy+1,:]
#a2topprwatMS = a2topprwatNS*0.0 -9999.

a1latMy  = a2latMy[iyMy:eyMy+1,:]
a1lonMy  = a2latMy[iyMy:eyMy+1,:]

a2toptbNS= a2toptbNS[iyMy:eyMy+1,:]
#***********************************
#- Read GMI 1C (Tb) ----
#-----------------------------------
a3tb = load_GMI_Tb13(Year,Mon,Day,oid)
a2tb = a3tb[iyTmp:eyTmp+1,xpos,:][:,::-1]
#***********************************
#- Read DPR data ----
#-- Cooresponding DPR yx --------
#srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
srcDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1xdpr = np.load(xPath)[iyTmp:eyTmp+1,xpos-83]
a1ydpr = np.load(yPath)[iyTmp:eyTmp+1,xpos-83]

#-- Read DPR (Ku) ----
#dprDir  = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/%04d/%02d/%02d'%(Year,Mon,Day)
if Year in [2014,2015]:
    dprDir  = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06/%04d/%02d/%02d'%(Year,Mon,Day)
else:
    dprDir  = workbaseDir + '/hk01/PMM/NASA/GPM.Ku/2A/V06/%04d/%02d/%02d'%(Year,Mon,Day)

dprPath = glob.glob(dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid))[0]


with h5py.File(dprPath, 'r') as h:
    a3dprzmNS = h['NS/PRE/zFactorMeasured'][:]
    #a3dprprNS = h['NS/SLV/precipRate'][:]
    a2latdpr = h['NS/Latitude'][:]   # test
    a2londpr = h['NS/Longitude'][:]  # test
    #a2cfbBin = h['NS/PRE/binClutterFreeBottom'][:]
    #a2surfBin= h['NS/PRE/binRealSurface'][:]


#-- test --
print 'latdpr-extracted'
print a2latdpr[a1ydpr,a1xdpr]
print 'londpr-extracted'
print a2londpr[a1ydpr,a1xdpr]
print 'shape'
print a2latdpr.shape
print
#- extract lower levels --
a3dprzmNS = a3dprzmNS[:,:,-88*2:]   # 125 m layers
#a3dprzmMS = a3dprzmMS[:,:,-88*2:]  # 125 m layers

#- average every 2-levels (125m -> 250m layer) --
a3dprzmNS = average_2ranges_3d(a3dprzmNS, miss=-9999.9)


##- Extract and average 9-grids Zm ---
#_,_,nztmp = a3dprzmNS.shape 
#a2dprzmNS = ave_9grids_3d(a3dprzmNS, a1ydpr, a1xdpr, miss=-9999.9).reshape(-1, nztmp)
#
#_,_,nztmp = a3dprprNS.shape 
#a2dprprNS = ave_9grids_3d(a3dprprNS, a1ydpr, a1xdpr, miss=-9999.9).reshape(-1, nztmp)


#- Extract Radar pixels that correspond GMI ---
#a2dprprNS = a3dprprNS[a1ydpr,a1xdpr]
a2dprzmNS = a3dprzmNS[a1ydpr,a1xdpr]
#- Screen dB less than 15dB 
#a2topzmMS = (ma.masked_less(a2topzmMS, 1500)*0.01).filled(miss)  # 50x250m leyers
a2topzmNS = (ma.masked_less(a2topzmNS, 1500)*0.01).filled(miss)  # 50x250m layers
#a2dprzmMS = ma.masked_less(a2dprzmMS,15).filled(miss)
a2dprzmNS = ma.masked_less(a2dprzmNS,15).filled(miss)[:,-50:]

#******************************************************
# Figure Zm and Tb
#******************************************************
fig   = plt.figure(figsize=(12,10))
ssize = 1

#nx    = a2prNS.shape[0]
#nbin  = a2prNS.shape[1]
nx    = a2dprzmNS.shape[0]
nbin  = a2dprzmNS.shape[1]
aspect1 = 0.16*nx/nbin
aspect2 = aspect1*5

zmin,  zmax  = 15, 45
tbmin, tbmax = 140, 300
tick_locs88 = [88, 80, 72, 64, 56, 48, 40, 32, 24, 16, 8, 0]
tick_lbls88 = [' 0',' 2',' 4',' 6',' 8','10','12','14','16','18','20','22']
tick_locs64 = [64, 56, 48, 40, 32, 24, 16, 8, 0]
tick_lbls64 = [' 0',' 2',' 4',' 6',' 8','10','12','14','16']
tick_locs50 = [2, 10, 18, 26, 34, 42]  # 0.25km x 50
tick_lbls50 = [12,10, 8, 6, 4, 2]  # 0.25km x 50



tick_locsTb = [None]
tick_lblsTb = [None]
#for i in range(6):
#for i in range(8):
for i in range(4):
    nx,nh = a2dprzmNS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'

    #--- Zm ---------------------------- 
    x0= 0.1
    h = 0.2
    w = 0.7

    #--- Tb ---------------------------- 
    if i==0:
        y0 = 0.02
        a2dat = ma.masked_less_equal(a2tb.T,0)
        vmin, vmax = tbmin, tbmax
        tick_locs, tick_lbls = tick_locsTb, tick_lblsTb
        aspect= aspect2
        stype = 'Obs Tb(K)'
        cbarlbl='Kelvin'

    elif i==1:
        y0 = 0.25
        a2dat = ma.masked_less_equal(a2toptbNS.T,0)
        vmin, vmax = tbmin, tbmax
        tick_locs, tick_lbls = tick_locsTb, tick_lblsTb
        aspect= aspect2
        stype = 'Top-Ranked Tb(K) Ret=NS'
        cbarlbl='Kelvin'

    elif i==2:
        y0 = 0.48
        print 'a2dprzmNS.shape', a2dprzmNS.shape
        a2dat = ma.masked_less_equal(a2dprzmNS.T,0)
        vmin, vmax = zmin, zmax
        tick_locs, tick_lbls = tick_locs50, tick_lbls50
        aspect= aspect1
        stype = 'Zm DPR-Ku'
        cbarlbl='dBZ (uncorrected)'

    elif i==3:
        y0 = 0.71
        a2dat = ma.masked_less_equal(a2topzmNS.T,0)
        vmin, vmax = zmin, zmax
        tick_locs, tick_lbls = tick_locs50, tick_lbls50
        aspect= aspect1
        stype = 'Top-Ranked Zm (Ku)'
        cbarlbl='dBZ (uncorrected)'


    ax = fig.add_axes([x0, y0, w, h])
    cax= fig.add_axes([x0+w+0.01, y0, 0.02, h*0.9])
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.grid()
    ax.set_title(stype+' '+stamp)
    ax.set_yticks(tick_locs)
    ax.set_yticklabels(tick_lbls)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl)

    if i in [0,1]:
        iscan1 = 0
        iscan2 = a2dat.shape[0]-1
        xoff = (iscan2 - iscan1)/14
        #xoff = 0
        ax.text(xoff,12, r'10.7V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff,11, r'10.7H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff,10, r'18.7V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 9, r'18.7H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 8, r'23.8V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 7, r'36.5V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 6, r'36.5H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 5, r'89.0V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 4, r'89.0H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 3, r'166V'         , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 2, r'166H'         , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 1, r'183.31$\pm$3V', family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        ax.text(xoff, 0, r'183.31$\pm$7V', family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')

##------------
outPath  = figDir + '/prof.zm.tb.%06d.y%04d-%04d.nrec%d.png'%(oid,iy,ey,DB_MAXREC)
plt.savefig(outPath)
print outPath
plt.clf()


