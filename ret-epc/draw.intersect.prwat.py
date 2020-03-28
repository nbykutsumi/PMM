# %%
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
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    srcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()



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
    #expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    xpos    = 100  # x-position for cross section
    #xpos    = 107  # x-position for cross section

    ## Colorado case
    #oid = 1574
    #Year,Mon,Day = 2014,6,8
    #iy, ey = -9999,-9999
    #clat    = 40    # SE.US case. oid = 003556
    #clon    = 256-360   # 2014/10/14  05:42:03 UTC
    #DB_MAXREC = 10000
    #DB_MINREC = 1000
    ##expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    #expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
    ##xpos    = 100  # x-position for cross section
    #xpos    = 115  # x-position for cross section


    ## Europe just to check batch-version
    #oid = 1927
    #Year,Mon,Day = 2014,7,1
    #iy, ey = 1764, 1784
    #clat    = 50.8
    #clon    = 28.1
    #DB_MAXREC = 10000
    ##DB_MAXREC = 20000


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
    BBox    = [[clat-dlatlon, clon-dlatlon],[clat+dlatlon,clon+dlatlon]]
    [[lllat,lllon],[urlat,urlon]] = BBox
    
    #xpos    = 100  # x-position for cross section
    #xpos    = 87  # x-position for cross section
else:
    print 'Too many standard input'
    print 'Exit'
    sys.exit()

miss = -9999.
#--------------------------------
#*********************************************
def rs2relip(a3prof, a2elev, miss_out=-9999.):
    ''' a3prof: bottom to top order in vertical axis '''
    vres = 500  # m
    tmpbottom = -1000  # m
    nl,nx,nz = a3prof.shape  # 36 layers (500m)
    tmpnz = 2*nz + 2  # 2=tmpbottom/500
    a2surfbin = ((a2elev-tmpbottom)/ vres).astype('int16')

    a3tmp = np.ones([nl,nx,tmpnz], float32) * miss_out 
    X,Y = np.meshgrid( np.arange(nx), np.arange(nl))
    #print X.shape, a3prof.shape, a2surfbin.shape
    #print ''
    a1x = X.flatten().astype('int32')
    a1y = Y.flatten().astype('int32')
    a1surfbin = a2surfbin.flatten()
    for iz in range(nz):
        a1k = a1surfbin + iz  # if iz=0 --> a1k=a1surfbin
        a3tmp[a1y,a1x,a1k] = a3prof[:,:,iz].flatten()

    return a3tmp[:,:,2:nz+2].astype('float32')
 

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

#*********************************************
def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''

    if ma.is_masked(a3in):
        a3in = a3in.filled(miss)  # 2019/12/02
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
    srcDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    a1xgmi2 = ma.masked_less(np.load(srcDir + '/Xpy.1.%06d.npy'%(oid)).flatten(), 0)
    a1ygmi2 = ma.masked_less(np.load(srcDir + '/Ypy.1.%06d.npy'%(oid)).flatten(), 0)

    a1flag  = a1xgmi2.mask + a1ygmi2.mask
    
    a1xgmi2 = a1xgmi2.filled(0) 
    a1ygmi2 = a1ygmi2.filled(0) 
    #-- Read HDF -----
    #gmiDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
    gmiDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
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

def ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale, lspecies=[0,2,3,4]):
    lsecies = [0,1,2,3,4]  #  
    nh = 28
    ny,nx = a2tIndex.shape
    a4out = empty([len(lspecies),ny, nx, nh], dtype='float32')

    for i,species in enumerate(lspecies):
        a1profScale= a3profScale[:,:,species].flatten()
        a1profNum  = a3profNum[:,:,species].flatten()
        a1tIndex   = a2tIndex.flatten()

        #-- Handle non-precipitation pixels --
        a1flag = ma.masked_equal(a1profNum, 0).mask
        a1profNum[a1flag] = 1
        a1tIndex[a1flag] = 1

        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:, a1tIndex-1, species]
        a2prof[a1flag,:] = 0.0
        a4out[i] = a2prof.reshape(ny,nx,nh)

    return a4out

def gprofLayerconversion(a3prof):
    ny,nx,nz = a3prof.shape
    a3outTop = zeros([ny,nx, (nz-20)*2],float32)
    a3outTop[:,:,0::2] = a3prof[:,:,20:]
    a3outTop[:,:,1::2] = a3prof[:,:,20:]
    return concatenate([a3prof[:,:,:20], a3outTop], axis=2)


def layer500m_to_250m(a3prof):
    ny,nx,nz = a3prof.shape
    a3out = zeros([ny,nx,nz*2], float32)
    a3out[:,:,0::2] = a3prof[:,:,:]
    a3out[:,:,1::2] = a3prof[:,:,:]
    return a3out


def load_gprof_prwat(gprofPath,miss_out=-9999.):
    with h5py.File(gprofPath,'r') as h:
        a4clusterProf= h['GprofDHeadr/clusterProfiles'][:]  # (profNumber, nlev, nT, nspecies) = (80, 28, 12, 5)
        a1hgtTopLayer= h['GprofDHeadr/hgtTopLayer'][:]  # (28,)
        species    = h['GprofDHeadr/speciesDescription'][:]
        species    = [''.join( map(chr, line) ) for line in species]

        '''
        #--- Species ---
        ['Rain Water Content\x00\x05\x00',
         'Cloud Water Content\x00!',
         'Ice Water Content\x00\x00\x00\x00',
         'Snow Water Content\x00  ',
         'Grauple/Hail Content\x00']

        #---- hgtTopLayer ----,miss_out=-9999.-
        #[ 0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5,  5. ,  5.5,
        6. ,  6.5,  7. ,  7.5,  8. ,  8.5,  9. ,  9.5, 10. , 11. , 12. ,
        13. , 14. , 15. , 16. , 17. , 18. ], dtype=float32)
        '''

        a2qFlag    = h['S1/qualityFlag'][:,83:137+1]  # (Y,X)
        a2sfcprecg = h['S1/surfacePrecipitation'][:,83:137+1]  # (Y,X)
        a2mlPrecip = h['S1/mostLikelyPrecipitation'][:,83:137+1] # (Y,X)
        a2tIndex   = h['S1/profileTemp2mIndex'][:,83:137+1] # (Y,X)  zero=missing value?
        a3profNum  = h['S1/profileNumber'][:,83:137+1,:] # (Y,X, nspecies)
        a3profScale= h['S1/profileScale'][:,83:137+1,:]  # (Y,X, nspecies)


    a4profg = ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale)
    a3profg = ma.masked_less(a4profg,0).sum(axis=0).filled(miss_out)

    a3profg = gprofLayerconversion(a3profg)

    #-- Shift to relative to elipsoid ----
    oid = int(gprofPath.split('.')[-3])
    elevDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/%04d/%02d/%02d'%(Year,Mon,Day)
    elevPath = elevDir + '/gtopo.%06d.npy'%(oid)
    a2elev   = np.load(elevPath)[:,83:137+1]
    a3profg = rs2relip(a3profg, a2elev, miss_out=-9999.)

    #a3tmp0 = np.array([[[0.5, 1, 1.5]]])
    #a2elev = np.array([[750]])
    #a3tmp1 = rs2relip(a3tmp0, a2elev, miss_out=-9999.)
    #print ''
    #print a3tmp0
    #print ''
    #print a3tmp1
    #sys.exit()


    #-------------------------------------

    a3profg = layer500m_to_250m(a3profg)
    a3profg = a3profg[:,:,::-1]   # reverse: Bottom to Top --> Top to Bottom
    return a3profg

#********************************
#- Read My data ----
#srcDir = '/home/utsumi/temp/out/my'
#outDir = '/home/utsumi/temp/out/my'
srcDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

a2prwatNS= np.load(srcDir + '/prwatprofNS.%s.npy'%(stamp))[:, xpos, :]
#a2prwatMS= a2prwatNS*0 -9999. # test

a2topprwatNS = a2prwatNS*0 - 9999.  # test
#a2topprwatNS = np.load(srcDir + '/top-prwatprofNS.%s.npy'%(stamp))[:, xpos, :]
#a2topprwatMS = a2topprwatNS*0 - 9999.  # test



print srcDir + '/prwatprofNS.%s.npy'%(stamp)
print a2prwatNS.shape
print a2prwatNS.max()

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

a2topprwatNS = a2topprwatNS[iyMy:eyMy+1,:]
#a2topprwatMS = a2topprwatMS[iyMy:eyMy+1,:]
#a2topprwatMS = a2topprwatNS*0.0 -9999.

a2prwatNS= a2prwatNS[iyMy:eyMy+1,:]
#a2prwatMS= a2prwatMS[iyMy:eyMy+1,:]
#a2prwatMS= a2prwatNS*0.0 -9999.

a1latMy  = a2latMy[iyMy:eyMy+1,:]
a1lonMy  = a2latMy[iyMy:eyMy+1,:]
print 'after extract'
print a2prwatNS
print a2prwatNS.max()
#sys.exit()
#***********************************
#- Read DPR data ----
#-- Cooresponding DPR yx --------
#srcDir = matchbaseDir + '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
#srcDir = matchbaseDir + '/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
srcDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1xdpr = np.load(xPath)[iyTmp:eyTmp+1,xpos-83]
a1ydpr = np.load(yPath)[iyTmp:eyTmp+1,xpos-83]

#-- Read DPRCMB (Ku) ----
cmbDir  = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
dprPath = glob.glob(cmbDir + '/2B.GPM.DPRGMI.*.%06d.V06A.HDF5'%(oid))[0]
with h5py.File(dprPath, 'r') as h:
    a3cmbprwatNS = h['NS/precipTotWaterCont'][:]
    a2cmbelev    = h['NS/Input/surfaceElevation'][:]
#- extract lower levels --
a3cmbprwatNS = a3cmbprwatNS[:,:,-50:]

#- Extract Radar pixels that correspond GMI ---
a2cmbprwatNS = a3cmbprwatNS[a1ydpr,a1xdpr]
a1cmbelev    = a2cmbelev[a1ydpr,a1xdpr]

#-- Read GPROF -----------
gprofDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid)
gprofPath= glob.glob(ssearch)[0]
a3prwatgprof = load_gprof_prwat(gprofPath)
print xpos, a3prwatgprof.shape
a2prwatgprof = a3prwatgprof[iyTmp:eyTmp+1,xpos-83,-50:]

#******************************************************
# Precipitation water content profile
#******************************************************
fig   = plt.figure(figsize=(12,10))
ssize = 1

nx    = a2prwatNS.shape[0]
nbin  = a2prwatNS.shape[1]
aspect1 = 0.16*nx/nbin
aspect2 = aspect1*5
#prmin, prmax = 0, 2.0
prmin, prmax = 0, 1.5
tick_locs22 = [0,8,16,24,32, 40, 48]
tick_lbls22 = [0,2,4, 6, 8, 10, 12]

for i in range(4):
    nx,nh = a2prwatNS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'

    #x0= 0.1
    #h = 0.14
    #w = 0.7

    x0= 0.1
    h = 0.2
    w = 0.7

    #--- Precip water content -------------- 
    if i==3:
        y0 = 0.02
        y0 = 0.71
        a2dat = ma.masked_less_equal(a2cmbprwatNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'CMB/NS precip water product'
        cbarlbl='g/m3'


    elif i==2:
        y0 = 0.25
        y0 = 0.48
        a2dat = ma.masked_less_equal(a2prwatNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=CMB/NS'
        cbarlbl='g/m3'


    elif i==1:
        y0 = 0.48
        y0 = 0.25
        a2dat = ma.masked_less_equal(a2topprwatNS.T,0)

        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip water(CMB/NS)'
        cbarlbl='g/m3'

    elif i==0:
        y0 = 0.71
        y0 = 0.02
        a2dat = ma.masked_less_equal(a2prwatgprof.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'GPROF precip water content'
        cbarlbl='g/m3'





    a2dat = a2dat[::-1,:] # Tot to Bottom -> Bottom to Top

    ax = fig.add_axes([x0, y0, w, h])
    cax= fig.add_axes([x0+w+0.01, y0, 0.03, h*0.9])
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')

    ax.grid()
    stime = '%04d/%02d/%02d '%(Year,Mon,Day)
    ax.set_title(stime + stype+' '+ expr)
    #plt.yticks(tick_locs, tick_lbls, fontsize=14)
    #plt.yticks([0,5],['a','b'])
    ax.set_yticks(tick_locs)
    ax.set_yticklabels(tick_lbls)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl)

    #-- Surface elevation line ------
    print a1cmbelev.shape, a2dat.shape
    x = np.arange(a2dat.shape[1])
    y = a1cmbelev / 250.  # convert to imshow axis scale
    ax.plot(x,y,'-',color='k')

##------------
outPath  = figDir + '/prof.precipTotWater.%s.%06d.y%04d-%04d.nrec%d.png'%(expr,oid,iy,ey,DB_MAXREC)
plt.savefig(outPath)
print outPath
plt.clf()





# %%
