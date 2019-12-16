import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import sys, os, glob

argv = sys.argv
if len(argv)==2:
    largvs = argvs[1].split()
    dargv = {}
    for argvs in largvs:
        key,param = argvs.split('=')
        dargv[key] = param

    oid    = int(dargv['oid'])
    Year   = int(dargv['Year'])
    Mon    = int(dargv['Mon'])
    Day    = int(dargv['Day'])
    iy     = int(dargv['iy'])
    ey     = int(dargv['ey'])
    #DB_MAXREC = int(dargv['DB_MAXREC'])
    xpos   = int(dargv['xpos'])

elif len(argv)==1: 

    # check storm top case, oid=16157, 2017/1/1
    oid = 16157
    Year,Mon,Day = 2017,1,1
    #iy, ey = 1012, 1022
    #iy, ey = 1007, 1027
    #iy, ey = 962, 1072
    #iy, ey = 927, 1107
    iy, ey = None, None
    clat    = 45   #
    clon    = 110  #
    dscan   = 3

    ## Africa case
    #oid = 2421
    #iy, ey = 2029, 2089
    #clat    = 14 # Africa case
    #clon    = 2  # -180 - +180
    
    ## SE.US case, oid=003556, 2014/10/14
    #oid = 3556
    #Year,Mon,Day = 2014,10,14
    ##iy, ey = 1012, 1022
    ##iy, ey = 1007, 1027
    ##iy, ey = 962, 1072
    #iy, ey = 927, 1107
    ##clat    = 34    # SE.US case. oid = 003556
    ##clon    = -86   # 2014/10/14  05:42:03 UTC
    ##DB_MAXREC = 20000
    
    ## QJRMS case, oid=012149, 2016/4/18
    #oid = 12149
    #iy, ey = 1038, 1098
    #clat    = 32. # QJRMS case, oid=012149
    #clon    = -94 # -180 - +180
    #DB_MAXREC = 20000
    #Year,Mon,Day = 2016,4,18
    
    
    xpos    = 100  # x-position for cross section
else:
    print 'Too many standard input'
    print 'Exit'
    sys.exit()

miss = -9999.
#--------------------------------
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
        print 'No matching scan in the target domain was found.'
        print 'Exit'
        sys.exit()



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
def ret_cc(a2ref, a2dat, miss):
    a2mask1 = ma.masked_less_equal(a2ref, miss).mask
    a2mask2 = ma.masked_less_equal(a2dat, miss).mask
    a2mask  = a2mask1 + a2mask2
    a2ref = ma.masked_where(a2mask, a2ref)
    a2dat = ma.masked_where(a2mask, a2dat)

    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    return a1corr


def taylor_index(a2ref, a2dat, miss):
    a2mask1 = ma.masked_less_equal(a2ref, miss).mask
    a2mask2 = ma.masked_less_equal(a2dat, miss).mask
    a2mask  = a2mask1 + a2mask2 
    a2ref = ma.masked_where(a2mask, a2ref)  
    a2dat = ma.masked_where(a2mask, a2dat)

    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    corrmax= 1.0

    S = 4*(1.0+a1corr)**4 /((a1stdref/a1stddat + a1stddat/a1stdref)**2) / (1.0+corrmax)**4
    S = ma.masked_invalid(S).filled(miss)
    return S


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
    srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    a1xgmi2 = ma.masked_less(np.load(srcDir + '/Xpy.1.%06d.npy'%(oid)).flatten(), 0)
    a1ygmi2 = ma.masked_less(np.load(srcDir + '/Ypy.1.%06d.npy'%(oid)).flatten(), 0)

    a1flag  = a1xgmi2.mask + a1ygmi2.mask
    
    a1xgmi2 = a1xgmi2.filled(0) 
    a1ygmi2 = a1ygmi2.filled(0) 
    #-- Read HDF -----
    gmiDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)
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

        a2tIndex   = h['S1/profileTemp2mIndex'][:,83:137+1] # (Y,X)  zero=missing value?
        a3profNum  = h['S1/profileNumber'][:,83:137+1,:] # (Y,X, nspecies)
        a3profScale= h['S1/profileScale'][:,83:137+1,:]  # (Y,X, nspecies)


    a4profg = ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale)
    a3profg = ma.masked_less(a4profg,0).sum(axis=0).filled(miss_out)

    a3profg = gprofLayerconversion(a3profg)
    a3profg = layer500m_to_250m(a3profg)
    a3profg = a3profg[:,:,::-1]   # reverse: Bottom to Top --> Top to Bottom
    return a3profg

#***********************************
#-- Read GPROF -----------
gprofDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid)
gprofPath= glob.glob(ssearch)[0]

with h5py.File(gprofPath,'r') as h:
    a2latgprof = h['/S1/Latitude'][:]
    a2longprof = h['/S1/Longitude'][:]
    a2surfgprof= h['/S1/surfacePrecipitation'][:]

    a2qflag    = h['S1/qualityFlag'][:]
    a2surftype = h['S1/surfaceTypeIndex'][:]

if iy==None:
    ycnt = ret_domain_cy(a2latgprof, a2longprof, clat,clon, dlatlon=10)
    iy, ey = ycnt - dscan, ycnt + dscan

a1latgprof = a2latgprof[iy:ey+1,xpos]
a1longprof = a2longprof[iy:ey+1,xpos]
a1surfgprof= a2surfgprof[iy:ey+1,xpos]

a1qflag    = a2qflag[iy:ey+1,xpos]
a1surftype = a2surftype[iy:ey+1,xpos]

a3prwatgprof = load_gprof_prwat(gprofPath)
print xpos, a3prwatgprof.shape
a2prwatgprof = a3prwatgprof[iy:ey+1,xpos-83,-50:]  # 250m layers


#for i in range(len(a1surfgprof)):
#    print i, a1qflag[i],a1surftype[i],'  ',a2prwatgprof[i]
#
#sys.exit()
#-- Cooresponding DPR yx --------
srcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1xdpr = np.load(xPath)[iy:ey+1,xpos-83]
a1ydpr = np.load(yPath)[iy:ey+1,xpos-83]

#-- Read DPRCMB (Ku) ----
cmbDir  = '/work/hk01/PMM/NASA/GPM.DPRGMI/2B/V06/%04d/%02d/%02d'%(Year,Mon,Day)
dprPath = glob.glob(cmbDir + '/2B.GPM.DPRGMI.*.%06d.V06A.HDF5'%(oid))[0]
with h5py.File(dprPath, 'r') as h:
    a3cmbprwatNS = h['NS/precipTotWaterCont'][:]
    a2surfcmb    = h['NS/surfPrecipTotRate'][:]

#- extract lower levels --
a3cmbprwatNS = a3cmbprwatNS[:,:,-50:]

#- Extract Radar pixels that correspond GMI ---
a2cmbprwatNS = a3cmbprwatNS[a1ydpr,a1xdpr]   # 250m layers
a1surfcmb    = a2surfcmb[a1ydpr,a1xdpr]


#*** Indices *************
a2cmbprwatNSTmp = a2cmbprwatNS[:,-8:] # Exclude lowest 2km
a2cmbprwatNSTmp = a2cmbprwatNS[:,-8:] # Exclude lowest 2km
#-- Taylor's Index -------
a1taylor = taylor_index(a2cmbprwatNS, a2prwatgprof, miss=-9999.)

#-- Correlation Coefficient --
a1cc  = ret_cc(a2cmbprwatNS, a2prwatgprof, miss=-9999.)

#-- RMSE -------
a2prwatgprof = ma.masked_less(a2prwatgprof, 0)
a2cmbprwatNS = ma.masked_less(a2cmbprwatNS, 0)
a1rmse= np.sqrt( np.square( a2prwatgprof - a2cmbprwatNS ).mean(axis=1) )

#-- Abs surface bias --
a1surfbias = np.abs( ma.masked_less(a1surfgprof,0) - ma.masked_less(a1surfcmb,0))

#******************************************************
# Precipitation water content profile
#******************************************************
stamp = '%04d/%02d/%02d %05d y%04d-%04d'%(Year,Mon,Day, oid, iy, ey)

fig   = plt.figure(figsize=(12,7))
ssize = 1
nx, nh  = a2cmbprwatNS.shape
aspect1 = 0.16*nx/nh
aspect2 = aspect1*5
#prmin, prmax = 0, 2.0
prmin, prmax = 0, 100
ytick_locs22 = [0,8,16,24,32, 40, 48]
ytick_lbls22 = [0,2,4, 6, 8, 10, 12]
xtick_locs   = range(nx)[::10]
for i in range(2):
    a1y   = arange(nh)*0.25
    a1x   = a1latgprof
    cmap  = 'jet'

    x0= 0.1
    h = 0.23
    w = 0.7

    #--- Precip water content -------------- 
    if i==0:
        y0 = 0.05
        a2dat = ma.masked_less_equal(a2cmbprwatNS.T,0)
        vmin, vmax = prmin, prmax
        ytick_locs, ytick_lbls = ytick_locs22, ytick_lbls22
        aspect= aspect1
        stype = 'CMB/NS precip water product'
        cbarlbl='g/m3'

    elif i==1:
        y0 = 0.32
        a2dat = ma.masked_less_equal(a2prwatgprof.T,0)
        vmin, vmax = prmin, prmax
        ytick_locs, ytick_lbls = ytick_locs22, ytick_lbls22
        aspect= aspect1
        stype = 'GPROF precip water content'
        cbarlbl='g/m3'

    a2dat = a2dat[::-1,:] # Tot to Bottom -> Bottom to Top

    ax = fig.add_axes([x0, y0, w, h])
    cax= fig.add_axes([x0+w+0.01, y0, 0.03, h*0.9])
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    ax.grid()
    ax.set_title(stype+' '+stamp)
    #plt.yticks(tick_locs, tick_lbls, fontsize=14)
    #plt.yticks([0,5],['a','b'])
    ax.set_yticks(ytick_locs)
    ax.set_yticklabels(ytick_lbls)
    ax.set_xticks(xtick_locs)
    ax.set_xlim([-1,nx])

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl)
    print 'prof',x0, w, ax

#-- Taylor and CC ---
y0  = 0.62
ax = fig.add_axes([x0, y0, w, h*0.2])

a1dat = ma.masked_equal(a1cc,-9999.)
ax.plot(a1dat,'-',color='k')

a1dat = ma.masked_less(a1taylor,0)
ax.plot(a1dat,'--',color='red')

ax.grid()
ax.set_title('CC(black) + TaylorS(red)')
ax.set_xticks(xtick_locs)
ax.set_ylim([0,1])
ax.set_xlim([-1,nx])
print 'surf',x0, w, ax

#-- RMSE -------
y0  = 0.75
ax = fig.add_axes([x0, y0, w, h*0.2])

a1dat = ma.masked_equal(a1rmse,-9999.)
ax.plot(a1dat,'-',color='k')
ax.set_ylim([0,2])
ax.set_xlim([-1,nx])
ax.grid()

a1dat = a1surfbias
ax2 = ax.twinx()
ax2.plot(a1dat,'--',color='red')
ax2.set_ylim([0,20])
ax2.set_xlim([-1,nx])

ax.set_title('RMSE(black: left) + abs surface bias(red: right)')
ax.set_xticks(xtick_locs)
print 'RMSE',x0, w, ax
#-- Surface Precipitation ---
y0  = 0.88
ax = fig.add_axes([x0, y0, w, h*0.2])

a1dat = ma.masked_less(a1surfcmb,0)
ax.plot(a1dat,'-',color='k')

a1dat = ma.masked_less(a1surfgprof,0)
ax.plot(a1dat,'--',color='red')

ax.grid()
ax.set_title('surfac recip CMB(black) + GPROF(red)')
ax.set_xticks(xtick_locs)
ax.set_ylim([0,20])
ax.set_xlim([-1,nx])
print 'surf',x0, w, ax

##------------
outDir   = '/tank/utsumi/hometemp/validprof/caseprof'
outPath  = outDir + '/prof.prwat.%05d.y%04d-%04d.png'%(oid, iy,ey)
plt.savefig(outPath)
print outPath
plt.clf()


