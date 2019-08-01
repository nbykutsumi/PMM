import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
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
    clat   = int(dargv['clat'])
    clon   = int(dargv['clon'])
    DB_MAXREC = int(dargv['DB_MAXREC'])
    dlatlon= int(dargv['dlatlon'])
    xpos   = int(dargv['xpos'])

elif len(argv)==1: 
    ## Africa case
    #oid = 2421
    #iy, ey = 2029, 2089
    #clat    = 14 # Africa case
    #clon    = 2  # -180 - +180
    
    ## SE.US case, oid=003556, 2014/10/14
    #oid = 3556
    #Year,Mon,Day = 2014,10,14
    ##iy, ey = 1012, 1022
    #iy, ey = 927, 1107
    #clat    = 34    # SE.US case. oid = 003556
    #clon    = -86   # 2014/10/14  05:42:03 UTC
    #DB_MAXREC = 20000

    # SW.Japan typhoon case, oid=019015, 2017/07/03
    oid = 19015
    Year,Mon,Day = 2017,7,3
    iy, ey = 1793, 1973
    clat    = 33    # SE.US case. oid = 003556
    clon    = 130   # 2014/10/14  05:42:03 UTC
    DB_MAXREC = 20000

    
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
    
    xpos    = 100  # x-position for cross section
else:
    print 'Too many standard input'
    print 'Exit'
    sys.exit()

miss = -9999.
#--------------------------------
def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''
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

        a2qFlag    = h['S1/qualityFlag'][:,83:137+1]  # (Y,X)
        a2sfcprecg = h['S1/surfacePrecipitation'][:,83:137+1]  # (Y,X)
        a2mlPrecip = h['S1/mostLikelyPrecipitation'][:,83:137+1] # (Y,X)
        a2tIndex   = h['S1/profileTemp2mIndex'][:,83:137+1] # (Y,X)  zero=missing value?
        a3profNum  = h['S1/profileNumber'][:,83:137+1,:] # (Y,X, nspecies)
        a3profScale= h['S1/profileScale'][:,83:137+1,:]  # (Y,X, nspecies)


    a4profg = ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale)
    a3profg = ma.masked_less(a4profg,0).sum(axis=0).filled(miss_out)

    a3profg = gprofLayerconversion(a3profg)
    a3profg = layer500m_to_250m(a3profg)
    a3profg = a3profg[:,:,::-1]   # reverse: Bottom to Top --> Top to Bottom
    return a3profg

#********************************
#- Read My data ----
srcDir = '/home/utsumi/temp/out/my'
outDir = '/home/utsumi/temp/out/my'
stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iy, ey, DB_MAXREC)

a2topprwatNS = np.load(srcDir + '/top-prwatprofNS.%s.npy'%(stamp))[:, xpos, :]
a2topprwatMS = a2topprwatNS*0 - 9999.  # test

a2prwatNS= np.load(srcDir + '/prwatprofNS.%s.npy'%(stamp))[:, xpos, :]
a2prwatMS= a2prwatNS*0 -9999. # test

a1latMy  = np.load(srcDir + '/lat.%s.npy'%(stamp))[:,xpos]
a1lonMy  = np.load(srcDir + '/lon.%s.npy'%(stamp))[:,xpos]

#***********************************
#- Read DPR data ----
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
#- extract lower levels --
a3cmbprwatNS = a3cmbprwatNS[:,:,-50:]

#- Extract Radar pixels that correspond GMI ---
a2cmbprwatNS = a3cmbprwatNS[a1ydpr,a1xdpr]


#-- Read GPROF -----------
gprofDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
ssearch  = gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid)
gprofPath= glob.glob(ssearch)[0]
a3prwatgprof = load_gprof_prwat(gprofPath)
print xpos, a3prwatgprof.shape
a2prwatgprof = a3prwatgprof[iy:ey+1,xpos-83,-50:]

#******************************************************
# Precipitation water content profile
#******************************************************
fig   = plt.figure(figsize=(12,12))
ssize = 1

nx    = a2prwatNS.shape[0]
nbin  = a2prwatNS.shape[1]
aspect1 = 0.16*nx/nbin
aspect2 = aspect1*5
prmin, prmax = 0, 2.0
tick_locs22 = [0,8,16,24,32, 40, 48]
tick_lbls22 = [0,2,4, 6, 8, 10, 12]

for i in range(5):
    nx,nh = a2prwatNS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'

    x0= 0.1
    h = 0.14
    w = 0.7

    #--- Precip water content -------------- 
    if i==0:
        y0 = 0.02
        a2dat = ma.masked_less_equal(a2cmbprwatNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'CMB/NS precip water product'
        cbarlbl='g/m3'

    elif i==1:
        y0 = 0.22
        a2dat = ma.masked_less_equal(a2prwatgprof.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'GPROF precip water content'
        cbarlbl='g/m3'



    elif i==2:
        y0 = 0.42
        a2dat = ma.masked_less_equal(a2prwatNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=CMB/NS'
        cbarlbl='g/m3'


    elif i==3:
        y0 = 0.62
        a2dat = ma.masked_less_equal(a2topprwatNS.T,0)

        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip water(CMB/NS)'
        cbarlbl='g/m3'

    elif i==4:
        y0 = 0.82
        a2dat = ma.masked_less_equal(a2topprwatMS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip water (CMB/MS)'
        cbarlbl='g/m3'


    a2dat = a2dat[::-1,:] # Tot to Bottom -> Bottom to Top
    ax = fig.add_axes([x0, y0, w, h])
    cax= fig.add_axes([x0+w+0.01, y0, 0.03, h*0.9])
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    ax.grid()
    ax.set_title(stype+' '+stamp)
    #plt.yticks(tick_locs, tick_lbls, fontsize=14)
    #plt.yticks([0,5],['a','b'])
    ax.set_yticks(tick_locs)
    ax.set_yticklabels(tick_lbls)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl)

##------------
outPath  = outDir + '/prof.precipTotWater.y%04d-%04d.nrec%d.png'%(iy,ey,DB_MAXREC)
plt.savefig(outPath)
print outPath
plt.clf()



'''
#******************************************************
# Precipitation profile
#******************************************************
fig   = plt.figure(figsize=(12,12))
ssize = 1

nx    = a2prNS.shape[0]
nbin  = a2prNS.shape[1]
aspect1 = 0.16*nx/nbin
aspect2 = aspect1*5
prmin, prmax = 0, 10
tick_locs22 = [22, 14, 6]
tick_lbls22 = [' 0',' 2',' 4']

for i in range(5):
    nx,nh = a2prNS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'

    x0= 0.1
    h = 0.14
    w = 0.7

    #--- Precip ---------------------------- 
    if i==0:
        y0 = 0.02
        a2dat = ma.masked_less_equal(a2dprprNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'DPR-Ku precip product'
        cbarlbl='mm/h'

    elif i==1:
        y0 = 0.22
        a2dat = ma.masked_less_equal(a2prNS.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=NS'
        cbarlbl='mm/h'

    elif i==2:
        y0 = 0.42
        a2dat = ma.masked_less_equal(a2prNScmb.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=NScmb'
        cbarlbl='mm/h'

    elif i==3:
        y0 = 0.62
        a2dat = ma.masked_less_equal(a2topprNS.T,0)*0.01

        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip (DPR/NS)'
        cbarlbl='mm/h'

    elif i==4:
        y0 = 0.82
        a2dat = ma.masked_less_equal(a2topprNScmb.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip (CMB/NS)'
        cbarlbl='mm/h'


    print ''
    print stype
    print a2dat.min(),a2dat.max()

    ax = fig.add_axes([x0, y0, w, h])
    cax= fig.add_axes([x0+w+0.01, y0, 0.03, h*0.9])
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.grid()
    ax.set_title(stype+' '+stamp)
    plt.yticks(tick_locs, tick_lbls, fontsize=14)

    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl)

##------------
outPath  = outDir + '/prof.precip.y%04d-%04d.nrec%d.png'%(iy,ey,DB_MAXREC)
plt.savefig(outPath)
print outPath
plt.clf()

'''




'''
for i in range(5):
    nx,nh = a2prNS.shape
    a1y   = arange(nh)*0.25
    a1x   = a1latMy
    cmap  = 'jet'

    #--- Precipitation retrievals ------
    if i==0:
        ax = fig.add_axes([0.05,0.05,0.15,0.22])
        #a2dat = ma.masked_less_equal(a2prNS.T,0)*0.01
        a2dat = ma.masked_less_equal(a2dprprNS.T,0)
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'DPR-Ku precip product'


    elif i==1:
        ax = fig.add_axes([0.05,0.3,0.15,0.22])
        a2dat = ma.masked_less_equal(a2prNS.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=NS'

    elif i==2:
        ax = fig.add_axes([0.05,0.6,0.15,0.22])
        a2dat = ma.masked_less_equal(a2prNScmb.T,0)*0.01
        im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=prmin, vmax=prmax)
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Ret=NScmb'


    #--- Top-Ranked Precipitation ------
    elif i==3:
        ax = fig.add_axes([0.25,0.3,0.15,0.22])
        a2dat = ma.masked_less_equal(a2topprNS.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip (DPR/NS)'

    elif i==4:
        ax = fig.add_axes([0.25,0.6,0.15,0.22])
        a2dat = ma.masked_less_equal(a2topprNScmb.T,0)*0.01
        vmin, vmax = prmin, prmax
        tick_locs, tick_lbls = tick_locs22, tick_lbls22
        aspect= aspect1
        stype = 'Top-Ranked Precip (CMB/NS)'

    #--- Zm ---------------------------- 
    elif i==5:
        ax = fig.add_axes([0.5,0.05,0.15,0.22])
        a2dat = ma.masked_less_equal(a2dprzmNS.T,0)
        vmin, vmax = zmin, zmax
        tick_locs, tick_lbls = tick_locs64, tick_lbls64
        aspect= aspect1
        stype = 'Zm DPR-Ku'

    elif i==6:
        ax = fig.add_axes([0.5,0.3,0.15,0.22])
        a2dat = ma.masked_less_equal(a2topzmNS.T,0)
        vmin, vmax = zmin, zmax
        tick_locs, tick_lbls = tick_locs64, tick_lbls64
        aspect= aspect1
        stype = 'Top-Ranked Zm (Ku)'

    elif i==7:
        ax = fig.add_axes([0.5,0.6,0.15,0.22])
        a2dat = ma.masked_less_equal(a2topzmMS.T,0)
        vmin, vmax = zmin, zmax
        tick_locs, tick_lbls = tick_locs64, tick_lbls64
        aspect= aspect1
        stype = 'Top-Ranked Zm (Ka)'

     #--- Tb ---------------------------- 
    elif i==8:
        ax = fig.add_axes([0.75,0.05,0.15,0.22])
        a2dat = ma.masked_less_equal(a2tb.T,0)
        vmin, vmax = tbmin, tbmax
        tick_locs, tick_lbls = tick_locsTb, tick_lblsTb
        aspect= aspect2
        stype = 'Obs Tb(K)'

    elif i==9:
        ax = fig.add_axes([0.75,0.3,0.15,0.22])
        a2dat = ma.masked_less_equal(a2toptbNS.T,0)
        vmin, vmax = tbmin, tbmax
        tick_locs, tick_lbls = tick_locsTb, tick_lblsTb
        aspect= aspect2
        stype = 'Top-Ranked Tb(K) Ret=NS'


    print ''
    print stype
    print a2dat.min(),a2dat.max()
    im = ax.imshow(a2dat, interpolation='none', aspect=aspect, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.grid()
    plt.yticks(tick_locs, tick_lbls, fontsize=14)
    print 'colorbar  i=',i
    plt.colorbar(im, orientation='horizontal')
    plt.title(stype+' '+stamp)

    if i in [8,9]:
        iscan1 = 0
        iscan2 = a2dat.shape[0]-1
        xoff = (iscan2 - iscan1)/14
        plt.text(xoff,12, r'10.7V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff,11, r'10.7H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff,10, r'18.7V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 9, r'18.7H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 8, r'23.8V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 7, r'36.5V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 6, r'36.5H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 5, r'89.0V'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 4, r'89.0H'        , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 3, r'166V'         , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 2, r'166H'         , family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 1, r'183.31$\pm$3V', family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')
        plt.text(xoff, 0, r'183.31$\pm$7V', family='Arial', va='center', ha='left', fontsize=10, fontweight='bold')



##------------
outPath  = outDir + '/prof.zm.tb.y%04d-%04d.nrec%d.png'%(iy,ey,DB_MAXREC)
plt.savefig(outPath)
print outPath
plt.clf()

'''
