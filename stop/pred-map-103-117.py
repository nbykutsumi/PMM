import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
#import matplotlib.gridspec as gridspec
from glob import glob
import tensorflow as tf
import numpy.ma as ma
import sys,os, glob
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
from mpl_toolkits.basemap import Basemap
import socket
import h5py

myhost = socket.gethostname()
if myhost=='DESKTOP-RVR6SB8':
    pmmbaseDir ='/mnt/j/PMM'
    hdfbaseDir =''
    figDir = '/mnt/c/ubuntu/fig'
elif myhost =='shui':
    pmmbaseDir ='/tank/utsumi/PMM'
    hdfbaseDir ='/work/hk01/PMM/NASA'
    figDir = '/home.rainbow/utsumi/public_html/tempfig/stop'
#elif myhost =='well':
#    pmmbaseDir ='/tank/utsumi/PMM'
#    figDir = '/mnt/c/ubuntu/fig'

else:
    print 'check myhost',myhost
#get_ipython().magic(u'matplotlib inline')
#******************************************
# Functions
#******************************************
def shift_array(ain=None, dy=None,dx=None,miss=-9999):
    ny,nx,nz = ain.shape
    aout = np.ones([ny,nx,nz]).astype(ain.dtype)*miss
    if   dy<=0: iy0=0; ey0=ny-abs(dy); iy1=abs(dy); ey1=ny
    elif dy> 0: iy0=abs(dy); ey0=ny; iy1=0; ey1=ny-abs(dy)
    if   dx<=0: ix0=0; ex0=nx-abs(dx); ix1=abs(dx); ex1=nx
    elif dx> 0: ix0=abs(dx); ex0=nx; ix1=0; ex1=nx-abs(dx)

    aout[iy0:ey0,ix0:ex0] = ain[iy1:ey1,ix1:ex1]
    return aout


def read_pc_coef(isurf):
    #*********************************
    # Read PC coefficient
    #*********************************
    #coefDir = '/work/hk01/utsumi/PMM/stop/data/coef'
    #coefDir = '/mnt/j/PMM/stop/data/coef'
    coefDir = pmmbaseDir + '/stop/data/coef'
    egvecPath = coefDir + '/egvec.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    egvalPath = coefDir + '/egval.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    varratioPath = coefDir + '/varratio.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    
    a2egvec = np.load(egvecPath)  # (n-th, ncomb)
    a1varratio = np.load(varratioPath)
    a1cumvarratio= np.cumsum(a1varratio)
    return a2egvec, a1varratio, a1cumvarratio

def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)

def unit(x):
    return ( x - np.min(x,0) )/( np.max(x,0) - np.min(x,0) )


def rmse(x,y):
    x = x.flatten()
    y = y.flatten()
    return np.sqrt((((x-y))**2).mean())
def Rmse(x,y):
    Min,Max=MinStop,MaxStop
    return np.sqrt( ( ( ((Max-Min)*x+Min).flatten()-((Max-Min)*y+Min).flatten() )**2 ).mean() )

def cc(x,y):
    return np.corrcoef( x.flatten(), y.flatten() )[0,1]


print 'Define functions'
#******************************************
# Main loop start
#******************************************
lorbit = [[2017,1,20,16452]]
restriction = 10
coef_b  = 5
ldy   = [-1,0,1]
ldx   = [-3,-2,-1,0,1,2,3]
#ldx   = [-2,-1,0,1,2]
ldydx = [[dy,dx] for dy in ldy for dx in ldx]
ntc1  = 9
ntc2  = 4
ncomb = (ntc1+ntc2)* len(ldydx)
imid  = int((len(ldy)*len(ldx)-1)/2)
lisurf = range(1,14)
#norb  =30
norb  =10
#lact = ['H','L','LT']
#act = 'LT'  # 'H', 'L', 'LT'
act = 'L'  # 'H', 'L', 'LT'

latmin,latmax = -90, 90
#BBox = [[28,-105],[46,-76]]
BBox = [[30,-90],[43,-76]]
[[lllat,lllon],[urlat,urlon]] = BBox

for (yyyy,mm,dd,oid) in lorbit:
    print oid,yyyy,mm,dd
    tcPath1 = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d/Tc.%06d.npy'%(yyyy,mm,dd,oid)
    tcPath2 = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d/TcS2.1.%06d.npy'%(yyyy,mm,dd,oid)
    
    stopPath= pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.Ku.V06A.heightStormTop/%04d/%02d/%02d/heightStormTop.1.%06d.npy'%(yyyy,mm,dd,oid)
    surfPath= pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.GMI.surfaceTypeIndex/%04d/%02d/%02d/surfaceTypeIndex.%06d.npy'%(yyyy,mm,dd,oid)
    latPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.GMI.Latitude/%04d/%02d/%02d/Latitude.%06d.npy'%(yyyy,mm,dd,oid)
    lonPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp103-117.GMI.Longitude/%04d/%02d/%02d/Longitude.%06d.npy'%(yyyy,mm,dd,oid)
    t2mPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(yyyy,mm,dd,oid)
    
    
    if 'H' in act:
        atc1 = np.load(tcPath1)
        atc2 = np.load(tcPath2)
        atc  = np.concatenate([atc1,atc2],axis=2)
    elif 'L' in act:
        atc  = np.load(tcPath1)
    
    if 'T' in act:
        at2mOrg = np.load(t2mPath)[:,103:117+1]
    
    a2stop   = np.load(stopPath)
    a2surf   = np.load(surfPath)
    a2lat    = np.load(latPath)
    a2lon    = np.load(lonPath)
    
    
    aTcOrg  = np.concatenate([ shift_array(atc, dy, dx) for dy in ldy for dx in ldx], axis=2)
    

    print 'predTc.shape',aTcOrg.shape
    nz = aTcOrg.shape[2]
    aTcOrg     = aTcOrg.reshape(-1, nz)
    asurfOrg   = a2surf.flatten()
    astopOrg   = a2stop.flatten()
    alatOrg    = a2lat.flatten()
    alonOrg    = a2lon.flatten()

    if 'T' in act:
        aTcOrg =   np.concatenate([aTcOrg, at2mOrg.reshape(-1,1)], axis=1)
    print aTcOrg.shape, alatOrg.shape, asurfOrg.shape,astopOrg.shape
    #******************************************
    ny,nx = atc.shape[:2]
    predStop = np.ones(ny*nx)*(-9999.)  # 1-D

    #******************************************
    # Screen surf type
    #******************************************
    for isurf in lisurf:
        a1flagsurf  = ma.masked_equal(asurfOrg, isurf).mask

        #******************************************
        # Screen invalid Tc
        #******************************************
        a1flagtc  = ma.masked_inside(aTcOrg, 50, 350).all(axis=1).mask

        #******************************************
        # Extract
        #******************************************
        a1flag = a1flagsurf * a1flagtc
        a1idx = np.arange(aTcOrg.shape[0]).astype('int32')
        a1idx = a1idx[a1flag]

         
        aTc   = aTcOrg[a1idx]
        astop = astopOrg[a1idx]
        
        #******************************************
        # Preprocess parameters 
        #******************************************
        preptype = 'act.%s.nynx.%dx%d.isurf.%d.Mon.%d.Lat.%d.%d'%(act, len(ldy),len(ldx),isurf,mm,latmin,latmax)
        paramDir = pmmbaseDir + '/stop/prep-param/%s'%(preptype)
     
        #*** Tc mean, std *******
        meanPath = paramDir + '/mean.Tc.npy'
        stdPath  = paramDir + '/std.Tc.npy'
        ameanTc = np.load(meanPath)
        astdTc  = np.load(stdPath)
        
        #*** PC coefficient (eigen vector)
        egvecPath = paramDir + '/egvec.npy'
        varratioPath = paramDir + '/varratio.npy'
        a2egvec   = np.load(egvecPath)
        a1varratio= np.load(varratioPath)
        
        #*** PC min, max ***************
        minPath = paramDir + '/pc.min.npy'
        maxPath = paramDir + '/pc.max.npy'
        MinPC   = np.load(minPath)
        MaxPC   = np.load(maxPath)
        
        #******************************************
        # PCA
        #******************************************
        aTc = (aTc-ameanTc)/astdTc
        #reduction = np.dot(predTc, a2egvec[:restriction,:].T)
        reduction = np.dot(aTc, a2egvec[:restriction,:].T)
        
        #******************************************
        # Normalize by min and max
        #******************************************
        MinStop = 0
        MaxStop = 32000
        
        aX = my_unit(reduction,MinPC,MaxPC)
        aY = my_unit(astop,MinStop,MaxStop)

        #******************************************
        # Checkpoint
        #******************************************
        expr= 'act.%s.surf%d.b%d.lat.%d.%d'%(act,isurf,coef_b,latmin,latmax)
        #ckptDir = '/work/hk01/utsumi/PMM/stop/ml-param-%d'%(act)
        ckptDir = pmmbaseDir + '/stop/ml-param/%s'%(expr)
        ckptPath= ckptDir + '/stop.%02d'%(isurf)
        #******************************************
        # Prediction
        #******************************************
        print aX.shape
        with tf.Session() as sess:
            saver = tf.train.import_meta_graph(ckptPath + '.meta')
            saver.restore(sess, tf.train.latest_checkpoint(ckptDir + '/'))
            graph = tf.get_default_graph()
            X      = graph.get_tensor_by_name('input:0')
            pred   = graph.get_tensor_by_name('pred:0')
            print X
            print pred
            out    = sess.run(pred, feed_dict={X:aX}).flatten()
        print out.shape

        prediction = (MaxStop - MinStop)*out + MinStop

        print a1idx.shape, prediction.shape
        print 'prediction.min(),max()',prediction.min(),prediction.max()
        
         
        predStop[a1idx] = prediction
        
       
    #******************************************
    # Reshape
    #******************************************
    predStop = predStop.reshape(ny,nx)

    #******************************************
    # Mask radar stop <=0
    #******************************************
    predStop = ma.masked_where(a2stop<=0, predStop)

    #******************************************
    # Mask predicted stop <=0
    #******************************************
    predStop = ma.masked_less_equal(predStop, 0)
    predStop = predStop / 1000. 

    a2stop   = ma.masked_less_equal(a2stop,0)/1000.
    #******************************************
    # Figure
    #******************************************
    fig = plt.figure(figsize=(10,5))
    vmin = 2
    vmax = 8
    #mycm = 'gist_ncar'
    mycm = 'jet'
    meridians = np.arange(-180,180, 2)
    parallels = np.arange(-90,90, 2)

    #******************************************
    # Prediction
    #******************************************
    ax  = fig.add_axes([0.05,0.15,0.4,0.75])
    M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

    im  = M.pcolormesh(a2lon, a2lat, predStop, vmin=vmin, vmax=vmax, cmap=mycm)
    M.drawcoastlines()

    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')

    #******************************************
    # Radar
    #******************************************
    ax  = fig.add_axes([0.55,0.15,0.4,0.75])
  
    M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

    im  = M.pcolormesh(a2lon, a2lat, a2stop, vmin=vmin, vmax=vmax, cmap=mycm)
    M.drawcoastlines()

    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d')
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')

    ##******************************************
    ## Boundary lines for radar
    ##******************************************
    #hdfDir = hdfbaseDir + '/GPM.Ku/2A/V06/%04d/%02d/%02d'%(yyyy,mm,dd)
    #hdfPath= glob.glob(hdfDir + '/2A.GPM.Ku.*.%06d.????.HDF5'%(oid))[0]

    #with h5py.File(hdfPath,'r') as h:
    #    a2lathdf = h['/NS/Latitude'][:]
    #    a2lonhdf = h['/NS/Longitude'][:]

    #a1lathdf0 = a2lathdf[:,0]
    #a1lonhdf0 = a2lonhdf[:,0]
    #a1lathdf1 = a2lathdf[:,1]
    #a1lonhdf1 = a2lonhdf[:,1]
    #M.plot(a1lonhdf0, a1lathdf0,color='k')
    #M.plot(a1lonhdf1, a1lathdf1,color='k')

    #******************************************
    # Colorbar
    #******************************************
    ax = fig.add_axes([0.1,0.07,0.8,0.04])
    plt.colorbar(im, orientation='horizontal', cax=ax) 
    #******************************************
    # Save
    #******************************************
    figPath = figDir + '/mnt/c/ubuntu/fig/fig.stop.map.png'
    figPath = figDir + '/fig.stop.map.png'
    pl.savefig(figPath)
    print figPath 
    
 
