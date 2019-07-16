import matplotlib
matplotlib.use('Agg')
from numpy import *
import h5py
import myfunc.util
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import glob

#Year,Mon,Day = 2018,5,1
Year,Mon,Day = 2018,5,22
ssearch = '/work/hk01/PMM/NASA/GPM.DPR/2A/V06/%04d/%02d/%02d/2A.GPM.DPR.V*.V06A.HDF5'%(Year,Mon,Day)
lsrcPath = glob.glob(ssearch)[:1]

iy,ey = 1000,1010

for srcPath in lsrcPath:
    #srcPath = '/work/hk01/PMM/NASA/GPM.DPR/2A/V06/2018/05/01/2A.GPM.DPR.V8-20180723.20180501-S014830-E032102.023700.V06A.HDF5'
    #srcPath = '/work/hk01/PMM/NASA/GPM.DPR/2A/V06/2018/05/22/2A.GPM.DPR.V8-20180723.20180522-S235622-E012857.024041.V06A.HDF5'
    print srcPath
    with h5py.File(srcPath,'r') as h:
        a2latns = h['/NS/Latitude'][:]
        a2lonns = h['/NS/Longitude'][:]

        a2latms = h['/MS/Latitude'][:]
        a2lonms = h['/MS/Longitude'][:]
        a2laths = h['/HS/Latitude'][:]
        a2lonhs = h['/HS/Longitude'][:]
    
        a2nsurfms = h['/MS/SLV/precipRateNearSurface'][:]
        a2nsurfhs = h['/HS/SLV/precipRateNearSurface'][:]
    
    fileName = srcPath.split('/')[-1] 
    oid = fileName.split('.')[-3] 
    ver = fileName.split('.')[-2]


    a2latns = a2latns[iy:ey+1] 
    a2lonns = a2lonns[iy:ey+1] 
    a2latms = a2latms[iy:ey+1] 
    a2lonms = a2lonms[iy:ey+1] 
    a2laths = a2laths[iy:ey+1]
    a2lonhs = a2lonhs[iy:ey+1]

    #--- Figure ---- 
    stitle  = 'DPR %s %04d/%02d/%02d rev=%s'%(ver,Year,Mon,Day,oid) 

    fig = plt.figure(figsize=(10,10))
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.scatter(a2lonns, a2latns, marker='+', edgecolor='k', facecolor='k', s=80, label='NS')
    ax.scatter(a2lonms, a2latms, marker='o', edgecolor='b', facecolor='none', s=80, label='MS')
    ax.scatter(a2lonhs, a2laths, marker='D', edgecolor='r', facecolor='none', s=80, label='HS')
    ax.legend(fontsize=30)
    plt.title(stitle, fontsize=20)
    figPath = '/home/utsumi/temp/MSandHS.%s.%04d-%02d-%02d.png'%(ver,Year,Mon,Day)
    plt.savefig(figPath)
    print figPath

