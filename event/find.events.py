import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
import glob
from datetime import datetime, timedelta
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys, os

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,1,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

wy    = 10
wx    = 21
wyh   = int(wy/2)
wxh   = int(wx/2)
cx    = 110 # GMI center angle bin (0,1,2,..)
cx2   = int((137-83+1)/2)

gmibaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
dprbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.9ave.precipRateESurface'
dprorgbaseDir='/work/hk01/PMM/NASA/GPM.Ku/2A/V06'

def draw_lines(ax):
    # lines
    ax.plot(a1edgelon_gmi0, a1edgelat_gmi0, '-', color='k')
    ax.plot(a1edgelon_gmi1, a1edgelat_gmi1, '-', color='k')

    ax.plot(a1edgelon_dpr0, a1edgelat_dpr0, '-', color='k')
    ax.plot(a1edgelon_dpr1, a1edgelat_dpr1, '-', color='k')

    ax.plot(a1lonBBox0, a1latBBox0, '-', color='k')
    ax.plot(a1lonBBox1, a1latBBox1, '-', color='k')
    ax.plot(a1lonBBox2, a1latBBox2, '-', color='k')
    ax.plot(a1lonBBox3, a1latBBox3, '-', color='k')

    # Grid lines
    M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.7, rotation=60)
    M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.7)

    # Coastlines
    M.drawcoastlines()

    return ax


for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    gmiDir = gmibaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    dprDir = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    dprorgDir=dprorgbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)

    
    ssearch = gmiDir + '/2A.GPM.GMI*.HDF5'
    lgmiPath = sort(glob.glob(ssearch))

    for gmiPath in lgmiPath:
        gnum = gmiPath.split('.')[-3]
        print gmiPath
        print gnum

        #- Read GMI precipitation --
        h = h5py.File(gmiPath,'r')
        a2gmiFull    = h['S1/surfacePrecipitation']
        a2latgmiFull = h['S1/Latitude']
        a2longmiFull = h['S1/Longitude']

        a2gmi    = a2gmiFull[:,83:137+1]
        a2latgmi = a2latgmiFull[:,83:137+1]
        a2longmi = a2longmiFull[:,83:137+1]

        #- Read GMI time info --
        grp  = 'S1/ScanTime'
        a1yyyy = h[grp+'/Year' ][:]
        a1mm   = h[grp+'/Month'][:]
        a1dd   = h[grp+'/DayOfMonth'][:]
        a1hh   = h[grp+'/Hour'][:]
        a1mn   = h[grp+'/Minute'][:]
        #ss   = h[grp+'/Second'][:]


        #- Read DPR precipitation (9-grid average) --
        dprPath = dprDir + '/precipRateESurface.%s.npy'%(gnum)
        a2dpr   = np.load(dprPath)
        ny,nx   = a2gmi.shape


        #- Read DPR original (only for figure) --
        ssearch    = dprorgDir + '/2A.GPM.Ku.*.%s.V06A.HDF5'%(gnum)
        dprorgPath = glob.glob(ssearch)[0]
        h = h5py.File(dprorgPath,'r')
        a2dprorg   = h['NS/SLV/precipRateESurface'][:]
        a2latdprorg   = h['/NS/Latitude'][:]
        a2londprorg   = h['/NS/Longitude'][:]

        #- Search events -

        for iy0 in range(0,ny-ny%wy-wy, wyh):
            dpr = ma.masked_less(a2dpr[iy0:iy0+wy, cx2-wxh:cx2+wxh], 0).mean()
            gmi = ma.masked_less(a2gmi[iy0:iy0+wy, cx2-wxh:cx2+wxh], 0).mean()
            if (type(dpr)==np.ma.core.MaskedConstant)or(type(gmi)==np.ma.core.MaskedConstant):
                continue

            #thpr = 3
            thpr = 5
            if ((dpr<thpr)and(gmi<thpr)):
                continue


            prcp_small = min(dpr,gmi)
            prcp_large = max(dpr,gmi)
            if prcp_small==0:
                rat ==9999
            else:
                rat = prcp_large/prcp_small

            if rat<2: continue


            #adpr = ma.masked_less(a2dpr[iy0:iy0+wy, cx2-wxh:cx2+wxh], 0)
            #agmi = ma.masked_less(a2gmi[iy0:iy0+wy, cx2-wxh:cx2+wxh], 0)
            #sys.exit()

            print '%.1f, %.1f, %.1f'%(gmi,dpr,rat)
            y0 = int(iy0+wy*0.5)
            x0 = int(cx+wx*0.5)
            clat = a2latgmiFull[y0,x0]
            clon = a2longmiFull[y0,x0]

            #*******************************
            #-- Draw snapshot --
            #*******************************
            fig = plt.figure(figsize=(15,15))
            dlat,dlon = 5,5
            #dlat,dlon = 1.5, 1.5
            lllat = clat -dlat
            urlat = clat +dlat
            lllon = clon -dlon
            urlon = clon +dlon

            dgrid = 1.0
            meridians = arange(int(lllon)-3,int(urlon)+3,dgrid)
            parallels = arange(int(lllat)-3,int(urlat)+3,dgrid)

            gmisize = 5

            vmin_precip = 0
            vmax_precip = 10
            #cmap_precip = 'gist_ncar'
            cmap_precip = 'nipy_spectral'
            #-- Edge lines --
            a1edgelat_gmi0 = a2latgmiFull[:,15]
            a1edgelat_gmi1 = a2latgmiFull[:,-15]
            a1edgelon_gmi0 = a2longmiFull[:,16]
            a1edgelon_gmi1 = a2longmiFull[:,-16]

            a1edgelat_dpr0 = a2latdprorg[:,0]
            a1edgelat_dpr1 = a2latdprorg[:,-1]
            a1edgelon_dpr0 = a2londprorg[:,0]
            a1edgelon_dpr1 = a2londprorg[:,-1]

            #-- BBox lines --
            #a1latBBox0 = a1latgmi[iy0:iy0+wy, cx2-wxh:cx+wxh]
            a1latBBox0 = a2latgmi[iy0:iy0+wy, cx2-wxh]
            a1latBBox1 = a2latgmi[iy0:iy0+wy, cx2+wxh-1]
            a1latBBox2 = a2latgmi[iy0,      cx2-wxh:cx2+wxh]
            a1latBBox3 = a2latgmi[iy0+wy-1, cx2-wxh:cx2+wxh]

            a1lonBBox0 = a2longmi[iy0:iy0+wy, cx2-wxh]
            a1lonBBox1 = a2longmi[iy0:iy0+wy, cx2+wxh-1]
            a1lonBBox2 = a2longmi[iy0,      cx2-wxh:cx2+wxh]
            a1lonBBox3 = a2longmi[iy0+wy-1, cx2-wxh:cx2+wxh]

            #-- Date --
            yyyy  = a1yyyy[iy0]
            mm    = a1mm[iy0]
            dd    = a1dd[iy0]
            hh    = a1hh[iy0]
            mn    = a1mn[iy0]
            sdate = '%04d-%02d-%02d UTC%02d:%02d %s'%(yyyy,mm,dd,hh,mn,gnum)

            #*****************
            #-- GMI -------
            a2fig = ma.masked_less_equal(a2gmiFull,0)
            a2lon = a2longmiFull
            a2lat = a2latgmiFull
            ax    = fig.add_axes([0.1,0.6,0.3,0.3])
            cax   = fig.add_axes([0.4+0.01,0.6,0.02,0.3])

            ax.set_title('GMI'+' '+sdate)

            M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

            im  = ax.scatter(a2lon, a2lat, c=a2fig, marker='o',s=gmisize, vmin=vmin_precip, vmax=vmax_precip, cmap=cmap_precip)

            #-- Borders, Grids, Coastlines, colorbars
            ax = draw_lines(ax)
            plt.colorbar(im, cax=cax)


            #*****************
            ##-- DPR-averaged --
            a2fig = ma.masked_less_equal(a2dpr,0)
            a2lon = a2longmi
            a2lat = a2latgmi
            ax    = fig.add_axes([0.5,0.6,0.3,0.3])
            cax   = fig.add_axes([0.8+0.01,0.6,0.02,0.3])
            ax.set_title('DPR-Ave'+' '+sdate)

            M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
            im  = ax.scatter(a2lon, a2lat, c=a2fig, marker='o',s=gmisize, vmin=vmin_precip, vmax=vmax_precip, cmap=cmap_precip)

            #-- Borders, Grids, Coastlines, colorbars
            ax = draw_lines(ax)
            plt.colorbar(im, cax=cax)

            #-- Save --
            figDir  = '/work/hk01/utsumi/PMM/event/fig/th%02dmmh/%04d/%02d'%(thpr,Year,Mon)
            figPath = figDir + '/snap.%04d.%02d.%02d.%s.Lat%04.1f.Lon%05.1f.png'%(Year,Mon,Day,gnum,clat,clon)
            util.mk_dir(figDir)
            plt.savefig(figPath)
            print figPath
            sys.exit()



             
