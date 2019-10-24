import matplotlib
matplotlib.use('Agg')
import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import string
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir      = '/home.rainbow/utsumi/public_html/tempfig/validprof'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    #figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/ret'
    figDir  = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************
lseason=['JJADJF','DJF','JJA']
#lseason=['JJADJF']
#lseason = [6,7]
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

ny,nx = 120,360
nz    = 25 
latmin= -60
lonmin= -180
lrettype = ['epc','gprof']
#lrettype = ['gprof']
#lvar = ['profpmw','profrad','top-profpmw']
lvar = ['profpmw','profrad']
lstype= ['sea','land','veg','snow','coast','all']
#lstype = ['sea']
lptype= ['conv','stra','all']
#lptype= ['conv']
#lprrange=[[0.5,999],[1,3],[8,12]]
lprrange=[[1,3],[8,12]]
#lprrange=[[1,3]]
#lprrange=[[0.5,999]]

#lregion = ['TRO','SUBN','MIDN']
lregion = ['TRO','MIDN']
#lregion = ['AMZ','CUS','EUS','TIB','NETP','SETP','NTA','STA','WTP','ETI','WMP','WMA','TAF','NEA','SEC','NIN']
dBBox = {
         'TRO':  [[-15,-180],[15,180]]
        ,'SUBN': [[15,-180],[30,180]]
        ,'MIDN': [[35,-180],[50,180]]
        ,'AMZ':  [[-5,-70],[5,-53]]
        ,'CUS':  [[35,-105],[45,-95]]
        ,'EUS':  [[30,-90],[40,-80]]
        ,'TIB':  [[30,-90],[40,-80]]
        ,'NETP':  [[0,-120],[10,-110]]
        ,'SETP':  [[-10,-120],[0,-110]]
        ,'NTA' : [[0,-35],[10,-25]]
        ,'STA' : [[-10,-35],[0,-25]]
        ,'WTP':  [[0,140],[10,150]]
        ,'ETI':  [[-5,85],[5,95]]
        ,'WMP':  [[30,143],[40,158]]
        ,'WMA':  [[30,-74],[40,-59]]
        ,'TAF':  [[0, 15],[15,30]]
        ,'NEA':  [[40,120],[50,130]]
        ,'SEC':  [[22,105],[37,115]]
        ,'NIN':  [[20,75],[25,85]]
        }

#*******************************
def draw_map(a2dat=None, lregion=None, trans=False):
    fig = plt.figure(figsize=(8,4))
    #fig = plt.figure(figsize=(10,5))
    ax  = fig.add_axes([0.08,0.1,0.8,0.8])
    a1lat = np.arange(-59.5,59.5+0.01,1.0)
    a1lon = np.arange(-179.5,179.5+0.01,1.0)
    X,Y = np.meshgrid(a1lon,a1lat)
    M = Basemap(resolution='l', llcrnrlat=-60, llcrnrlon=-180, urcrnrlat=60, urcrnrlon=180, ax=ax)
    if a2dat is not None:
        im  = M.pcolormesh(X, Y, a2dat, vmin=vmin, vmax=vmax, cmap=mycm)
        cax = fig.add_axes([0.89,0.21,0.02,0.55])
        cbar= plt.colorbar(im, orientation='vertical', cax=cax)
        cax.tick_params(labelsize=13)
    
    #plt.title(stitle, fontsize=15)
    M.drawcoastlines(linewidth=1)

    M.drawmeridians(np.arange(-180,180+1,30), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50)
    M.drawparallels(np.arange(-60,60+1,30), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    #-- region boxes --
    for region in lregion:
        [[lat0,lon0],[lat1,lon1]] = dBBox[region]
        print lat0,lon0,lat1,lon1
        M.plot([lon0,lon0], [lat0,lat1], '-',c='r')
        M.plot([lon1,lon1], [lat0,lat1], '-',c='r')
        M.plot([lon0,lon1], [lat0,lat0], '-',c='r')
        M.plot([lon0,lon1], [lat1,lat1], '-',c='r')
    #------------------
    figPath = '/home/utsumi/temp/ret/map.region.png'
    plt.savefig(figPath, transparent=trans)
    print figPath
    plt.clf()


def calc_cc(x,y,axis):
    ny,nx,nz = x.shape
    xm = x.mean(axis=axis).reshape(ny,nx,1)
    ym = y.mean(axis=axis).reshape(ny,nx,1)
    A  = ((x-xm)*(y-ym)).sum(axis=axis)
    B  = ((x-xm)**2).sum(axis=axis)
    C  = ((y-ym)**2).sum(axis=axis)
    return A/( np.sqrt(B*C) )

def calc_rmse(x,y,axis):
    ny,nx,nz = x.shape
    return np.sqrt(((x-y)**2).sum(axis=axis)/nz)


def ret_lYM(seson):
    if season=='JJADJF':
        lYM = util.ret_lYM([2014,6],[2014,8]) + util.ret_lYM([2014,12],[2015,2])
    elif season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
    elif type(season)==int:
        if season <6:
            lYM = [[2015,season]]
        else:
            lYM = [[2014,season]]
    elif season =='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    else:
        print 'check season',season
        sys.exit()
    return lYM


#*******************************
#-- Elevation --------
a2orog = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#-- Freezing level ---
srcDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
d2freez = {}

dsumtmp = {}
dnumtmp = {}
for [Year,Mon] in [[2014,6],[2014,7],[2014,8],[2014,12],[2015,1],[2015,2]]:
    dsumtmp[Year,Mon]= np.load(srcDir + '/zeroDegAltituderad.sum.pr0.5.%04d%02d.sp.one.npy'%(Year,Mon))
    dnumtmp[Year,Mon]= np.load(srcDir + '/zeroDegAltituderad.num.pr0.5.%04d%02d.sp.one.npy'%(Year,Mon))

#-- JJADJF ----
lYMTmp=[[2014,6],[2014,7],[2014,8],[2014,12],[2015,1],[2015,2]]
asumtmp = np.array([dsumtmp[Year,Mon] for (Year,Mon) in lYMTmp])
anumtmp = np.array([dnumtmp[Year,Mon] for (Year,Mon) in lYMTmp])
asumtmp = ma.masked_less(asumtmp,0).sum(axis=0)
anumtmp = ma.masked_less(anumtmp,0).sum(axis=0)
d2freez['JJADJF'] = (ma.masked_where(anumtmp==0, asumtmp)/anumtmp)

#-- JJA ----
lYMTmp=[[2014,6],[2014,7],[2014,8]]
asumtmp = np.array([dsumtmp[Year,Mon] for (Year,Mon) in lYMTmp])
anumtmp = np.array([dnumtmp[Year,Mon] for (Year,Mon) in lYMTmp])
asumtmp = ma.masked_less(asumtmp,0).sum(axis=0)
anumtmp = ma.masked_less(anumtmp,0).sum(axis=0)
d2freez['JJA'] = (ma.masked_where(anumtmp==0, asumtmp)/anumtmp)

#-- DJF ----
lYMTmp=[[2014,12],[2015,1],[2015,2]]
asumtmp = np.array([dsumtmp[Year,Mon] for (Year,Mon) in lYMTmp])
anumtmp = np.array([dnumtmp[Year,Mon] for (Year,Mon) in lYMTmp])

asumtmp = ma.masked_less(asumtmp,0).sum(axis=0)
anumtmp = ma.masked_less(anumtmp,0).sum(axis=0)
d2freez['DJF'] = (ma.masked_where(anumtmp==0, asumtmp)/anumtmp)

#*** Draw region map *******************
draw_map(a2dat=None, lregion=lregion, trans=True)

#*******************************


for season in lseason:
    lYM = ret_lYM(season)

    ltype = [(stype,ptype,prrange) for stype in lstype
                                  for ptype in lptype
                                  for prrange in lprrange]
    
    for (stype,ptype,prrange) in ltype:
        thpr0,thpr1 = prrange
        dvar = {}
        dstd = {}
        dnum = {}
        dnumprof = {}
        for rettype in lrettype:
            for var in lvar:
                key = (rettype,var)

                if (rettype=='gprof')and(var=='top-profpmw'): continue
                if (rettype=='gprof')and(var=='profrad'): continue

    
                #** Initialize ******
                if nz ==1:
                    a2sum = zeros([ny,nx],float32)
                    a2num = zeros([ny,nx], int32)
                else:
                    a3sum = zeros([ny,nx,nz], float32)
                    a2num = zeros([ny,nx], int32)
                    a3ss  = zeros([ny,nx,nz], float32)
                    a3numprof = zeros([ny,nx,nz], int32)
    
                for Year,Mon in lYM:
                    if rettype =='epc': 
                        outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
        
                    elif rettype =='gprof':
                        outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
        
                    else:
                        print 'check retype',rettype
        
                    util.mk_dir(outDir)

                    #stamp  = 'pr%.1f.%04d%02d'%(thpr,Year,Mon)

                    stamp  = 's-%s.p-%s.pr-%.1f-%.1f.%04d%02d'%(stype,ptype,thpr0,thpr1,Year,Mon)

     
                    sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                    numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                    sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
                    numprofPath= outDir  + '/%s.numprof.%s.sp.one.npy' %(var,stamp)
                
                    a3sumTmp = np.load(sumPath)[:,:,:nz]
                    a2numTmp = np.load(numPath)
                    a3ssTmp  = np.load(sum2Path)
                    a3numprofTmp = np.load(numprofPath)
    
                    a3sum = a3sum + a3sumTmp
                    a2num = a2num + a2numTmp
                    a3ss  = a3ss  + a3ssTmp
                    a3numprof= a3numprof + a3numprofTmp


    
                dvar[key] = a3sum / a2num.reshape(ny,nx,1)
                dvar[key] = ma.masked_invalid(dvar[key])
    
                dstd[key] = np.sqrt( ma.masked_invalid( (a3ss - (a3sum**2)/a3numprof)/a3numprof )) 
                dnum[key] = a2num
                dnumprof[key] = a3numprof 
        #*** Draw *******************
        stampOut  = 's-%s.p-%s.pr-%.1f-%.1f.%s'%(stype,ptype,thpr0,thpr1,season)
        for region in lregion:
            [[lat0,lon0],[lat1,lon1]] = dBBox[region]
            y0 = int(floor(lat0 - latmin))
            y1 = int(floor(lat1 - latmin))
            x0 = int(floor(lon0 - lonmin))
            x1 = int(floor(lon1 - lonmin))

            #-- Freezing level ---
            #freezJJA = d2freez['JJA'][y0:y1+1,x0:x1+1]


            #-------------------
            n = dnum['epc','profrad'][y0:y1+1,x0:x1+1].sum() 
            #-- mean profile ---             
            fig = plt.figure(figsize=(2.5,5))
            ax  = fig.add_axes([0.25,0.15,0.65,0.7])
    
            a1y = 0.25 + np.arange(nz) * 0.5 # [km]
            a1rad = dvar['epc',  'profrad'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            a1pmw = dvar['epc',  'profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            a1gpr = dvar['gprof','profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            #if rettype=='epc':
            #    a1top = dvar['top-profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            
            ax.plot( a1rad, a1y, '-', c='k', linewidth=2, label='CMB') 
            ax.plot( a1pmw, a1y, '-', c='k', linewidth=1, label='EPC')
            ax.plot( a1gpr, a1y, '--', c='k', linewidth=1, label='GPROF')
            #if rettype=='epc':
            #    ax.plot( a1top, a1y, '--', c='k', linewidth=1, label='PMW 1st')
    
            stitle = '%s %s %s %s'%(region, stype, ptype, season)
            stitle = stitle + '\n'+'%.1f-%.1fmm/h N=%d'%(thpr0,thpr1, n)
            plt.title(stitle)
            plt.legend()   

            plt.xlabel('precipitation water (g/m3)')
            plt.ylabel('hight (above sea level) (km)')
            figPath = figDir + '/prof.%s.%s.png'%(stampOut,region)
            plt.savefig(figPath)
            print figPath
    
            #-- Standard deviation---             
            fig = plt.figure(figsize=(2.5,6))
            ax  = fig.add_axes([0.25,0.15,0.65,0.7])
    
            a1y = 0.25 + np.arange(nz) * 0.5 # [km]
            a1rad = dstd['epc',  'profrad'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            a1pmw = dstd['epc',  'profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            a1gpr = dstd['gprof','profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            #if rettype=='epc':
            #    a1top = dstd['top-profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            
            ax.plot( a1rad, a1y, '-', c='k', linewidth=2, label='CMB') 
            ax.plot( a1pmw, a1y, '-', c='k', linewidth=1, label='EPC')
            ax.plot( a1gpr, a1y, '--', c='k', linewidth=1, label='GPROF')
            #if rettype=='epc':
            #    ax.plot( a1top, a1y, '--', c='k', linewidth=1, label='PMW 1st')
    
            stitle = '%s %s %s %s'%(region, stype, ptype, season)
            stitle = stitle + '\n'+'%.1f-%.1fmm/h N=%d'%(thpr0,thpr1, n)

            plt.title(stitle)
            plt.legend()   

            plt.xlabel('Standard deviation (g/m3)', fontsize=12)
            plt.ylabel('hight (above sea level) (km)', fontsize=12)
            figPath = figDir + '/profstd.%s.%s.png'%(stampOut,region)
            plt.savefig(figPath)
            print figPath


            #-- count profile --- 
            fig = plt.figure(figsize=(2.5,6))
            ax  = fig.add_axes([0.25,0.15,0.65,0.7])
    
            a1y = 0.25 + np.arange(nz) * 0.5 # [km]
            a1rad = dnumprof['epc',  'profrad'][y0:y1+1,x0:x1+1].sum(axis=(0,1))
            a1pmw = dnumprof['epc',  'profpmw'][y0:y1+1,x0:x1+1].sum(axis=(0,1))
            a1gpr = dnumprof['gprof','profpmw'][y0:y1+1,x0:x1+1].sum(axis=(0,1))
            #if rettype=='epc':
            #    a1top = dstd['top-profpmw'][y0:y1+1,x0:x1+1].mean(axis=(0,1))
            
            ax.plot( a1rad, a1y, '-', c='k', linewidth=2, label='CMB') 
            ax.plot( a1pmw, a1y, '-', c='k', linewidth=1, label='PMW')
            ax.plot( a1gpr, a1y, '--', c='k', linewidth=1, label='GPROF')
            #if rettype=='epc':
            #    ax.plot( a1top, a1y, '--', c='k', linewidth=1, label='PMW 1st')
    
            stitle = '%s %s %s %s'%(region, stype, ptype, season)
            stitle = stitle + '\n'+'%.1f-%.1fmm/h N=%d'%(thpr0,thpr1, n)

            plt.title(stitle)
            plt.legend()   

            plt.xlabel('Standard deviation (g/m3)', fontsize=12)
            plt.ylabel('hight (above sea level) (km)', fontsize=12)
            figPath = figDir + '/profnum.%s.%s.png'%(stampOut,region)
            plt.savefig(figPath)
            print figPath
    
 
       
