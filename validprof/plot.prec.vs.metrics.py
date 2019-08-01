import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *

Year = 2017
lseason = ['ALL']
#lseason = ['DJF','MAM','JJA','SON']
lregion = ['High','Mid','Low'][::-1]
ny,nx = 120,360
srcDir = '/tank/utsumi/validprof/maperror'
#lvar = ['precbias','precbrat','profS','profrmse']
lvar = ['precrad','profS','precbias']
#llandsea = ['sea','land','coast']
llandsea = ['sea','land']

flag_byregion = True
flag_bybias   = False
#llandsea = ['land']
#a2precbias = np.load(srcDir + '/ave.precbias.%04d%02d.npy'%(Year,Mon))
#a2precbrat = np.load(srcDir + '/ave.precbrat.%04d%02d.npy'%(Year,Mon))
#a2profS    = np.load(srcDir + '/ave.profS.%04d%02d.npy'%(Year,Mon))
#a2profrmse = np.load(srcDir + '/ave.profrmse.%04d%02d.npy'%(Year,Mon))

#-- Land Sea Mask ---
a2landfrac = np.load('/tank/utsumi/validprof/const/landfrac.sp.one.120x360.npy')
a2landsea = ma.masked_less(a2landfrac,0.2).filled(0)  # Sea
a2landsea = ma.masked_inside(a2landfrac,0.2,0.8).filled(2)  # Coast
a2landsea = ma.masked_greater(a2landfrac,0.8).filled(1)  # Land

a2lon, a2lat = np.meshgrid( arange(-180+0.5,180-0.5+1,1), arange(-60+0.5,60-0.5+1,1))
da2regionmask = {}
da2regionmask['Low'] = ma.masked_outside(a2lat, -25,25).mask
da2regionmask['Mid'] = ma.masked_outside(np.abs(a2lat), 25,50).mask
da2regionmask['High'] = ma.masked_outside(np.abs(a2lat), 50,60).mask

#--------------------
for season in lseason:

    da2dat = {}
    for ivar, var in enumerate(lvar):
        lMon = util.ret_lmon(season)

        if var == 'precbias':
            a2sumrad = np.zeros([ny,nx],float32)
            a2sumpmw = np.zeros([ny,nx],float32)
            a2numrad = np.zeros([ny,nx],int32)
            a2numpmw = np.zeros([ny,nx],int32)

        else:
            a2sum = np.zeros([ny,nx],float32)
            a2num = np.zeros([ny,nx],int32)

        for Mon in lMon:
            if var=='precbias':
                a2varrad = np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precrad', Year,Mon))
                a2varpmw = np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precpmw', Year,Mon))
                a2numradTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%('precrad', Year,Mon))
                a2numpmwTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%('precpmw', Year,Mon))

                a2sumrad = a2sumrad + a2varrad * a2numradTmp
                a2sumpmw = a2sumpmw + a2varpmw * a2numpmwTmp
                a2numrad = a2numrad + a2numradTmp
                a2numpmw = a2numpmw + a2numpmwTmp

         
            else:     
                a2var = np.load(srcDir + '/ave.%s.%04d%02d.npy'%(var, Year,Mon))
                a2numTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%(var, Year,Mon))
                a2sum = a2sum + a2var * a2numTmp
                a2num = a2num + a2numTmp

       
        if var =='precbias':
            a2precrad = ma.masked_where(a2numrad==0, a2sumrad)/a2numrad
            a2precpmw = ma.masked_where(a2numpmw==0, a2sumpmw)/a2numpmw
            a2dat = a2precpmw - a2precrad
        else:
            a2dat = ma.masked_where(a2num==0, a2sum)/a2num
       
 
        da2dat[var] = a2dat

    #************************************************************
    #--- Intensity vs similarity (classified by latitude) -----
    #************************************************************
    for var in ['profS']:
        if flag_byregion==False:
            continue
        for landsea in llandsea:
            if landsea=='sea':
                a2lsmask = ma.masked_not_equal(a2landsea,0).mask
            elif landsea=='land':
                a2lsmask = ma.masked_not_equal(a2landsea,1).mask
            elif landsea=='coast':
                a2lsmask = ma.masked_not_equal(a2landsea,2).mask

            fig = plt.figure(figsize=(6,6))    
            ax  = fig.add_axes([0.15,0.15,0.75,0.75])
            for region in lregion:
                a2mask = a2lsmask + da2regionmask[region]
                a2precTmp = ma.masked_where(a2mask, da2dat['precrad'])
                a2datTmp  = ma.masked_where(a2mask, da2dat[var])

                color1 = {'Low':'salmon','Mid':'darkseagreen','High':'royalblue'}[region]
                color2 = {'Low':'brown','Mid':'darkgreen','High':'blue'}[region]
                #-- average lines ----
                a1cnt = np.arange(0.125,2.0+0.001,0.125)
                a1ave = ones(len(a1cnt))*-9999.
                for icnt, cnt in enumerate(a1cnt):
                    a2maskTmp = ma.masked_outside(a2precTmp, cnt-0.125, cnt+0.125).mask
                    a1ave[icnt] = ma.masked_where(a2maskTmp, a2datTmp).mean()
   
                a1ave = ma.masked_less(a1ave,0) 
                #---------------------


                ax.scatter(a2precTmp, a2datTmp, s=8, label=region, color=color1)
                ax.plot(a1cnt, a1ave, '-', color=color2, linewidth=4)
            ax.set_xlim([0,2])
            ax.set_ylim([0,0.7])

            ax.set_ylabel('Similarity (S-Index)',fontsize=20)
            ax.set_xlabel('Mean precipitation intensity (mm/day)',fontsize=20)
            plt.legend(fontsize=30)
            legend = ax.legend(frameon=True, markerscale=4, fontsize=20)


            plt.title(landsea, fontsize=20)
            figDir = '/tank/utsumi/hometemp/validprof'
            figPath= figDir + '/plot.byRegion.prev.vs.%s.%s.png'%(var,landsea)
            plt.savefig(figPath)
            plt.clf()
            print figPath

    #************************************************************
    #--- Intensity vs similarity (classified by bias level) -----
    #************************************************************
    for var in ['profS']:
        if flag_bybias==False:
            continue

        for landsea in llandsea:
            if landsea=='sea':
                a2lsmask = ma.masked_not_equal(a2landsea,0).mask
            elif landsea=='land':
                a2lsmask = ma.masked_not_equal(a2landsea,1).mask
            elif landsea=='coast':
                a2lsmask = ma.masked_not_equal(a2landsea,2).mask


            fig = plt.figure(figsize=(6,6))    
            ax  = fig.add_axes([0.15,0.15,0.75,0.75])
            for bias in ['Low','Mid','High']:
                #--- Bias mask -------
                a2bias = da2dat['precbias']
                if bias=='Low':
                    a2biasmask = ma.masked_outside(a2bias, -0.1, 0.1).mask
                elif bias=='Mid':
                    a2biasmask = ma.masked_outside(np.abs(a2bias),0.1,0.3).mask
                elif bias=='High':
                    a2biasmask = ma.masked_less(np.abs(a2bias),0.3).mask

                a2mask = a2lsmask + a2biasmask
                a2precTmp = ma.masked_where(a2mask, da2dat['precrad'])
                a2datTmp  = ma.masked_where(a2mask, da2dat[var])

                print bias,a2datTmp.min(), a2datTmp.max()
                color1 = {'Low':'salmon','Mid':'darkseagreen','High':'royalblue'}[bias]
                color2 = {'Low':'brown','Mid':'darkgreen','High':'blue'}[bias]
                #-- average lines ----
                a1cnt = np.arange(0.125,2.0+0.001,0.125)
                a1ave = ones(len(a1cnt))*-9999.
                for icnt, cnt in enumerate(a1cnt):
                    a2maskTmp = ma.masked_outside(a2precTmp, cnt-0.125, cnt+0.125).mask
                    print bias,icnt, 'count', (~a2maskTmp).sum()
                    if (~a2maskTmp).sum() > 20:
                        a1ave[icnt] = ma.masked_where(a2maskTmp, a2datTmp).mean()
                    else:
                        a1ave[icnt] = None
   
                a1ave = ma.masked_less(a1ave,0) 
                #---------------------


                ax.scatter(a2precTmp, a2datTmp, s=8, label=bias + ' bias', color=color1)
                ax.plot(a1cnt, a1ave, '-', color=color2, linewidth=4)
            ax.set_xlim([0,2])
            ax.set_ylim([0,0.7])

            ax.set_ylabel('Similarity (S-Index)',fontsize=20)
            ax.set_xlabel('Mean precipitation intensity (mm/day)',fontsize=20)
            plt.legend(fontsize=30)
            legend = ax.legend(frameon=True, markerscale=4, fontsize=20)


            plt.title(landsea, fontsize=20)
            figDir = '/tank/utsumi/hometemp/validprof'
            figPath= figDir + '/plot.byBias.prec.vs.%s.%s.png'%(var,landsea)
            plt.savefig(figPath)
            plt.clf()
            print figPath



