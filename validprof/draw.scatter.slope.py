import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket
import pickle
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'

else:
    print 'check myhost'
    sys.exit()
#*******************************
calcflag = True
#calcflag = False
iYM  = [2014,6]
eYM  = [2014,6]
lYM  = util.ret_lYM(iYM,eYM)
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]

DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

miss_out= -9999.
lrettype = ['gprof']
#lrettype = ['epc']
#lrettype = ['epc','gprof']
#lstype   = ['all','sea','land','veg','snow','coast']
#lstype   = ['all','sea','veg']
lstype   = ['sea','veg']
#lstype   = ['sea']
#lptype   = ['conv','stra']
#lptype   = ['all','conv','stra']
lptype   = ['all']
#lph = ['A','L','H']
#lph = ['L','H']
lph = ['A']

thpr = 0
dvminmax = {'A':[-1.2,1.2], 'H':[-1.2,1.2], 'L':[-0.5,0.5]}
dvnummax = {}
for rettype in lrettype:
    for Year,Mon in lYM:
        if calcflag != True: continue
        eDay   = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,1)
        eDTime = datetime(Year,Mon,eDay)
        lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
    
        #****************************
        # Initialize
        #----------------------------    
        a1pmw  = np.array([])
        a1rad  = np.array([])
        a1stype= np.array([])
        a1ptype= np.array([])
        a1stop = np.array([])
        a1freez= np.array([])
        #----------------------------    
    
        for DTime in lDTime:
            print rettype, DTime
            Year,Mon,Day = DTime.timetuple()[:3]
            if [Year,Mon,Day] in lskipdates:
                continue
            if rettype in ['epc']:
                srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
            elif rettype=='gprof':
                srcDir  = tankbaseDir + '/utsumi/PMM/validprof/pair/gprof/%04d/%02d/%02d'%(Year,Mon,Day)
            else:
                print 'check rettype',rettype
                sys.exit()
        
            ssearch = srcDir + '/slopepmw.??????.npy'
            lprofPath = np.sort(glob.glob(ssearch))
    
            if len(lprofPath)==0:
                print 'No files'
                print ssearch
                #sys.exit() 
                continue
            for profPath in lprofPath:
                oid = int(profPath.split('.')[-2])    
        
                a1pmwTmp  = np.load(srcDir + '/slopepmw.%06d.npy'%(oid)) # g/m3/km
                a1radTmp  = np.load(srcDir + '/sloperad.%06d.npy'%(oid)) # g/m3/km
                a1stypeTmp= np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))
                a1ptypeTmp= np.load(srcDir + '/typePreciprad.%06d.npy'%(oid))
                #a1stopTmp = np.load(srcDir + '/stoprad.%06d.npy'%(oid))
                a1stopTmp = np.load(srcDir + '/stop-profrad.%06d.npy'%(oid))
                a1freezTmp= np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid))


                a1dprx = np.load(srcDir + '/dprx.%06d.npy'%(oid))
                a1precpmw = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
                a1precrad = np.load(srcDir + '/precrad.%06d.npy'%(oid))

                #---- Screening ----
                #a1flagp1 = ma.masked_greater_equal(a1precpmw, thpr).mask
                a1flagp1 = ma.masked_greater(a1precpmw, thpr).mask
                #a1flagp2 = ma.masked_greater_equal(a1precrad, thpr).mask
                a1flagp2 = ma.masked_greater(a1precrad, thpr).mask
                a1flagp  = a1flagp1 + a1flagp2

                a1flagx  = ma.masked_inside(a1dprx, 14, 34).mask  # nadir=24

                a1flag   = a1flagp * a1flagx

                if a1flag.sum()==0: continue

                a1pmwTmp   = a1pmwTmp[a1flag]
                a1radTmp   = a1radTmp[a1flag]
                a1stypeTmp = a1stypeTmp[a1flag]
                a1ptypeTmp = a1ptypeTmp[a1flag]
                a1stopTmp  = a1stopTmp[a1flag]
                a1freezTmp = a1freezTmp[a1flag]

                #-------------------
                a1pmw   = np.concatenate([a1pmw, a1pmwTmp])
                a1rad   = np.concatenate([a1rad, a1radTmp])
                a1stype = np.concatenate([a1stype, a1stypeTmp])
                a1ptype = np.concatenate([a1ptype, a1ptypeTmp])
                a1stop  = np.concatenate([a1stop,  a1stopTmp]) 
                a1freez = np.concatenate([a1freez, a1freezTmp]) 


        #-- Surface type flag -------------------
        a1flagsea  = ma.masked_equal(a1stype,1).mask
        a1flagveg  = ma.masked_inside(a1stype,3,7).mask
        a1flagsnow = ma.masked_inside(a1stype,8,11).mask
        a1flagcoast= ma.masked_equal(a1stype,13).mask
        a1flagland = a1flagveg + a1flagsnow

        #-- Precip type ----------
        a1ptype = ma.masked_less(a1ptype,0).astype('int16')
        a1stra  = (a1ptype %10).astype('int16')
        a1conv  = (a1ptype%100-a1stra).astype('int16')/10
        a1other = (a1ptype/100).astype('int16')
        a1all   = (a1stra + a1conv + a1other).astype(float32)
        a1stra  = a1stra / a1all
        a1conv  = a1conv / a1all
        #-- Precipitation type flag -------------
        a1flagconv = ma.masked_greater(a1conv,0.6).mask
        a1flagstra = ma.masked_greater(a1stra,0.6).mask


        #-- Histograms log-scale -------
        lkey = [(stype,ptype,ph)
                for stype in lstype
                for ptype in lptype
                for ph    in lph
                ]

        for key in lkey:
            (stype,ptype,ph) = key

            #-- Surface type ---
            if stype == 'sea':
                a1flagstype = a1flagsea
            elif stype=='veg':
                a1flagstype = a1flagveg
            elif stype=='snow':
                a1flagstype = a1flagsnow
            elif stype=='coast':
                a1flagstype = a1flagcoast
            elif stype=='land':
                a1flagstype = a1flagland
            elif stype=='all':
                a1flagstype = np.array([True]*len(a1stype))
            else:
                print 'check stype',stype
                sys.exit()
 
            if type(a1flagstype) is np.bool_:
                a1flagstype = np.array([a1flagstype]*len(a1stype))
        
            #-- Precipitation type --
            if ptype =='conv':
                a1flagptype = a1flagconv
            elif ptype=='stra':
                a1flagptype = a1flagstra
            elif ptype=='all':
                a1flagptype = np.array([True])
            else:
                print 'check ptype',ptype
                sys.exit()

            #-- Precipitation height ------------
            a1fpd = ma.masked_less(a1stop,0) - ma.masked_less(a1freez,0)  # freezing precipitation depth
            if ph =='L':
                a1flagph = ma.masked_less(a1fpd, -0.5).mask  # 500m lower than FL
            elif ph=='H':
                a1flagph = ma.masked_greater(a1fpd, 0.5).mask # 500m higher than FL
            elif ph=='A':
                a1flagph = np.array([True]*len(a1fpd))

            else:
                print 'check ph',ph
                sys.exit()


            #-- Screen ---------------
            vmin, vmax = dvminmax[ph]

            a1flag = a1flagstype * a1flagptype * a1flagph

            a1pmwTmp = a1pmw[a1flag]
            a1radTmp = a1rad[a1flag]

            bins  = np.arange(vmin,vmax+0.0001, 0.025)
            H,xedges,yedges = np.histogram2d(a1radTmp, a1pmwTmp, bins = bins)
            H = H.T
        
            #*******************************
            # Save
            #-------------------------------
            stamp =  '%s.s-%s.p-%s.ph-%s.%04d.%02d'%(rettype,stype, ptype, ph, Year, Mon)
            
            pickleDir  = '/home/utsumi/temp/ret/pickle'
            util.mk_dir(pickleDir)
            
            histoPath  = pickleDir + '/histo.slope.%s.bfile'%(stamp)
            binsPath   = pickleDir + '/bins.slope.%s.npy'%(stamp)
            if calcflag == True:
                with open(histoPath, 'wb') as f:
                    pickle.dump(H, f)
            
                np.save(binsPath, bins)

                print histoPath




    #*******************************
    # Load
    #-------------------------------
    lkey = [(stype,ptype,ph)
            for stype in lstype
            for ptype in lptype
            for ph    in lph
            ]

    for key in lkey:
        (stype,ptype,ph) = key

        vmin, vmax = dvminmax[ph]

        H = None
        for Year,Mon in lYM:
            pickleDir  = '/home/utsumi/temp/ret/pickle'

            stamp =  '%s.s-%s.p-%s.ph-%s.%04d.%02d'%(rettype,stype, ptype, ph, Year, Mon)

            histoPath  = pickleDir + '/histo.slope.%s.bfile'%(stamp)
            binsPath   = pickleDir + '/bins.slope.%s.npy'%(stamp)

            with open(histoPath, 'r') as f:
                Htmp = pickle.load(f)
            bins = np.load(binsPath)

            if H is None:
                H = Htmp
            else:
                H = H + Htmp

        ny,nx = H.shape
        nxh = int(nx/2)
        nyh = int(ny/2)

        nall= float(H.sum())

        r1 = H[nyh:, nxh:].sum() / nall
        r2 = H[nyh:, :nxh].sum() / nall
        r3 = H[:nyh, :nxh].sum() / nall
        r4 = H[:nyh, nxh:].sum() / nall
        
        #-- Figure density plot ----
        stampout = '%s.s-%s.p-%s.ph-%s'%(rettype,stype,ptype,ph)
        X,Y = np.meshgrid(bins,bins)

        if rettype==lrettype[0]:
            dvnummax[key] = np.percentile(H,99.9)

        fig = plt.figure(figsize=[6,6])
        ax  = fig.add_axes([0.15,0.13,0.68,0.68])
        im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=dvnummax[key])

        #-- plot 1:1 line
        #ax.plot(array([vmin,vmax]),array([vmin,vmax]),'-',color='k',linewidth=0.5)

        #-- plot x=0 and y=0 lines --
        ax.plot(array([vmin,vmax]),array([0,0]),'-',color='k',linewidth=1)
        ax.plot(array([0,0]),array([vmin,vmax]),'-',color='k',linewidth=1)


        #-- axis labels ---
        if ph in ['L']:
            lticks = [-0.5, 0, 0.5]
        else:
            lticks = [-1,0,1]
         
        lticklabels = lticks

        ax.set_xticks(lticks)
        ax.set_xticklabels(lticklabels, fontsize=16)
        ax.set_yticks(lticks)
        ax.set_yticklabels(lticklabels, fontsize=16)

        ax.set_xlabel('CMB [g/m3/km]', fontsize=22)
        ax.set_ylabel('%s [g/m3/km]'%(rettype), fontsize=22)
        ax.set_ylim([vmin,vmax])
        ax.set_xlim([vmin,vmax])

        ntot = H.sum()
        rettypeout= str.upper(rettype)
        if stype=='sea':
            stypeout = 'Ocean'
        else:
            stypeout = str.capitalize(stype)
        ptypeout = str.capitalize(ptype)
        phout = {'A':'All', 'L':'Warm', 'H':'Cold'}[ph]
        #plt.title('%s'%(stampout), fontsize=18)
        plt.title('%s %s %s %s N=%d'%(rettypeout, ptypeout, stypeout, phout, ntot), fontsize=18)

        #-- Text for Fractions ----
        tsize = 30
        mycm = 'r'
        plt.text( 0.7, 0.9,  '%.2f'%(r1), transform=ax.transAxes, fontsize=tsize, color=mycm)
        plt.text( 0.1, 0.9,  '%.2f'%(r2), transform=ax.transAxes, fontsize=tsize, color=mycm)
        plt.text( 0.1, 0.03, '%.2f'%(r3), transform=ax.transAxes, fontsize=tsize, color=mycm)
        plt.text( 0.7, 0.03, '%.2f'%(r4), transform=ax.transAxes, fontsize=tsize, color=mycm)

        cax = fig.add_axes([0.84,0.15,0.02, 0.6])
        cbar=plt.colorbar(im, orientation='vertical', cax=cax)
        cbar.ax.tick_params(labelsize=16)

        #figPath= figDir + '/scatter.%s.%s.%s.%s.%s.png'%(exprTmp,region,obstype,surftype,season)
        figDir = '/home/utsumi/temp/ret'
        figPath= figDir + '/scatter.slope.%s.png'%(stampout)
        util.mk_dir(figDir)

        plt.savefig(figPath)
        print figPath

