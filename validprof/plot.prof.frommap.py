from numpy import *
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import myfunc.util as util
import socket
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir      = '/home.rainbow/utsumi/public_html/tempfig/validprof'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/validprof'

else:
    print 'check myhost'
    sys.exit()
#*******************************
ny,nx = 120,360
nz  = 25
llatlon = [[36,145],[36,-140],[36,-90],[36,-20],[36,20]]
#llatlon = [[55,30],[55,60],[55,120],[55,145],[55,-140],[55,-120],[55,-90],[55,-20]]
#llatlon = [[30+0.5,140+0.5]]
#llatlon = [[10,145],[10,-140],[10,-70],[10,-20],[10,20]]
drad = 0 # degree
lat0 = -60
lon0 = -180
lseason = ['JJA']
#lrettype= ['epc','rad','gprof','epc-top']
lrettype= ['rad','epc','gprof']
#lrettype= ['epc','rad']
#lrettype = ['epc']

#thnum = 10
thnum = 0
#-- Land fraction ----
a2landfrac = np.load(tankbaseDir + '/utsumi/validprof/const/landfrac.sp.one.120x360.npy')

#-- Elevation --------
a2orog = np.load(tankbaseDir + '/utsumi/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#---------------------
for season in lseason:

    if season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
        #lYM = util.ret_lYM([2014,6],[2014,6])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])


    for (lat,lon) in llatlon:
        y0 = int(floor((lat-drad-lat0)/1.0))
        x0 = int(floor((lon-drad-lon0)/1.0))
        if drad !=0:
            y1 = int(floor((lat+drad-lat0-0.01)/1.0))
            x1 = int(floor((lon+drad-lon0-0.01)/1.0))
        else:
            y1 = int(floor((lat+drad-lat0)/1.0))
            x1 = int(floor((lon+drad-lon0)/1.0))


        daprof = {}
        dnum   = {}
        for rettype in lrettype:
            a3sum = None
            a3num = None
            a3ss  = None 
    
            for Year,Mon in lYM:
                srcDir = tankbaseDir + '/utsumi/validprof/mapprof/%s'%(rettype)
                sumPath= srcDir + '/prof.sum.%04d%02d.sp.one.npy'%(Year,Mon)
                numPath= srcDir + '/prof.num.%04d%02d.sp.one.npy'%(Year,Mon)
                nprofPath= srcDir + '/prof.numprof.%04d%02d.sp.one.npy'%(Year,Mon)
                ssPath = srcDir + '/prof.sum2.%04d%02d.sp.one.npy'%(Year,Mon)
    
                a3sumTmp= np.load(sumPath)[:,:,:nz]
                a3ssTmp = np.load(ssPath) [:,:,:nz]
                a3nprofTmp= np.load(nprofPath)[:,:,:nz]
                a2numTmp= np.load(numPath)
    
                if a3sum is None:
                    a3sum = a3sumTmp
                    a3ss  = a3ssTmp
                    a2num = a2numTmp
                    a3nprof=a3nprofTmp
                else:
                    a3sum = a3sum + ma.masked_less(a3sumTmp,0).filled(0.0)
                    a3ss  = a3ss  + a3ssTmp
                    a2num = a2num + a2numTmp
                    a3nprof=a3nprof + a3nprofTmp
     
        
                asum = a3sum[y0:y1+1,x0:x1+1,:].sum(axis=(0,1))
                anum = a2num[y0:y1+1,x0:x1+1].sum()
                aprof =  ma.masked_invalid(asum / anum)
    
                anprof= a3nprof[y0:y1+1,x0:x1+1,:].sum(axis=(0,1))
                aprof = ma.masked_where(anprof<thnum, aprof)
                #print a3sum[y0:y1+1,x0:x1+1,:]
                #print ''
                #print a2num[y0:y1+1,x0:x1+1].shape
                #print ''

            daprof[rettype] = aprof
            dnum  [rettype] = anum

        #** Draw figure ****
        landfrac = a2landfrac[y0:y1+1,x0:x1+1].mean()
        elev     = a2orog[y0:y1+1,x0:x1+1].mean()
        print lat,lon,elev
        a1y  = arange(0,nz*0.5,0.5)       
        vmax = None
 
        fig = plt.figure(figsize=(2.5,5))
        ax  = fig.add_axes([0.12,0.1,0.8,0.7])
        
        ax.plot(daprof['rad'], a1y, '-', linewidth=2, color='k', label='DPR(%d)'%dnum['rad'])
        ax.plot(daprof['epc'], a1y, '-', linewidth=1, color='k', label='EPC(%d)'%dnum['epc'])
        ax.plot(daprof['gprof'], a1y, '--', linewidth=1, color='k', label='GPROF(%d)'%dnum['gprof'])
       
        ax.legend() 
        plt.xlim([0,vmax])
        plt.ylim([0,a1y.max()])
    
        if lat>0: slat='%dN'%(lat)
        else:     slat='%dS'%(abs(lat))
        if lon>0: slon='%dE'%(lon)
        else:     slon='%dW'%(abs(lon))
        stitle = '%s (%s,%s)'%(season, slat,slon)
        stitle = stitle + '\n' + 'LF=%.2f (%dm)'%(landfrac, elev)
        plt.title(stitle) 
        
        util.mk_dir(figDir) 
        figPath= figDir + '/prof.%s.lat%03d.lon%04d.png'%(season,lat,lon)
        plt.savefig(figPath)
        print figPath
        
