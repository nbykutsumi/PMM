import matplotlib
matplotlib.use('Agg')
from numpy import *
import socket
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import myfunc.util as util
import numpy as np
import calendar

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
lseason=['JJ']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

ny,nx = 120,360
rettype = 'epc'
lthpr = [0.5]

#-- Elevation --------
a2orog = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

a2landfrac = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/landfrac.sp.one.120x360.npy')
#*******************************
# Functions
#-------------------------------
def ret_lYM(season):
    if season=='JJA':
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
    if season=='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    return lYM
#-------------------------------


for thpr in lthpr:
    for season in lseason:
        lYM = ret_lYM(season)

        dvar = {}
        #************************
        # Storm top
        #------------------------
        for var in ['stoprad','top-stoppmw']:
            #----------------
            #** Initialize ******
            a2sum = zeros([ny,nx],float32)
            a2num = zeros([ny,nx], int32)
            for Year,Mon in lYM:
                outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
    
                stamp  = 'pr%.1f.%04d%02d'%(thpr,Year,Mon)
                sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
    
                a2sumTmp = np.load(sumPath)
                a2numTmp = np.load(numPath)
    
                a2sum = a2sum + a2sumTmp
                a2num = a2num + a2numTmp
    
    
            dvar[var] = a2sum / a2num.reshape(ny,nx)
            dvar[var] = ma.masked_invalid(dvar[var])

        #************************
        # Precip
        #------------------------
        for rettype in ['epcNScmb','dpr']:
            a2sum = zeros([ny,nx], float32)
            a2num = zeros([ny,nx], int32)
            for Year,Mon in lYM:
                if rettype in ['epcNScmb','epcNS']:
                    srcDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s.%s'%(rettype, expr)
                else:
                    srcDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s'%(rettype)
    
                iDay = 1
                #eDay = 15
                eDay = calendar.monthrange(Year,Mon)[1]
                for Day in range(iDay,eDay+1):
    
                    sumPath= srcDir  + '/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
                    numPath= srcDir  + '/prec.num.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
    
                    a2sumTmp = np.load(sumPath)
                    a2numTmp = np.load(numPath)
    
                    a2sum = a2sum + a2sumTmp
                    a2num = a2num + a2numTmp
    
            dvar[rettype] = (ma.masked_where(a2num==0, a2sum) / a2num).filled(0.0)


        #****************************
        # Scatter
        #****************************
        a2maskorog = ma.masked_greater(a2orog,2000).mask
        a2masksea  = ma.masked_less(a2landfrac,0.9).mask
        a2mask = a2maskorog + a2masksea       
 
        a2precbias = dvar['epcNScmb'] - dvar['dpr']
        a2stopbias = dvar['top-stoppmw'] - dvar['stoprad']

        a2precbias = ma.masked_where( a2mask, a2precbias)
        a2stopbias = ma.masked_where( a2mask, a2stopbias)

        plt.scatter(a2stopbias, a2precbias, s=0.3, color='k')
        #plt.imshow(a2stopbias,origin='lower')
        figPath = figDir + '/scatter.stopbias-precbias.pr%.1f.%s.png'%(thpr,season)
        plt.savefig(figPath)
        print figPath       


 
         
