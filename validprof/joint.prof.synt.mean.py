from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket
import numpy as np
#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir      = '/home/utsumi/temp/ret'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir      = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************
lseason=['ALL']
lregion = ['AMZ','CUS','EUS','TIB','NETP','SETP','NTP','STP','WTP','ETI','WMP','WMA','TAF','NEA']
lstop   = ['All','High','Mid','Low']
lthpr   = [0.1, 10]
iy =30 # top
ey = -30
ix =10
ex = -10
#*** 3-panels ********
for season in lseason:
    for region in lregion:
        i = -1
        ddat = {}
        for thpr in lthpr:
            for stop in lstop:
                i=i+1
                idx_db0, idx_db1 = 1, 24388
                figPath = figDir + '/prof.synt.pr%.1f.%s.%s.sh-%s.%05d-%05d.png'%(thpr, region, season, stop, idx_db0, idx_db1)
                if os.path.exists(figPath):
                    iimg = Image.open(figPath)
                    a2array = asarray(iimg)[iy:ey, ix:ex]
                else:
                    a2array = ddat[0]*0 + 255

                ddat[i] = a2array
        
        a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=1)
        a2line1 = concatenate([ddat[4],ddat[5],ddat[6],ddat[7]],axis=1)
        a2oarray = np.concatenate([a2line0, a2line1], axis=0)
        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.prof.synt.%s.%s.png'%(region,season)
        oimg.save(outPath)
        print outPath
    
 
    
