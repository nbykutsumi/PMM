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
    figDir      = '/home/utsumi/temp/multi'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir      = '/home/utsumi/temp/multi'

else:
    print 'check myhost'
    sys.exit()
#*******************************
lsensor = ['AMSR2','SSMIS','ATMS','MHS']
#lsensor = ['AMSR2']
ltrange  = [(0,33), (20,20+33), (26,26+33)]
lstype = ['sea','veg']
lregion= ['tro','mid']
lptype = ['stra','conv']

lregt = [['mid',(0,7)],['mid',(20,20+7)],['tro',(26,26+7)]]

iy =1 # top
ey = -1
ix =0
ex = -25
lvar = ['ave','num']
#**********
for sensor in lsensor:
    for var in lvar:
        i = -1
        ddat = {}
        for regt in lregt:
            region,(t0,t1) = regt
            for ptype in lptype:
                for stype in lstype:
                    i=i+1
                    stamp ='%s-s-%s.p-%s.r-%s.t-%d-%d'%(sensor,stype,ptype,region,t0,t1)
                    figPath = figDir + '/prof.%s.%s.png'%(stamp,var)
                    iimg = Image.open(figPath)
                    if (i+1)%4.==0:
                        a2array = asarray(iimg)[iy:ey, ix:]
                    else:
                        a2array = asarray(iimg)[iy:ey, ix:ex]

                    ddat[i] = a2array

                    print stamp
                    print ddat[i].shape

        ddat[-1] = a2array *0 + 255


        a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=1)
        a2line1 = concatenate([ddat[4],ddat[5],ddat[6],ddat[7]],axis=1)
        a2line2 = concatenate([ddat[8],ddat[9],ddat[10],ddat[11]],axis=1)
        a2oarray = np.concatenate([a2line0, a2line1, a2line2], axis=0)

        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.prof.%s.%s.png'%(sensor,var)
        oimg.save(outPath)
        print outPath
            
 
    
