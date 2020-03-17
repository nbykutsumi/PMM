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
lregfl = [['mid',(0,1)],['mid',(3,4)],['tro',(4,5)]]
lstype = ['sea','veg']
lptype= ['stra','conv']

iy =10 # top
ey = -1
ix =20
ex = -15
lvar = ['ave','num']
rel = 'SF'
#**********
for var in lvar:
    i = -1
    ddat = {}
    for regfl in lregfl:
        region,(fz0,fz1) = regfl
        for ptype in lptype:
            for stype in lstype:
                i=i+1
                stamp ='rel%s.s-%s.p-%s.r-%s.fz-%d-%d'%(rel,stype,ptype,region,fz0,fz1)
                figPath = figDir + '/prof.%s.%s.png'%(stamp,var)
                iimg = Image.open(figPath)
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
    outPath = figDir + '/joint.prof.%s.png'%(var)
    oimg.save(outPath)
    print outPath
            
 
    
