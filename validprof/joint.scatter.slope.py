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
lrettype = ['epc','gprof']
lstype = ['sea','veg']
#lstype = ['sea','veg','coast','snow']
#lptype= ['stra','conv']
lptype = ['all']
lph    = ['H','L']
#lph    = ['A','H','L']
#lph   = ['A']
#lptype= ['conv']

iy =50 # top
ey = -40
ix =2
ex = -15

i = -1
ddat = {}
for rettype in lrettype:
    for ph in lph:
        for stype in lstype:
            for ptype in lptype:
                i=i+1
                stamp = '%s.s-%s.p-%s.ph-%s'%(rettype,stype,ptype,ph)

                figPath = figDir + '/scatter.slope.%s.png'%(stamp)
                iimg = Image.open(figPath)
                a2array = asarray(iimg)[iy:ey, ix:ex]
        
                ddat[i] = a2array

                print stamp
                print ddat[i].shape

ddat[-1] = a2array *0 + 255


if i==7:
    a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=0)
    a2line1 = concatenate([ddat[4],ddat[5],ddat[6],ddat[7]],axis=0)
    a2oarray = np.concatenate([a2line0, a2line1], axis=1)

elif i==15:
    a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=1)
    a2line1 = concatenate([ddat[4],ddat[5],ddat[6],ddat[7]],axis=1)
    a2line2 = concatenate([ddat[8],ddat[9],ddat[10],ddat[11]],axis=1)
    a2line3 = concatenate([ddat[12],ddat[13],ddat[14],ddat[15]],axis=1)
    a2oarray = np.concatenate([a2line0, a2line1, a2line2, a2line3], axis=0)



else:
    print 'check i',i
    sys.exit()

oimg    = Image.fromarray(a2oarray)
outPath = figDir + '/joint.scatter.slope.png'
oimg.save(outPath)
print outPath
        
 
    
