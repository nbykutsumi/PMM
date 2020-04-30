from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket

#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir      = '/home/temp/ret'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir      = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************

iy =30 # top
ey = -50
ix =1
ex = -1


llkey = [
         ['epc.stoprad', 'epc.dstop' ,'gprof-shift.dstop'],
]

lthwat = [0.05, 0.2]

#*** 3-panels ********
for thwat in lthwat:
    for lkey in llkey:
        for var2 in ['ave','std']:
            ddat = {}
            for i,key in enumerate(lkey):
                rettype, var = key.split('.')
                figDir = '/home/utsumi/temp/ret'
                figPath= figDir + '/map.%s.%s.th-%.3f.%s.png'%(rettype,var,thwat,var2)
                iimg = Image.open(figPath)
                if i==len(lkey)-1:
                    print 'last',key,var2
                    a2array = asarray(iimg)[iy:-1, ix:ex]
                else:
                    a2array = asarray(iimg)[iy:ey, ix:ex]
                ddat[i] = a2array
            
            a2oarray = concatenate([ddat[0],ddat[1],ddat[2]],axis=0)
            oimg    = Image.fromarray(a2oarray)
            outPath = figDir + '/joint.map.orbit.%s.th-%.3f.%s.png'%(var,thwat,var2)
            oimg.save(outPath)
            print outPath
        
 
    
