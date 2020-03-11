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
ey = -60
ix =10
ex = -10

llkey = [['epc.peakhrad','epc.dpeakh','epc.dpeakh'],
         ['epc.peakvrad','epc.dpeakv','epc.dpeakv'],
         ['epc.stoprad', 'epc.dstop' ,'epc.dstop'],
         ['epc.condrad', 'epc.dcond' ,'epc.dcond'],
]


#*** 3-panels ********
for lkey in llkey:
    for var2 in ['ave','std']:
        ddat = {}
        for i,key in enumerate(lkey):
            rettype, var = key.split('.')
            figDir = '/home/utsumi/temp/ret'
            figPath= figDir + '/map.%s.%s.%s.png'%(rettype,var,var2)
            iimg = Image.open(figPath)
            if i==len(lkey)-1:
                print 'last',key,var2
                a2array = asarray(iimg)[iy:-1, ix:ex]
            else:
                a2array = asarray(iimg)[iy:ey, ix:ex]
            ddat[i] = a2array
        
        a2oarray = concatenate([ddat[0],ddat[1],ddat[2]],axis=0)
        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.map.orbit.%s.%s.png'%(var,var2)
        oimg.save(outPath)
        print outPath
    
 
    
