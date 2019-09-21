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
    figDir      = '/home.rainbow/utsumi/public_html/tempfig/validprof'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/validprof'

else:
    print 'check myhost'
    sys.exit()
#*******************************
lseason=['JJA']

iy =30 # top
ey = -30
ix =10
ex = -10
#*** 2-panels ********
#lmetrix=['rmse','cc']
#lrettype = ['epc','gprof']
#for season in lseason:
#    for metrix in lmetrix:
#        ddat = {}
#        for i,rettype in enumerate(lrettype):
#            figPath= figDir + '/mmap.%s.prof.%s.%s.png'%(metrix,rettype,season)
#            iimg = Image.open(figPath)
#            a2array = asarray(iimg)[iy:ey, ix:ex]
#            ddat[i] = a2array
#    
#        a2oarray = concatenate([ddat[0],ddat[1]],axis=1)
#        oimg    = Image.fromarray(a2oarray)
#        outPath = figDir + '/joint.map.%s.%s.png'%(metrix,season)
#        oimg.save(outPath)
#        print outPath
    

#*** 3-panels ********
lmetrix=['peakh-asl','peakh-agl','peakw','totwat']
lrettype = ['dpr','epc','gprof']
for season in lseason:
    for metrix in lmetrix:
        ddat = {}
        for i,rettype in enumerate(lrettype):
            figPath= figDir + '/mmap.%s.prof.%s.%s.png'%(metrix,rettype,season)
            iimg = Image.open(figPath)
            a2array = asarray(iimg)[iy:ey, ix:ex]
            ddat[i] = a2array
    
        a2oarray = concatenate([ddat[0],ddat[1],ddat[2]],axis=0)
        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.map.%s.%s.png'%(metrix,season)
        oimg.save(outPath)
        print outPath
    
 
    
