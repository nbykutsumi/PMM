from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket

lobs = ['mrms','cmb']
lsurftype= ['ocean','vegetation','coast','snow']
lrettype = ['glb.v03.minrec1000.maxrec10000','gprof']
lseason = ['ALL']
dregion = {'mrms':'US', 'cmb':'GLB'}

myhost = socket.gethostname()
if myhost == 'shui':
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()

iy =0   # top
ey =-1
ix =40
ex =-10

for season in lseason:
    for obs in lobs:
        region = dregion[obs]
        for rettype in lrettype:
            i = -1
            ddat = {}
            for surftype in lsurftype:
                i = i+1
                figPath = figDir + '/scatter.%s.%s.%s.%s.%s.png'%(region,obs,rettype,surftype,season)
                iimg    = Image.open(figPath)
                if i==0:
                    a2array = asarray(iimg)[iy:ey, 0:-1]
                else:
                    a2array = asarray(iimg)[iy:ey, ix:ex]

                ddat[i] = a2array


            a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=1)
            a2oarray= a2line0 
            oimg    = Image.fromarray(a2oarray)
            outPath = figDir + '/joint.scatter.%s.%s.%s.%s.png'%(region,obs,rettype,season)
            oimg.save(outPath)
            print outPath
