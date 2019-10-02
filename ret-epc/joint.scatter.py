from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket

lseason = ['JJA']
lsurftype= ['ocean','vegetation','coast']
lrettype = ['NS','MS','NScmb','MScmb','GPROF']
expr = 'glb.minrec1000.maxrec10000'
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
ix =0
ex =-1

for season in lseason:
    i = -1
    ddat = {}
    for rettype in lrettype:
        for surftype in lsurftype:
            i = i+1
            figPath = figDir + '/scatter.%s.%s.%s.%s.png'%(expr,rettype,surftype,season)
            iimg    = Image.open(figPath)
            a2array = asarray(iimg)[iy:ey, ix:ex]
            ddat[i] = a2array    


    a2line0 = concatenate([ddat[0],ddat[1],ddat[2]],axis=0)
    a2line1 = concatenate([ddat[3],ddat[4],ddat[5]],axis=0)
    a2line2 = concatenate([ddat[6],ddat[7],ddat[8]],axis=0)
    a2line3 = concatenate([ddat[9],ddat[10],ddat[11]],axis=0)
    a2line4 = concatenate([ddat[12],ddat[13],ddat[14]],axis=0)

    a2oarray= concatenate([a2line0,a2line1,a2line2,a2line3,a2line4],axis=1)

    oimg    = Image.fromarray(a2oarray)
    outPath = figDir + '/joint.scatter.log.%s.png'%(season)
    oimg.save(outPath)
    print outPath
