from numpy import *
import os, sys
from PIL import Image
import myfunc.util as util
import socket
from datetime import datetime, timedelta

expr = 'glb.v03.minrec1000.maxrec10000'
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
iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,15)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    figPath = {}
    figPath[0]= figDir + '/mmap.uncond.prec.dpr.%04d%02d%02d.png'%(Year,Mon,Day)
    figPath[1]= figDir + '/mmap.uncond.prec.%s.epcNScmb.%04d%02d%02d.png'%(expr,Year,Mon,Day)
    figPath[2]= figDir + '/mmap.uncond.prec.gprof.%04d%02d%02d.png'%(Year,Mon,Day)
    figPath[3]= figDir + '/mmap.uncond.difprec.epcNScmb.%s.%04d%02d%02d.png'%(expr,Year,Mon,Day)
    figPath[4]= figDir + '/mmap.uncond.difprec.gprof.%04d%02d%02d.png'%(Year,Mon,Day)

    ddat = {}
    for i in range(5):
        iimg    = Image.open(figPath[i])
        a2array = asarray(iimg)[iy:ey, ix:ex]
        ddat[i] = a2array    

    ddat[-1]= a2array*0 + 255
    a2line0 = concatenate([ddat[0],ddat[1],ddat[2]],axis=1)
    a2line1 = concatenate([ddat[-1],ddat[3],ddat[4]],axis=1)

    a2oarray= concatenate([a2line0,a2line1],axis=0)

    oimg    = Image.fromarray(a2oarray)
    outPath = figDir + '/joint.mmap.uncond.daily.%04d%02d%02d.png'%(Year,Mon,Day)
    oimg.save(outPath)
    print outPath
