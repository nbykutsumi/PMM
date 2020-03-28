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
lvar = ['stop','cond','dstop','dcond','cc']
lvar_byself = ['stop','cond']
lvar_byrad  = ['dstop','dcond','cc']
var2 = 'ave'
for var in lvar:
    #**** 3x2 ******
    if var in lvar_byself:
        ddat = {}
        ldattype = ['cmb','epc','gprof-shift']
        i= -1
        for idattype,dattype in enumerate(ldattype):
            for ptype in ['stra','conv']:
                i = i+1
                if dattype=='cmb':
                    rettype='epc'
                    var1   = var+'rad'
                    bywhat ='rad'
                else:
                    rettype=dattype
                    var1   = var+'pmw'
                    bywhat ='pmw'
    
                figDir = '/home/utsumi/temp/ret'
                figPath= figDir + '/map.%s.%s.by-%s-%s.%s.png'%(rettype,var1,bywhat,ptype,var2)
                iimg = Image.open(figPath)
                if idattype==len(ldattype)-1:
                    a2array = asarray(iimg)[iy:-1, ix:ex]
                else:
                    a2array = asarray(iimg)[iy:ey, ix:ex]
                ddat[i] = a2array
            
        line0 = np.concatenate([ddat[0],ddat[1]],axis=1)
        line1 = np.concatenate([ddat[2],ddat[3]],axis=1)
        line2 = np.concatenate([ddat[4],ddat[5]],axis=1)
        a2oarray  = np.concatenate([line0,line1,line2],axis=0)
        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.map.orbit.%s.%s.by-ptype.png'%(var1,var2)
        oimg.save(outPath)
        print outPath

    #**** 2x2 ******
    if var in lvar_byrad:
        ddat = {}
        i= -1
        lrettype = ['epc','gprof-shift']
        for irettype,rettype in enumerate(lrettype):
            for ptype in ['stra','conv']:
                i = i+1
                figDir = '/home/utsumi/temp/ret'
                figPath= figDir + '/map.%s.%s.by-rad-%s.%s.png'%(rettype,var,ptype,var2)
                iimg = Image.open(figPath)
                if idattype==len(ldattype)-1:
                    a2array = asarray(iimg)[iy:-1, ix:ex]
                else:
                    a2array = asarray(iimg)[iy:ey, ix:ex]
                ddat[i] = a2array
            
        line0 = np.concatenate([ddat[0],ddat[1]],axis=1)
        line1 = np.concatenate([ddat[2],ddat[3]],axis=1)
        a2oarray  = np.concatenate([line0,line1],axis=0)
        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.map.orbit.%s.%s.by-ptype.png'%(var,var2)
        oimg.save(outPath)
        print outPath
 
    
