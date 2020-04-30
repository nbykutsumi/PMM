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
lsensor = ['GMI','AMSR2','SSMIS','ATMS','MHS']
#lstype = ['ocean','vegetation','coast','snow']
lstype = ['ocean','vegetation','coast']
lrettype= ['EPC','GPROF']
#lrettype= ['GPROF']
lthwat   = [0.033, 0.05, 0.2]
coefflag = 'nocoef'  # 'nocoef', 'wcoef'
iy =1 # top
ey = -1
ix =0
ex = -5
#**********
dsatname = {'GMI':'GPM', 'AMSR2':'GCOMW','SSMIS':'ALLSATE','ATMS':"ALLSATE",'MHS':'ALLSATE'}

for thwat in lthwat:
    for stype in lstype:
        ddat = {}
        i = -1
        for rettype in lrettype:
            for sensor in lsensor:
                i=i+1
                satname= dsatname[sensor]
                figPath= figDir + '/scatter.height.prf.%s.smp1000.%s.%s.%s.th-%.3f.png'%(sensor,coefflag,rettype,stype, thwat)
                iimg = Image.open(figPath)
                a2array = asarray(iimg)[iy:ey, ix:ex]

                ddat[i] = a2array

        ddat[-1] = a2array *0 + 255

        a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3],ddat[4]],axis=1)
        a2line1 = concatenate([ddat[5],ddat[6],ddat[7],ddat[8],ddat[9]],axis=1)
        a2oarray = np.concatenate([a2line0, a2line1], axis=0)

        oimg    = Image.fromarray(a2oarray)
        outPath = figDir + '/joint.height.%s.th-%.3f.png'%(stype, thwat)
        oimg.save(outPath)
        print outPath
            
 
    
