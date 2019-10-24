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
lseason=['JJA','DJF','JJADJF']
#lseason = [6,7]
#lregion = ['AMZ','CUS','EUS','TIB','NETP','SETP','NTA','STA','WTP','ETI','WMP','WMA','TAF','NEA','SEC','NIN']
#lregion = ['CUS','SEC','AMZ','TAF','WMP','WMA','WTP','NTA']
lregion = ['TRO','MIDN']
#lstype= ['sea','land','veg','snow','coast','all']
lstype = ['sea','veg','snow']
lptype= ['conv','stra']
#lptype= ['conv']
#lprrange=[[0.5,999],[1,3],[8,12]]
lprrange=[[1,3],[8,12]]

iy =30 # top
ey = -30
ix =10
ex = -10
lvar = ['prof','profstd','profnum']
#*** 3-panels ********
for var in lvar:
    for prrange in lprrange:
        for season in lseason:
            thpr0, thpr1 = prrange
            i = -1
            ddat = {}
            for region in lregion:
                for stype in lstype:
                    for ptype in lptype:
                        i=i+1
                        stamp = 's-%s.p-%s.pr-%.1f-%.1f.%s'%(stype,ptype,thpr0,thpr1,season)

                        figPath = figDir + '/%s.%s.%s.png'%(var,stamp,region)
                        iimg = Image.open(figPath)
                        a2array = asarray(iimg)[iy:ey, ix:ex]
                
                        ddat[i] = a2array

            ddat[-1] = a2array *0 + 255
            if i==15:
                a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3],ddat[4],ddat[5],ddat[6],ddat[7]],axis=1)
                a2line1 = concatenate([ddat[8],ddat[9],ddat[10],ddat[11],ddat[12],ddat[13],ddat[14],ddat[15]],axis=1)
                a2oarray = np.concatenate([a2line0, a2line1], axis=0)

            elif i==11:
                a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3],ddat[4],ddat[5]],axis=1)
                a2line1 = concatenate([ddat[6],ddat[7],ddat[8],ddat[9],ddat[10],ddat[11]],axis=1)
                a2oarray = np.concatenate([a2line0, a2line1], axis=0)


            elif i==7:
                a2line0 = concatenate([ddat[0],ddat[1],ddat[2],ddat[3]],axis=1)
                a2line1 = concatenate([ddat[4],ddat[5],ddat[6],ddat[7]],axis=1)
                a2oarray = np.concatenate([a2line0, a2line1], axis=0)
            else:
                print 'check i',i
                sys.exit()

            oimg    = Image.fromarray(a2oarray)
            outPath = figDir + '/joint.%s.%s.p-%.1f-%.1f.png'%(var,season,thpr0,thpr1)
            oimg.save(outPath)
            print outPath
        
 
    
