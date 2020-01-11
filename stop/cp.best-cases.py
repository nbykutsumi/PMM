import numpy as np
import pylab as pl
import glob

import numpy.ma as ma
import sys,os,shutil
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
#from sklearn.decomposition import PCA
import calendar
import socket
import pickle
import random
#get_ipython().magic(u'matplotlib inline')


hostname = socket.gethostname()
if hostname == 'shui':
    #stopbaseDir= '/work/hk01/utsumi/PMM/stop'
    stopbaseDir= '/tank/utsumi/PMM/stop'
    #figDir = '/home.rainbow/utsumi/public_html/tempfig/stop'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    stopbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'
    #figDir = '/home/utsumi/mnt/lab_home_rainbow/utsumi/public_html/tempfig/stop'
    figDir = '/home/utsumi/temp/stop'

else:
    stopbaseDir= '/mnt/j/PMM/stop'
    figDir = '/mnt/c/ubuntu/fig'


Year  = 2017
season = 0

argvs = sys.argv
if len(argvs)>1:
    print argvs
    lsurf = argvs[1:]
else:
    #lsurf = ['ocean','seaice','vege','snow','swater','coast','siedge']
    #lsurf = ['vege']
    lsurf = ['siedge']
    #lsurf = ['seaice','vege']
    #lsurf = ['snow','swater','coast','siedge']
    #lsurf = ['vege']    

expr = 'best01'
act  = 'HTQZ'

lsurfTmp = []
for x in lsurf:
    if x.isdigit(): x = int(x)
    lsurfTmp.append(x)
lsurf = lsurfTmp


dcase={
    'ocean':['mse09-01-10'],
    'seaice':['mse10-01-10'],
    'vege':['mse09-01-08'],
    'snow':['mse09-01-11'],
    'swater':['mse10-01-04'],
    'coast':['mse09-01-07'],
    'siedge':['mse10-01-14']
    }

season = 0

for surf in lsurf:
    [exprTmp]  = dcase[surf]

    icpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(exprTmp,act, season) 
    ocpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season) 
    util.mk_dir(ocpDir)
    #-- Check point files --
    ssearch = icpDir + '/*%s*'%(surf)
    lsrcPath = glob.glob(ssearch)
    for srcPath in lsrcPath:
        print srcPath
        shutil.copy(srcPath, ocpDir)    


    print '' 
    print surf 
    print icpDir
    print ocpDir



