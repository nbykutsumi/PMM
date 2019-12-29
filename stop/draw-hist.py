import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pylab as pl
from glob import glob
import numpy.ma as ma
import sys,os
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
#from sklearn.decomposition import PCA
import calendar
import socket
import pickle
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
#season = 1
#lisurf = range(1,14+1)
#lisurf = range(2,14+1)
#lisurf = range(4,14+1)
lisurf = [7]
saveprep = 1
savemodel= 1
#restmodel = 1
restmodel = 0
expr  = 'a01'
#lact = ['H','L','LT']
#lact = ['L','LT']
lact = ['LTQZ']
#lact = ['L']
#***********************************************************
# Functions
#***********************************************************
#***********************************************************
# Main loop start
#***********************************************************

for act in lact:
    cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)
    
    for isurf in lisurf:
        histPath = cpDir + '/hist-%s-%s-s%02d.pickle'%(expr,act, isurf)


        with open(histPath,'rb') as f:
            d = pickle.load(f)

        print d
        epochs = 5
        print d.keys()
        adif_t = d['mean_absolute_error']
        adif_v = d['val_mean_absolute_error']
        x      = (np.arange(len(adif_t))+1) * epochs

        fig = plt.figure(figsize=(5,3))
        ax  = fig.add_axes([0.15,0.15,0.7,0.7])
        ax.plot(x, adif_t, '-',color='k',label = 'train')
        ax.plot(x, adif_v, '--',color='k',label = 'validation')

        ax.set_ylim([0,2000])
        ax.set_ylabel('mean abs. error')
        ax.set_xlabel('epoch')
        ax.legend()

        stitle = '%s-%s surf=%d'%(expr, act , isurf)
        plt.title(stitle)

        figPath = figDir + '/hist.%s-%s-surf%d.png'%(expr, act,isurf)

        plt.savefig(figPath)
        print figPath

