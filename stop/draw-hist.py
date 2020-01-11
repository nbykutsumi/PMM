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
lsurf = ['vege']
#lsurf = ['ocean']
saveprep = 1
savemodel= 1
#restmodel = 1
restmodel = 0
expr  = 'mse04-01'
#lact = ['H','L','LT']
#lact = ['L','LT']
#lact = ['LTQZ']
lact = ['HTQZ']
#lact = ['L']
#***********************************************************
# Functions
#***********************************************************
#***********************************************************
# Main loop start
#***********************************************************

for act in lact:
    cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)
    
    for surf in lsurf:
        histPath = cpDir + '/hist-%s-%s-s%s.pickle'%(expr,act, surf)


        with open(histPath,'rb') as f:
            d = pickle.load(f)

        print d
        epochs = 5
        print d.keys()
        #adif_t = d['mean_absolute_error']
        #adif_v = d['val_mean_absolute_error']
        adif_t = d['mae']
        adif_v = d['val_mae']

        x      = (np.arange(len(adif_t))+1) * epochs

        fig = plt.figure(figsize=(5,3))
        ax  = fig.add_axes([0.15,0.15,0.7,0.7])
        ax.plot(x, adif_t, '-',color='k',label = 'train')
        ax.plot(x, adif_v, '--',color='k',label = 'validation')

        ax.set_ylim([0,4000])
        ax.set_ylabel('mean abs. error')
        ax.set_xlabel('epoch')
        ax.legend()

        stitle = '%s-%s surf=%s'%(expr, act , surf)
        plt.title(stitle)

        figPath = figDir + '/hist.%s-%s-surf-%s.png'%(expr, act,surf)

        plt.savefig(figPath)
        print figPath

