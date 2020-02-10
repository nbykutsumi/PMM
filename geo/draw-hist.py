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
    tankbaseDir = '/tank'
    figDir = '/home/utsumi/temp/geo'
elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    figDir = '/home/utsumi/temp/geo'

else:
    stopbaseDir= '/mnt/j/PMM/stop'
    figDir = '/mnt/c/ubuntu/fig'


Year  = 2017
#expr  = 'pnt'
expr  = 'cnv'
act   = '04'
#***********************************************************
# Functions
#***********************************************************
#***********************************************************
# Main loop start
#***********************************************************
cpDir = tankbaseDir + '/utsumi/PMM/himawari/obt.ptype/cp/%s-%s'%(expr,act)
histPath = cpDir + '/hist.pickle'

with open(histPath,'rb') as f:
    d = pickle.load(f)

print d
epochs = 10
print d.keys()
a_t = d['accuracy']
a_v = d['val_accuracy']

x      = (np.arange(len(a_t))+1) * epochs

fig = plt.figure(figsize=(5,3))
ax  = fig.add_axes([0.15,0.15,0.7,0.7])
ax.plot(x, a_t, '-',color='k',label = 'train')
ax.plot(x, a_v, '--',color='k',label = 'validation')

ax.set_ylim([0.75,0.8])
ax.set_ylabel('accuracy')
ax.set_xlabel('epoch')
ax.legend()

stitle = '%s-%s'%(expr, act)
plt.title(stitle)

figPath = figDir + '/hist.%s-%s.png'%(expr, act)

plt.savefig(figPath)
print figPath

