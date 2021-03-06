#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

import numpy as np
import pylab as pl
from sklearn.model_selection import train_test_split
import matplotlib.gridspec as gridspec
from glob import glob
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, initializers
from tensorflow.keras.models import model_from_json

import numpy.ma as ma
import sys,os
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
    tankbaseDir = '/tank'
    figDir = '/home/utsumi/temp/stop'

elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    figDir = '/home/utsumi/temp/stop'

else:
    stopbaseDir= '/mnt/j/PMM/stop'
    figDir = '/mnt/c/ubuntu/fig'

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,29)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

lch = [['tir',1],['tir',2],['tir',3],['tir',4],['tir',5],['tir',6],['tir',7],['tir',8],['tir',9],['tir',10],['sir',1],['sir',2]]

nygeo,nxgeo = 7, 7
cygeo,cxgeo = int(nygeo/2), int(nxgeo/2)
miss  = -9999.
argvs = sys.argv
if len(argvs)>1:
    print argvs
    lsurf = argvs[1:]

lDTimeSkip = util.ret_lDTime(datetime(2017,9,25),datetime(2017,9,29),timedelta(days=1))

#lunits = [96,96,64,64]
#lunits = [96,96,64]
lunits = [64,64,64]
#lunits = [64,64]
#expr = 'pnt'
expr = 'cnv'
act  = '04'
nloop = 100
#clipnorm = 1.0e-3
#clipnorm = 1.0e-2
savemodel = 1
restmodel = 1
onlypred  = 1
#***********************************************************
# Functions
#***********************************************************
def append_history(history, histPath):
    hist = history.history
    hist_new = {}
    if os.path.exists(histPath):
        with open(histPath, 'rb') as f:
            hist_old = pickle.load(f)

        keys = history.history
        for key in keys:
            hist_new[key] = hist_old[key] + hist[key]
            
    else:
        hist_new = hist

 
    with open(histPath, 'wb') as f:
        pickle.dump(hist_new, f) 



def split2batchs(a,bsize):
    n=len(a)/bsize + 1
    return [a[i*bsize:(i+1)*bsize] for i in range(n)] 

def build_model_Conv2D():
    nblocks = 3
    for i in range(nblocks):
        model = tf.keras.Sequential()
        model.add( layers.Conv2D(12, kernel_size=(3, 3), strides=(1, 1),
                 padding='same',
                 ))
    
        model.add( layers.Activation('relu'))
        model.add( layers.Conv2D(12, kernel_size=(3, 3), strides=(1, 1),
                 padding='same',
                 ))
        model.add( layers.BatchNormalization())
        model.add( layers.Activation('relu'))
        model.add( layers.MaxPooling2D(pool_size=(2, 2), strides=(2, 2), padding='same'))


    model.add( layers.Flatten())
    model.add( layers.Dense(64) )
    model.add( layers.BatchNormalization())
    model.add( layers.Activation('relu'))

    model.add(layers.Dense(1, activation='sigmoid'))
    #optimizer = tf.keras.optimizers.Adam(clipnorm=clipnorm)
    optimizer = tf.keras.optimizers.Adam()
    #optimizer = tf.keras.optimizers.RMSprop(0.001)

    loss = 'binary_crossentropy'
    model.compile(optimizer=optimizer
                  ,loss = loss
                  ,metrics=['accuracy']
                 )

    return model 



def build_model_2d(lunits):
    nlayers = len(lunits)

    model = tf.keras.Sequential()
    for i in range(nlayers):
        mid_units = lunits[i]

        model.add( layers.Dense(mid_units) )
        model.add( layers.BatchNormalization())
        model.add( layers.Activation('relu'))


    model.add(layers.Dense(1, activation='sigmoid'))
    optimizer = tf.keras.optimizers.Adam(clipnorm=clipnorm)
    #optimizer = tf.keras.optimizers.RMSprop(0.001)

    loss = 'binary_crossentropy'
    model.compile(optimizer=optimizer
                  ,loss = loss
                  ,metrics=['accuracy']
                 )

    return model 


    
def rmse(x,y):
    x = x.flatten()
    y = y.flatten()
    return np.sqrt((((x-y))**2).mean())
def Rmse(x,y):
    Min,Max=stopmin,stopmax
    return np.sqrt( ( ( ((Max-Min)*x+Min).flatten()-((Max-Min)*y+Min).flatten() )**2 ).mean() )
def cc(x,y):
    return np.corrcoef( x.flatten(), y.flatten() )[0,1]
def sort(x):
    return np.sort(x.flatten())
def unit(x):
    return (x-np.min(x,0))/(np.max(x,0)-np.min(x,0))


def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)

#def unit(x):
#    return ( x - np.min(x,0) )/( np.max(x,0) - np.min(x,0) )


def split_data(X,Y, trainfrac=0.9, seed=None):

    a1idx = np.arange(X.shape[0]).astype('int32')
    if seed is not None:
        np.random.seed(seed)

    np.random.shuffle(a1idx)

    nrec = X.shape[0]
    ntrain = int(nrec*trainfrac)
    a1idxt = a1idx[:ntrain]
    a1idxv = a1idx[ntrain:]

    validX = X[a1idxv]
    validY = Y[a1idxv]
    trainX = X[a1idxt]
    trainY = Y[a1idxt]

    return trainX, trainY, validX, validY

def shuffle_data(X,Y, seed=None):

    a1idx = np.arange(X.shape[0]).astype('int32')
    if seed is not None:
        np.random.seed(seed)

    np.random.shuffle(a1idx)

    X = X[a1idx]
    Y = Y[a1idx]

    return X,Y




print 'Define functions'

def read_x_image(lDTime, lch):
    a4out = None
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/%04d/%02d/%02d/*'%(Year,Mon,Day)
        lsrcDir = glob(ssearch)
        
        for srcDir in lsrcDir:
            a4tmp = None
            for [ch,chnum] in lch:
                srcPath = srcDir + '/%s.%02d.npy'%(ch, chnum)
                ain = np.load(srcPath)
                nl,ny,nx = ain.shape
                ain = ain.reshape(nl,ny,nx,1)
                if a4tmp is None:
                    a4tmp = ain
                else:
                    a4tmp = np.concatenate([a4tmp, ain],axis=3)

            #-- Extend ----
            if a4out is None:
                a4out = a4tmp
            else:
                a4out = np.concatenate([a4out, a4tmp], axis=0) 
    return a4out




def read_x_point(lDTime, lch):
    a2out = None
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/%04d/%02d/%02d/*'%(Year,Mon,Day)
        lsrcDir = glob(ssearch)
        
        for srcDir in lsrcDir:
            a2tmp = None
            for [ch,chnum] in lch:
                srcPath = srcDir + '/%s.%02d.npy'%(ch, chnum)
                a1in = np.load(srcPath)[:,cygeo,cxgeo]
                a2in = a1in.reshape(-1,1)
                if a2tmp is None:
                    a2tmp = a2in
                else:
                    a2tmp = np.concatenate([a2tmp, a2in],axis=1)

            #-- Extend ----
            if a2out is None:
                a2out = a2tmp
            else:
                a2out = np.concatenate([a2out, a2tmp], axis=0) 
    return a2out


def read_y(lDTime):
    a1out = None
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/himawari/obt.ptype/%04d/%02d/%02d/*'%(Year,Mon,Day)
        lsrcDir = glob(ssearch)
        
        for srcDir in lsrcDir:
            srcPath = srcDir + '/ptype.npy'
            a1in = np.load(srcPath)
            #-- Extend ----
            if a1out is None:
                a1out = a1in
            else:
                a1out = np.concatenate([a1out, a1in], axis=0) 
    #-- convert ---
    ''' 1: Strat. 2: Conv. 3: Others '''
    a1out = (a1out/10000000).astype('int16')
    return a1out

#***********************************************************
# Read parameter
#***********************************************************
paraDir = tankbaseDir + '/utsumi/PMM/himawari/para'
paraPath= paraDir + '/norm-para.csv'
f=open(paraPath,'r'); lines=f.readlines(); f.close()
dave = {}
dstd = {}
for line in lines[1:]:
    line = line.strip().split(',')
    ch, chnum, ave, std, vmin, vmax = line 
    chnum = int(chnum)
    dave[ch,chnum] = float(ave)
    dstd[ch,chnum] = float(std)

a1ave = np.array([dave[ch,chnum] for (ch,chnum) in lch])
a1std = np.array([dstd[ch,chnum] for (ch,chnum) in lch])
#*********************************************
# Read data
#-----------------
if expr=='pnt':
    trainx = read_x_point(lDTime, lch)
    trainy = read_y(lDTime)
    
    #-- Data screening --
    a1flagx = ma.masked_not_equal(trainx, miss).mask.all(axis=1)
    a1flagy = ma.masked_not_equal(trainy, miss).mask
    
    a1flag = a1flagx * a1flagy
    
    trainx = trainx[a1flag]
    trainy = trainy[a1flag]
    
    #-- Normalize --
    trainx = (trainx -a1ave.reshape(1,-1))/a1std.reshape(1,-1)


elif expr=='cnv':
    trainx = read_x_image(lDTime, lch)
    trainy = read_y(lDTime)

    #-- Data screening --
    a1flagx = ma.masked_not_equal(trainx, miss).mask.all(axis=(1,2,3))
    a1flagy = ma.masked_not_equal(trainy, miss).mask
    
    a1flag = a1flagx * a1flagy
    
    trainx = trainx[a1flag]
    trainy = trainy[a1flag]
    
    #-- Normalize --
    trainx = (trainx -a1ave.reshape(1,1,1,-1))/a1std.reshape(1,1,1,-1)


else:
    print 'check expr=',expr
    sys.exit()


#-- set convective to 1 and strat/others to 0 ---
trainy = ma.masked_equal(trainy,2).mask.astype('int32')

#-- Split test data ------
trainx, trainy, testx, testy = split_data(trainx,trainy,trainfrac=0.9, seed=0)

#-- Correct imbalance ----
nrec = trainx.shape[0]
a1keepflag = np.array([True]*nrec)
a1idxstrat = ma.masked_where(trainy !=0, np.arange(nrec).astype('int32')).compressed()

a1idxremove= np.random.choice(a1idxstrat, int(len(a1idxstrat)*0.6), replace=False)

a1keepflag[a1idxremove] = False

trainx = trainx[a1keepflag]
trainy = trainy[a1keepflag]

#***************************************
# Builad model
#***************************************
cpDir = tankbaseDir + '/utsumi/PMM/himawari/obt.ptype/cp/%s-%s'%(expr,act)
util.mk_dir(cpDir)
#cpPath= cpDir + '/cp.ckpt'
#restcpPath = cpDir + '/cp-093.ckpt'
restcpPath = tf.train.latest_checkpoint(cpDir)

if expr=='pnt':
    model = build_model_2d(lunits)
elif expr=='cnv':
    model = build_model_Conv2D()

if savemodel==1:
    modelPath = cpDir + '/model.json'
    json_string = model.to_json()
    with open(modelPath, 'w') as f:
        f.write(json_string)

if (restmodel ==1)or(onlypred==1): 
    print '**************************'
    print 'Prediction only'
    print '**************************'
    #latest = tf.train.latest_checkpoint(cpDir)
    print ''
    print 'Load weights'
    print restcpPath
    print ''
    model.load_weights(restcpPath)


#***********************************************************
# Only prediction
#***********************************************************
if onlypred ==1:
    pred = model.predict(testx)  # 0-1
    predy= ma.masked_less(pred,0.5).filled(0)
    predy= ma.masked_greater_equal(predy,0.5).filled(1)
    predy=predy.astype('int32')

    Cc  = ((testy==1)&(predy==1)).sum()
    Cs  = ((testy==1)&(predy==0)).sum()
    Sc  = ((testy==0)&(predy==1)).sum()
    Ss  = ((testy==0)&(predy==0)).sum()

    print 'CC\tSc'
    print 'Cs\tSs'
    print Cc,Sc
    print Cs,Ss
    print pred
    print predy
    print ''
    print testy
    sys.exit()
#***********************************************************
# Training Main loop start
#***********************************************************
for iloop in range(nloop):
    #-- Shuffle -----
    trainx, trainy = shuffle_data(trainx, trainy)
    #----------------
    cpPath= cpDir + '/cp-%03d.ckpt'%(iloop)

    cp_callback = tf.keras.callbacks.ModelCheckpoint(cpPath, save_weights_only=True, verbose=1, monitor='accuracy', save_best_only=True, mode='min')
    
    early_stop = tf.keras.callbacks.EarlyStopping(monitor='accuracy', patience=3)

    EPOCHS = 10
    
    history = model.fit(
        x=trainx, y=trainy,
        batch_size=128,
        #epochs=EPOCHS, validation_data=(validX, validY),
        epochs=EPOCHS, validation_split=0.2,
        shuffle=True,
        callbacks = [cp_callback, early_stop],
        verbose=1)

    histPath = cpDir + '/hist.pickle'
    if restmodel ==0:
        if (iloop==0)and(os.path.exists(histPath)):
            os.remove(histPath)
    append_history(history, histPath)
 

