import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
import epcfunc
import math, time
#get_ipython().magic(u'matplotlib inline')


hostname = socket.gethostname()
if hostname == 'shui':
    tankbaseDir = '/tank'
    stopbaseDir= '/tank/utsumi/PMM/stop'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    stopbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'
    figDir = '/home/utsumi/temp/stop'

else:
    print 'check hostname',hostname
    sys.exit()


epcbaseDir = tankbaseDir + '/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117'

Year  = 2017

#argvs = sys.argv
lidx_db  = range(29*29*29)
onlypred = 0
savemodel= 1
restmodel = 1
#restmodel = 0
#restmodel = 1

nloop = 20

#expr = 'epc-mse-01-01'
#expr = 'epc-mse-01-02' # BS=128
#expr = 'epc-mse-01-02'  # BS=32
#expr = 'epc-mse-01-03'  # BS=32, Dropout
expr = 'epc-mse-02-00' # BS=32, Dropout, [128,96,64,64,64,64]
#lact = ['L','LT']
#lact = ['LTQZ']
#lact = ['HTQZ']
#lact = ['HTQZ']
#lact = ['HTQZ']
lact = ['HTQZ']
#lact = ['H']


daddvar = {'T':['t2m'], 'Q':['tqv'], 'Z':'gtopo'}
dvarminmax = {'t2m':[200,350],'tqv':[0,120],'gtopo':[-100,8800]}
stopmin = 0
stopmax = 32000
#clipnorm = 1.
#clipnorm = 0.01
clipnorm = 1e-3

dncol ={'H'  : 12
       ,'HT' : 12 + 2
       ,'HTQ': 13 + 2
       ,'HTZ': 13 + 2
       ,'HTQZ':13 + 3
       } 

#lunits = [96,96,64,64]
#lunits = [96,96,64]
#lunits = [64,64,64]
#lunits = [64,64,64]
lunits = [128,96,64,64,64,64]
#lunits = [64,64]


#***********************************************************
# Functions
#***********************************************************
def ret_a1idx(nrec, trainfrac=0.8):
    aidx_all = np.arange(nrec).astype('int32')
    np.random.seed(0)
    np.random.shuffle(aidx_all)
    nrec_train = int(math.floor(nrec*trainfrac))
    aidx_train = aidx_all[:nrec_train]
    if nrec_train ==nrec:
        aidx_test = np.array([])
    else:
        aidx_test = aidx_all[nrec_train:]

    return aidx_train, aidx_test

def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return np.array(lout)



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


def build_model_2d(ncol, lunits):

    nlayers = len(lunits)

    model = tf.keras.Sequential()
    #model.add( layers.InputLayer(input_shape=(ncol,)) )
    #model.add( layers.InputLayer(input_shape=(192,)) )
    for i in range(nlayers):
        mid_units = lunits[i]
        print i,'mid_units=',mid_units

        #model.add( layers.Dense(mid_units, activation='relu',kernel_initializer=initializers.he_normal()) )
        #model.add( layers.Dropout(rate=0.3) )

        model.add( layers.Dense(mid_units, kernel_initializer=initializers.he_normal()) )
        model.add( layers.BatchNormalization())
        model.add( layers.Activation('relu'))
        model.add( layers.Dropout(0.4))


    model.add(layers.Dense(1))
    optimizer = tf.keras.optimizers.Adam(clipnorm=clipnorm)
    #optimizer = tf.keras.optimizers.RMSprop(0.001)


    if 'wtpdf' in expr:
        params = dpdfparam[surf]
        maxpdf = dpdfmax[surf]
        loss   = loss_pdf(params, maxpdf, b=40)
    elif 'wtmape' in expr:
        params = dpdfparam[surf]
        maxpdf = dpdfmax[surf]
        loss   = loss_wtmape(params, maxpdf)
    elif 'male' in expr:
        loss   = loss_mean_abs_log_error

    elif 'mape' in expr:
        loss   = 'mape'
    elif 'msle' in expr:
        loss   = 'mean_squared_logarithmic_error'
    elif 'mse' in expr:
        loss   = 'mse'
    else:
        loss   = 'mse'

    model.compile(optimizer=optimizer
                  ,loss = loss
                  ,metrics=['mae']
                 )


    return model 




def Figure(Label, Prediction, bins, figPath):
    Min,Max = stopmin, stopmax
    #recover_testY = (Max-Min)*Label.flatten()      + Min
    #recover_pred  = (Max-Min)*Prediction.flatten() + Min

    recover_testY = Label.flatten()
    recover_pred  = Prediction.flatten()



    pl.figure(figsize=(15,15))
    gs = gridspec.GridSpec(2,2, width_ratios=[1,1], height_ratios=[1,1])
    
    pl.subplot(gs[0,:])
    pl.plot(recover_testY/1000., c='r', label ='Observation')
    pl.plot(recover_pred /1000., c='b', label ='Prediction')
    pl.ylabel('height(km)')
    pl.ylim([0,18])
    pl.legend()
    pl.title('%s act=%s'%(expr, act))
    print('RMSE:'     , np.round(rmse(Label, Prediction) , 4))
    print('real RMSE:', np.round(Rmse(Label, Prediction) , 4))
    print('CC:'       , np.round(  cc(Label, Prediction) , 4))
    
    pl.subplot(gs[2]) # values prediction and testY are between -4 and 4
    aa = recover_pred
    bb = recover_testY
    interval           = np.array([ Min + (Max - Min)/bins*i for i in range(bins+1) ])
    interval1          = np.array([ Min + (Max - Min)/bins*i for i in range(bins+1) ])
    revised_interval   = interval[:-1]  + (Max - Min)/(2*bins)
    revised_interval1  = interval1[:-1] + (Max - Min)/(2*bins)
    cumulative_number  = []
    cumulative_number1 = []
    for i in range(bins):
        cumulative_number.append(  (aa < interval[i+1] ).sum() - (aa < interval[i] ).sum() )
        cumulative_number1.append( (bb < interval1[i+1]).sum() - (bb < interval1[i]).sum() )
    pl.plot(revised_interval/.1000          , cumulative_number   , color='green', alpha=0.5, label='Prediction')    
    pl.fill_between(revised_interval/.1000  , cumulative_number, 0, color='green', alpha=0.5)
    pl.plot(revised_interval1/.1000         , cumulative_number1  , color='red'  , alpha=0.5 ,label='Observation')    
    pl.fill_between(revised_interval1/.1000 ,cumulative_number1, 0, color='red'  , alpha=0.5)
    pl.ylabel('number of samples')
    pl.xlabel('height(km)')
    pl.legend() 
    pl.title('Distribution')
    pl.legend()

    #*** 2D histogram **********
    H,xbnd,ybnd = np.histogram2d(recover_testY/1000, recover_pred/1000, bins=[np.arange(0,20,0.5), np.arange(0,20,0.5)])
    H = ma.masked_equal(H,0)
    X,Y = np.meshgrid(xbnd,ybnd)
    pl.subplot(gs[3])
    pl.pcolormesh(X,Y,H.T, cmap='jet')
    pl.axis([0,18,0,18])
    pl.xticks([0,5,10,15])
    pl.yticks([0,5,10,15])
    pl.plot([0,35],[0,35],'k')
    pl.xlabel('Observation(km)')
    pl.ylabel('Prediction(km)')
    pl.title('CC=%.2f'%(corr))
    pl.grid()
    pl.colorbar()
    

    pl.savefig(figPath)
    pl.clf()
    print figPath


    
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
def load_var(varName, lidat):
    adat = np.concatenate([np.load(epcbaseDir + '/batch.%s/%02d.npy'%(varName, idat)) for idat in lidat], axis=0)
    return adat

def load_data_2d(lidat):
    aepc  = load_var('epc', lidat)
    astop = load_var('Ku_NS_heightStormTop', lidat)
    print aepc.shape, aepc.min(), aepc.max()


    lvar  = []
    dvarx = {}
    if 'T' in act:
        lvar.append('t2m')
        dvarx['t2m'] = load_var('t2m',lidat)
   
    if 'Q' in act:
        lvar.append('tqv')
        dvarx['tqv'] = load_var('tqv',lidat)

    if 'Z' in act:
        lvar.append('gtopo')
        dvarx['gtopo'] = load_var('gtopo',lidat)


    #****************************************************
    # Screen invalid data
    #****************************************************
    a1flag = ma.masked_equal(load_var('flagInvalidTc', lidat), 0).mask  # invalid=1
    a1flag = a1flag * ~ma.masked_invalid(aepc).mask.any(axis=1)
    a1flag = a1flag * ~ma.masked_invalid(astop).mask
    a1flag = a1flag * ma.masked_inside(astop, stopmin, stopmax).mask

    for var in lvar:
        varmin, varmax = dvarminmax[var]
        a1flag = a1flag * ~ma.masked_invalid(dvarx[var]).mask 
        a1flag = a1flag * ma.masked_inside(dvarx[var],varmin, varmax).mask
        
    
    aepc  = aepc [a1flag]
    astop = astop[a1flag]
    for var in lvar:
        dvarx[var] = dvarx[var][a1flag]

    #***********************************************************
    # Normalize
    #***********************************************************
    aepc  = my_unit(aepc, a1pc_min, a1pc_max)
    #astop = my_unit(astop, stopmin, stopmax)
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        dvarx[var] = my_unit(dvarx[var], varmin, varmax)

    print 1,aepc.shape
    #***********************************************************
    # Stack
    #***********************************************************
    for var in lvar:
        aepc= np.hstack([aepc, dvarx[var].reshape(-1,1)])

    print 2,aepc.shape


  
    print 3,aepc.shape
    return aepc, astop


def split_data(X,Y, trainfrac=0.8):
    nrec = X.shape[0]
    ntrain = int(nrec*trainfrac)
    trainX = X[:ntrain]
    trainY = Y[:ntrain]
    validX = X[ntrain:]
    validY = Y[ntrain:]
    return trainX, trainY, validX, validY

def shuffle_data(X,Y):
    aidxTmp = np.arange(X.shape[0])
    np.random.seed(int(time.time()))
    np.random.shuffle(aidxTmp)
    X = X[aidxTmp]
    Y = Y[aidxTmp]
    print aidxTmp
    print aidxTmp.shape
    return X, Y

print 'Define functions'
#***********************************************************
# Read nrec list
#***********************************************************
nrecDir = epcbaseDir + '/list'
nrecPath= nrecDir + '/nrec.csv'

f=open(nrecPath,'r'); lines=f.readlines(); f.close()
a1nrec = []
for line in lines:
    a1nrec.append( int(line.strip().split(',')[1]) )
a1nrec = np.array(a1nrec)

#***********************************************************
# Read EPC parameter
#***********************************************************
#-- Read EPC range files --
coefDir  = tankbaseDir + '/utsumi/PMM/EPCDB/EPC_COEF/GMI'
rangePath = coefDir + '/PC_MIN_MAX_29.txt'
a2pc_edge = read_table(rangePath)

a1pc_min = a2pc_edge[:,0]
a1pc_max = a2pc_edge[:,-1]


##-- Read PC ave and std file --
#pcavePath  = coefDir + '/ave_pc.txt'
##pcavePath  = '/home/utsumi/bin/ENSPR/ave_pc.txt'
#a2pc_avestd= read_table(pcavePath)
#a1pc_ave   = a2pc_avestd[:,1]
#a1pc_std   = a2pc_avestd[:,2]

#***********************************************************
# Main loop start
#***********************************************************

for act in lact:
    #testX, testY = load_data_2d([8,9])

    cpDir = stopbaseDir + '/cp/%s-%s'%(expr,act)
    util.mk_dir(cpDir)
    cpPath= cpDir + '/cp.ckpt'
    
    ncol  = dncol[act]  
    model = build_model_2d(ncol, lunits)
    #model = build_model_cnn([ny,nx,nc])
    if savemodel==1:
        modelPath = cpDir + '/model.json'
        json_string = model.to_json()

        #if os.path.exists(modelPath):
        #    print 'model exists'
        #    print modelPath
        #    sys.exit()
            
 
    with open(modelPath, 'w') as f:
        f.write(json_string)

    if (restmodel ==1)or(onlypred==1): 
        #latest = tf.train.latest_checkpoint(cpDir)
        if os.path.exists(cpPath):
            print ''
            print 'Load weights'
            print cpPath
            print ''
            model.load_weights(cpPath)
        else:
            print ''
            print 'No weight filea'
            print cpPath
            print ''
   
    #************************
    if onlypred==1:
        print ''
        print ''
        print '***********************'
        print 'Only prediction'
        print '***********************'
        print ''
        print ''
        testX, testY = load_data_2d([8,9])

        pred = model.predict(testX)
        corr = np.corrcoef(testY.flatten(), pred.flatten())[0,1]
        print corr
        figPath = figDir + '/train.%s-%s.png'%(expr,act)
        Figure(testY, pred, 50, figPath)
        sys.exit()    
    #************************
    #************************
    for iloop in range(nloop):
        lidat = np.arange(8)
        np.random.seed(int(time.time()))
        np.random.shuffle(lidat)
        for iloop2 in range(8):
            trainX, trainY = load_data_2d([iloop2])
            if iloop2==7:
                figX, figY = load_data_2d([0])
            else:
                figX, figY = load_data_2d([iloop2+1])

            nrecfig = int(figX.shape[0])
            nrecfig = int(nrecfig*0.2)
            figX = figX[:nrecfig]
            figY = figY[:nrecfig]
 
            print 'trainX.shape=',trainX.shape
            trainX, trainY = shuffle_data(trainX, trainY)
            print 'trainX.shape=',trainX.shape
            #****************************************************
            # Training
            #****************************************************
            cp_callback = tf.keras.callbacks.ModelCheckpoint(cpPath, save_weights_only=True, verbose=1, monitor='val_loss', save_best_only=True, mode='min')
        
            early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)
        
            EPOCHS = 1 
 
            history = model.fit(
                x=trainX, y=trainY,
                #batch_size=128,
                batch_size=32,
                #epochs=EPOCHS, validation_data=(validX, validY),
                epochs=EPOCHS, validation_split=0.2,
                shuffle = True,
                callbacks = [cp_callback, early_stop],
                verbose=1)
        
            #history = model.fit(
            #    x=trainX, y=trainY,
            #    #batch_size=128,
            #    batch_size=256,
            #    epochs=EPOCHS, validation_split=0.2,
            #    callbacks = [cp_callback, early_stop],
            #    verbose=1)
                
        
            histPath = cpDir + '/hist-%s-%s.pickle'%(expr,act)
            if restmodel ==0:
                if (iloop==0)and(os.path.exists(histPath)):
                    
                    os.remove(histPath)
        
            append_history(history, histPath)
        
        
            #*******************************************************
            # Figure
            #*******************************************************
            model.load_weights(cpPath)
        
            pred = model.predict(figX)
            #expr= 'a-%s.s-%02d'%(act,surf)
            corr = np.corrcoef(figY.flatten(), pred.flatten())[0,1]
            print corr
            figPath = figDir + '/train.%s-%s.png'%(expr,act)
            Figure(figY, pred, 50, figPath)
        
            #***************************************************

