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
import optuna
#from optuna.integration import KerasPruningCallback
from optuna.integration.tfkeras import TFKerasPruningCallback

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


argvs = sys.argv

Year  = 2017
season = 0
#season = 9
if len(argvs)>1:
    print argvs
    lisurf = map(int,argvs[1:])
else:

    #lisurf = range(1,14+1)
    #lisurf = range(2,5+1)
    #lisurf = range(6,9+1)
    #lisurf = range(10,13+1)
    #lisurf = range(14,14+1)
    lisurf = [1]

savemodel= 0
restmodel = 0
#restmodel = 1
#expr  = 'a01'
expr  = 'opt01'
#lact = ['H','L','LT']
#lact = ['L','LT']
lact = ['LTQZ']
#lact = ['L']
daddvar = {'T':['t2m'], 'Q':['tqv'], 'Z':'gtopo'}
dvarminmax = {'t2m':[200,350],'tqv':[0,120],'gtopo':[-100,8800]}
dncol ={'H': 3*7*13
       ,'L': 3*7*9
       ,'LT':3*7*9 + 1
       ,'LTQ':3*7*9 + 2
       ,'LTZ':3*7*9 + 2
       ,'LTQZ':3*7*9 + 3
       } 

dnc   ={'H': 13
       ,'L': 9
       ,'LT':9 + 1
       ,'LTQ':9 + 2
       ,'LTZ':9 + 2
       ,'LTQZ':9 + 3
       } 


tcmin, tcmax = 50, 350
stopmin = 0
stopmax = 32000

ldy   = [-1,0,1]
ldx   = [-3,-2,-1,0,1,2,3]
ny = len(ldy)
nx = len(ldx)
cy = (ny-1)/2
cx = (nx-1)/2

lDTimeSkip = util.ret_lDTime(datetime(2017,9,25),datetime(2017,9,29),timedelta(days=1))
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




def build_model_2d_opt(trial):

    nlayers = trial.suggest_int('nlayers', 2,5)
    
    model = keras.Sequential()
    for i in range(nlayers):
        mid_units = int(trial.suggest_discrete_uniform('mid_units_%d'%(i), 32,96,32))

        #model.add( layers.Dense(mid_units, kernel_initializer=initializers.he_normal()) )
        #model.add( layers.BatchNormalization())
        #model.add( layers.Activation('relu'))

        model.add( layers.Dense(mid_units, activation='relu',kernel_initializer=initializers.he_normal()) )
        model.add( layers.Dropout(rate=0.3) )

    model.add(layers.Dense(1))
    #optimizer = tf.keras.optimizers.RMSprop(0.001)
    optimizer = tf.keras.optimizers.Adam()

    model.compile(loss='mse'
                  ,optimizer=optimizer
                  ,metrics=['mae', 'mse']
                 )
    return model 



def build_model_2d():
    model = keras.Sequential([
         layers.Dense(64, activation='relu', input_shape=[ncol])
        ,layers.Dense(64, activation='relu')
        ,layers.Dense(1)
    ])
    optimizer = tf.keras.optimizers.RMSprop(0.001)

    model.compile(loss='mse'
                  ,optimizer=optimizer
                  ,metrics=['mae', 'mse']
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
    pl.title('isurf=%d act=%s box=%dx%d'%(isurf, act,  len(ldy),len(ldx)))
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
    
    #pl.scatter(recover_testY/1000, recover_pred/1000,s=3)
    #pl.plot(np.arange(18000)/1000,np.arange(18000)/1000,c='black',linestyle = ':')
    #pl.axis([0,18,0,18])
    #pl.xticks([0,5,10,15])
    #pl.yticks([0,5,10,15])
    #pl.xlabel('Observation(km)')
    #pl.ylabel('Prediction(km)')
    #pl.title('Correlation')
    #pl.grid()

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



def read_Tc_4d(lDTime=None, ldydx=None, isurf=None, samplerate=None, ch='LH'):
    #a2tc = deque([])
    a4tc = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        #a2tcTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            #srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcDir = stopbaseDir + '/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcPath1=srcDir + '/Tc1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            srcPath2=srcDir + '/Tc2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            if not os.path.exists(srcPath1):
                print 'No file',srcPath1
                continue

            if ch=='H':
                atc1 = np.load(srcPath1)
                atc2 = np.load(srcPath2)
                atc  = np.c_[atc1, atc2]

            elif ch=='L':
                atc  = np.load(srcPath1)
            else:
                print 'check ch',ch
                sys.exit()

            if atc.shape[0]==0:
                continue

            if idydx ==0:
                nrec, nc = atc.shape
                a4tcTmp = np.empty([nrec, ny, nx, nc], 'float32')

            a4tcTmp[:,cy+dy,cx+dx,:] = atc


        if a4tcTmp is None:
            continue
        else:
            a4tcTmp = np.array(a4tcTmp)

        #**********************
        # Resample
        #**********************
        if samplerate is not None:
            np.random.seed(0)  # Do not change !!
            aidx = np.random.choice(range(a4tcTmp.shape[0]), int(a4tcTmp.shape[0]*samplerate), replace=False)

            a4tcTmp = a4tcTmp[a1idx,:] 
        #**********************
        a4tc.extend(a4tcTmp)

    return np.array(a4tc)


def read_Tc_2d(lDTime=None, ldydx=None, isurf=None, samplerate=None, ch='LH'):
    a2tc = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        a2tcTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            #srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcDir = stopbaseDir + '/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcPath1=srcDir + '/Tc1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            srcPath2=srcDir + '/Tc2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            if not os.path.exists(srcPath1):
                print 'No file',srcPath1
                continue

            if ch=='H':
                atc1 = np.load(srcPath1)
                atc2 = np.load(srcPath2)
                atc  = np.c_[atc1, atc2]
            elif ch=='L':
                atc  = np.load(srcPath1)
            else:
                print 'check ch',ch
                sys.exit()

            if atc.shape[0]==0:
                continue

            if a2tcTmp is None:
                a2tcTmp = atc
            else:
                a2tcTmp = np.c_[a2tcTmp, atc]

        if a2tcTmp is None:
            continue
        else:
            a2tcTmp = np.array(a2tcTmp)

        #**********************
        # Resample
        #**********************
        if samplerate is not None:
            np.random.seed(0)  # Do not change !!
            aidx = np.random.choice(range(a2tcTmp.shape[0]), int(a2tcTmp.shape[0]*samplerate), replace=False)

            a2tcTmp = a2tcTmp[aidx,:]            
        #**********************
        a2tc.extend(a2tcTmp)

    return np.array(a2tc)

def read_var_collect(varName=None, lDTime=None, ldydx=None, isurf=None, samplerate=None):
    a2var = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        a2varTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            #srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcDir = stopbaseDir + '/data/%s/%04d/%02d/%02d'%(varName,Year,Mon,Day)
            srcPath=srcDir + '/%s.%ddy.%ddx.%02dsurf.npy'%(varName,dy,dx,isurf)
            if not os.path.exists(srcPath):
                print 'No file',srcPath
                continue
            avar = np.load(srcPath)

            if avar.shape[0]==0:
                continue

            if a2varTmp is None:
                a2varTmp = avar
            else:
                a2varTmp = np.c_[a2varTmp, avar]

        #**********************
        # Resample
        #**********************
        if samplerate is not None:
            np.random.seed(0)  # Do not change !!
            aidx = np.random.choice(range(a2varTmp.shape[0]), int(a2varTmp.shape[0]*samplerate), replace=False)

            a2varTmp = a2varTmp[aidx] 
        #**********************


        if a2varTmp is None:
            continue
        else:
            a2varTmp = np.array(a2varTmp)
        #**********************
        a2var.extend(a2varTmp)
    return np.array(a2var)


def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)

#def unit(x):
#    return ( x - np.min(x,0) )/( np.max(x,0) - np.min(x,0) )




def load_data_2d(lDTime, shuffleflag=True, samplerate=None):
    dtrainx = {}
    lvar    = []
    #ldx   = [-2,-1,0,1,2]
    imid  = int((len(ldy)*len(ldx)-1)/2)
    ldydx = [[dy,dx] for dy in ldy for dx in ldx]


    #****************************************************
    # Read data
    #****************************************************
    if 'H' in act:  ch= 'H'
    elif 'L' in act: ch= 'L'
    else: print 'check act',act; sys.exit()

    trainTc   = read_Tc_2d(lDTime, ldydx, isurf, ch=ch, samplerate=samplerate)
    trainStop = read_var_collect('stop', lDTime, [[0,0]], isurf, samplerate=samplerate)
   
  
    if 'T' in act:
        lvar.append('t2m')
        dtrainx['t2m'] = read_var_collect('t2m', lDTime, [[0,0]], isurf, samplerate=samplerate)
   
    if 'Q' in act:
        lvar.append('tqv')
        dtrainx['tqv'] = read_var_collect('tqv', lDTime, [[0,0]], isurf, samplerate=samplerate)

    if 'Z' in act:
        lvar.append('gtopo')
        dtrainx['gtopo'] = read_var_collect('tqv', lDTime, [[0,0]], isurf, samplerate=samplerate)

    #****************************************************
    # Screen invalid data
    #****************************************************
    #a1flag = ma.masked_inside(trainTc, tcmin, tcmax).all(axis=1).mask
    a1flag = ma.masked_inside(trainTc, tcmin, tcmax).mask.all(axis=1)
    #a1flag = a1flag * ~ma.masked_invalid(trainTc).all(axis=1).mask
    a1flag = a1flag * ~ma.masked_invalid(trainTc).mask.any(axis=1)
    a1flag = a1flag * ~ma.masked_invalid(trainStop).mask

    
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        
        a1flag = a1flag * ma.masked_inside(dtrainx[var],varmin, varmax).mask
        a1flag = a1flag * ~ma.masked_invalid(dtrainx[var]).mask
        
    
    trainTc  = trainTc  [a1flag]
    trainStop= trainStop[a1flag]
  
    for var in lvar:
        dtrainx[var] = dtrainx[var][a1flag]
 
    #print 'After Tc screening'
    #print trainTc.shape, trainStop.shape
    #print 'trainTc.min, max=',trainTc.min(), trainTc.max()
    #print 'stop.min, max=', trainStop.min(), trainStop.max()
    
    
    #***********************************************************
    # Normalize
    #***********************************************************
    
    trainTc = my_unit(trainTc, tcmin, tcmax)          
    
    #trainStop = my_unit(trainStop, stopmin, stopmax)
    
     
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        dtrainx[var] = my_unit(dtrainx[var], varmin, varmax)
        #print var, trainTc.shape,dtrainx[var].shape
    
    nrec = trainStop.shape[0]
    for var in lvar:
        trainTc = np.hstack([trainTc, dtrainx[var].reshape(-1,1)])
    #print trainTc.shape        

    if shuffleflag is True:
        a1idx = np.arange(trainTc.shape[0]).astype('int32')
        np.random.shuffle(a1idx)

        trainTc   = trainTc[a1idx]
        trainStop = trainStop[a1idx]


    return trainTc, trainStop





def load_data_4d(lDTime):
    dtrainx = {}
    lvar    = []
    #ldx   = [-2,-1,0,1,2]
    imid  = int((len(ldy)*len(ldx)-1)/2)
    ldydx = [[dy,dx] for dy in ldy for dx in ldx]


    #****************************************************
    # Read data
    #****************************************************
    if 'H' in act:  ch= 'H'
    elif 'L' in act: ch= 'L'
    else: print 'check act',act; sys.exit()


    #trainTc   = read_Tc(lDTime, ldydx, isurf, ch=ch)
    trainTc   = read_Tc_4d(lDTime, ldydx, isurf, ch=ch)
    
    
    trainStop = read_var_collect('stop', lDTime, [[0,0]], isurf)
   

   
    if 'T' in act:
        lvar.append('t2m')
        dtrainx['t2m'] = read_var_collect('t2m', lDTime, [[0,0]], isurf)
   
    if 'Q' in act:
        lvar.append('tqv')
        dtrainx['tqv'] = read_var_collect('tqv', lDTime, [[0,0]], isurf)

    if 'Z' in act:
        lvar.append('gtopo')
        dtrainx['gtopo'] = read_var_collect('tqv', lDTime, [[0,0]], isurf)

    print ''
    print trainTc.shape,trainStop.shape 
    print ''
 
 
    
    #****************************************************
    # Screen invalid data
    #****************************************************
    #a1flag = ma.masked_inside(trainTc, tcmin, tcmax).all(axis=1).mask
    #
    #a1flag = a1flag * ~ma.masked_invalid(trainTc).all(axis=1).mask
    #
    #a1flag = a1flag * ~ma.masked_invalid(trainStop).mask

    a1flag = ma.masked_inside(trainTc, tcmin, tcmax).all(axis=(1,2,3)).mask
   
    print 'a',a1flag.shape 
    a1flag = a1flag * ~ma.masked_invalid(trainTc).all(axis=(1,2,3)).mask
    
    print 'b',a1flag.shape 
    print 'stop',trainStop.shape
    a1flag = a1flag * ~ma.masked_invalid(trainStop).mask
 

    
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        
        a1flag = a1flag * ma.masked_inside(dtrainx[var],varmin, varmax).mask
        a1flag = a1flag * ~ma.masked_invalid(dtrainx[var]).mask
        
    
    trainTc  = trainTc  [a1flag]
    trainStop= trainStop[a1flag]
  
    for var in lvar:
        dtrainx[var] = dtrainx[var][a1flag]
 
    print 'After Tc screening'
    print trainTc.shape, trainStop.shape
    print 'trainTc.min, max=',trainTc.min(), trainTc.max()
    print 'stop.min, max=', trainStop.min(), trainStop.max()
    
    
    #***********************************************************
    # Normalize
    #***********************************************************
    
    trainTc = my_unit(trainTc, tcmin, tcmax)          
    
    #trainStop = my_unit(trainStop, stopmin, stopmax)
    
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        dtrainx[var] = my_unit(dtrainx[var], varmin, varmax)
        #print var, trainTc.shape,dtrainx[var].shape
    
    nrec = trainStop.shape[0]
    for var in lvar:
        #trainTc = np.hstack([trainTc, dtrainx[var].reshape(-1,1)])
        avar = np.ones([nrec,ny,nx,1]) * dtrainx[var].reshape(-1,1,1,1)
        trainTc = np.concatenate([trainTc, avar], axis=3)
        #print var, trainTc.shape, dtrainx[var].reshape(-1,1).shape    
    #print trainTc.shape        

    return trainTc, trainStop





def objective(trial):
    keras.backend.clear_session()

    ncol  = dncol[act]  
    model = build_model_2d_opt(trial)
    
    if savemodel==1:
        modelPath = cpDir + '/model.json'
        json_string = model.to_json()
        with open(modelPath, 'w') as f:
            f.write(json_string)

    if restmodel ==1: 
        latest = tf.train.latest_checkpoint(cpDir)
    
        if latest is not None:
            print 'Load weights'
            model.load_weights(latest)


    #************************
    if   season == 0:
        lDTime = util.ret_lDTime(datetime(Year,1,1), datetime(Year,12,31),timedelta(days=1))
    
    else:
        Mon  = season
        eDay = calendar.monthrange(Year,Mon)[1]
        lDTime = util.ret_lDTime(datetime(Year,Mon,1), datetime(Year,Mon,eDay), timedelta(days=1))


    lDTime = [x for x in lDTime if not x in lDTimeSkip]
    #************************


    if isurf==1:
        bsize = 10
        samplerate = 0.5
    else:
        bsize = 80
        samplerate = None
    np.random.shuffle(lDTime)
    llDTime = split2batchs(lDTime,bsize)


    for klDTime,lDTime in enumerate(llDTime):

        lDTime_train, lDTime_test = train_test_split(lDTime, train_size=0.8, test_size=0.2)
    
        trainX, trainY = load_data_2d(lDTime_train, shuffleflag=True, samplerate=samplerate)
        #****************************************************
        # Training
        #****************************************************
        # Set test data aside
    
        #cpPath= cpDir + '/cp-s%02d-{epoch:04d}.ckpt'%(isurf)
        #cpPath= cpDir + '/cp-s%02d-{val_loss:.2f}.ckpt'%(isurf)
        cpPath= cpDir + '/cp-s%02d.ckpt'%(isurf)
        #cp_callback = tf.keras.callbacks.ModelCheckpoint(cpPath, save_weights_only=True, verbose=1, monitor='val_loss', save_best_only=True, mode='min')
        cp_callback = TFKerasPruningCallback(trial, 'val_acc')
    
        early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
   
    
        if isurf ==1: 
            EPOCHS = 2
        elif isurf==3:
            EPOCHS = 5
        else:
            EPOCHS = 5

        print 'trainX.shape',trainX.shape
        print 'start fit'
        history = model.fit(
            x=trainX, y=trainY,
            batch_size=128,
            epochs=EPOCHS, validation_split=0.2,
            #callbacks = [cp_callback, early_stop],
            callbacks = [cp_callback],
            verbose=1)

        #histPath = cpDir + '/hist-%s-%s-s%02d.pickle'%(expr,act, isurf)
        #if restmodel ==0:
        #    if (klDTime==0)and(os.path.exists(histPath)):
        #        
        #        os.remove(histPath)

        #append_history(history, histPath)

    #score = history['val_loss']
    testX, testY = load_data_2d(lDTime_test)
    score = model.evaluate(testX, testY, verbose=0)[1]
    return score




print 'Define functions'
#***********************************************************
# Main loop start
#***********************************************************

for act in lact:
    cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)
    util.mk_dir(cpDir)
    
    for isurf in lisurf:
       
        ncol  = dncol[act]  

        study = optuna.create_study(
                direction='minimize',
                pruner=optuna.pruners.MedianPruner()
            )

        study.optimize(objective, n_trials=50)

        pruned_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.PRUNED]
        complete_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.COMPLETE]
        print 'Study statistics:'
        print  'Number of finished trials: ', len(study.trials)
        print  'Number of pruned trials: ', len(pruned_trials)
        print  'Number of complete trials: ', len(complete_trials)


        print "Best trial:"
        trial = study.best_trial
        print "  Value: ", trial.value
        print "  Params: "
        for key, value in trial.params.items():
            print 'key',key, 'value',value

        now = datetime.now()
        studyPath = cpDir + '/study-surf%d-'%(isurf)+ now.strftime('%Y%m%d_%H%M%S') + '.pickle'
    
        with open(studyPath, 'wb') as f:
            pickle.dump(study, f)
    
        print studyPath

 


