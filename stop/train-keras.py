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
from tensorflow.keras import layers

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
lMon   = [0]
#lisurf = [4]
lisurf = range(1,14+1)
coef_b = 5
saveprep = 1
savemodel= 1
restmodel = 0
expr  = 'a01'
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
tcmin, tcmax = 50, 350
stopmin = 0
stopmax = 32000

ldy   = [-1,0,1]
ldx   = [-3,-2,-1,0,1,2,3]

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
        pickle.dump(hist_new, histPath) 



def split2batchs(a,bsize):
    n=len(a)/bsize + 1
    return [a[i*bsize:(i+1)*bsize] for i in range(n)] 

def build_model(ncol):
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

def Figure(Label, Prediction, bins):
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

    figPath = figDir + '/train.%s.png'%(expr)
    pl.savefig(figPath)
    pl.clf()
    print figPath


def FFN(TraX, TraY, TesX, TesY, learning_rate, epochs, batch_size, dim, act): 
    #ckptDir = stopbaseDir + '/stop/ml-param-%d'%(act)
    ckptDir = stopbaseDir + '/ml-param/%s'%(expr)
    util.mk_dir(ckptDir)
    ckptPath= ckptDir + '/stop.%02d'%(isurf)


    fn1 = tf.nn.sigmoid
    fn2 = tf.nn.relu
    def fn3(x):
        return x/(1+np.abs(x))
    ac  = [fn1,fn3,fn1,fn3,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1] # number of entry = len(dim) - 2
    total_batch = int(len(TraX)/batch_size) + 1
    Xdata = [ TraX[i*batch_size:(i+1)*batch_size] for i in range(total_batch) ]
    Ydata = [ TraY[i*batch_size:(i+1)*batch_size] for i in range(total_batch) ]
    
    tf.reset_default_graph()
    X = tf.placeholder(tf.float32, [None, TraX.shape[1]], 'input')
    Y = tf.placeholder(tf.float32, [None, TraY.shape[1]])
        
    W = [ tf.Variable(tf.random_normal([dim[i], dim[i+1]]), name='w%d'%(i)) for i in range(len(dim) - 1) ]
    b = [ tf.Variable(tf.random_normal([dim[i+1]]), name='b%d'%(i))         for i in range(len(dim) - 1) ]
    A = [ X ]
    for i in range(len(dim) - 2):
        A.append(ac[i](tf.matmul(A[-1], W[i]) + b[i]))
    A.append(tf.matmul(A[-1], W[-1]) + b[-1])  
    if act == 0:
        cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1]) ))) 
    elif act == 1:
        cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1])*error_function(Y,label,100,12,1,coef_b) ))) 
    else:
        cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1])*my_error_func(Y, coef_poly )))) 

    gogo = tf.train.AdamOptimizer(learning_rate).minimize(cost)
    real = tf.placeholder(tf.float32, [None, TraY.shape[1]])
    pred = tf.placeholder(tf.float32, [None, TraY.shape[1]])
    rmse = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(real - pred))))
    prediction=tf.add(tf.matmul(A[-2], W[-1]), b[-1], name='pred')  # for Save
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())


    #*** Saver ********
    saver = tf.train.Saver(max_to_keep=3)
    #******************

    for epoch in range(epochs):    
        for i in range(total_batch):
            feed1 = {X:Xdata[i], Y:Ydata[i]}
            sess.run(gogo, feed_dict = feed1)
            training_error = sess.run(cost, feed_dict = feed1)
            prediction     = sess.run(A[-1], feed_dict = {X:TesX})
            test_error     = sess.run(rmse, feed_dict = {real:TesY, pred:prediction})
        if epoch % 10 == 0:    
            print('Training Error:',training_error,'and','Testing Error:', test_error)

            #**** Save **********************
            if savemodel ==1:
                sv = saver.save(sess, ckptPath)
            #******************************** 
            
    return prediction



def mk_coef_polyfit(a1obs, degree, coef_b):
    # a1obs must be in range of [0,1]
    nbins = 100
    a1bnd = np.arange(nbins+1).astype('float32')/(nbins)
    frequency,_ = np.histogram(a1obs, bins=a1bnd)
    g = (coef_b - 1)*(-frequency/float(frequency.max()) +1)+1
    x = 0.5*(a1bnd[:-1]+a1bnd[1:])
    coef = np.polyfit(x,g,deg=degree)  # highest degree coef first.
    print 'coef.shape=',coef.shape
    return coef[::-1]

def my_error_func(a1obs, coef):
    degree = len(coef)-1
    for i in range(degree+1):
        if i==0:
            y = coef[i]*(a1obs**i)
        else:
            y = y + coef[i]*(a1obs**i)
    return y


    
def error_function(x, data, bins, degree, Min, Max):
    vec       = ((data - np.min(data,0))/float(np.max(data,0)-np.min(data,0)))
    interval  = [ i/float(bins) for i in range(bins + 1)]
    frequency = np.array([ ((vec<=interval[i+1]).sum() - (vec<interval[i]).sum())/float(len(vec)) for i in range(bins) ])
    xx        = np.arange(bins)/float(bins - 1)
    mat       = np.concatenate([(xx**i).reshape(-1,1) for i in range(degree)], axis=1)
    #print mat
    #print frequency
    coef      = np.dot(np.linalg.inv(np.dot(mat.T,mat)), np.dot(mat.T, frequency))
    poly      = 1 - sum([coef[i]*(x**i) for i in range(degree)])
    values    = 1 - sum([coef[i]*(xx**i) for i in range(degree)])
    M, N      = np.max(values), np.min(values)
    return (Max - Min)/float(M - N)*(poly - N) + Min 
def uniformize(reduction, label):
    mat = np.concatenate([reduction, label], axis=1)
    temporary = []
    for i in range(mat.shape[1]):
        a = np.arange(len(mat)).reshape(-1,1)
        b = np.concatenate([a,mat[:,[i]]], axis=1)
        c = b[b[:,1].argsort()]
        c[:,1] = np.arange(len(mat))/(len(mat)-1)
        d = c[c[:,0].argsort()]
        temporary.append(d[:,1])
    input_data = (np.array(temporary).T)[:, :-1]
    target     = (np.array(temporary).T)[:,[-1]]
    return input_data, target
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
def f_act(x,label):
    degree = 12
    #y_val = np.sort(unit(label.flatten()))
    y_val = np.sort(my_unit(label.flatten(),stopmin,stopmax))
    X     = (np.arange(len(y_val))/float(len(y_val)-1) )
    mat   = np.concatenate([(X**i).reshape(-1,1) for i in range(degree)], axis=1)
    coef  = np.dot(np.linalg.inv(np.dot(mat.T,mat)), np.dot(mat.T, y_val))
    poly  = sum([coef[i]*(x**i) for i in range(degree)])
    return poly



def read_Tc(lDTime=None, ldydx=None, isurf=None, samplerate=None, ch='LH'):
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

            a2tcTmp = a2tcTmp[a1idx,:]            
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

            a2varTmp = a2varTmp[a1idx,:] 
        #**********************


        if a2varTmp is None:
            continue
        else:
            a2varTmp = np.array(a2varTmp)
        #**********************
        a2var.extend(a2varTmp)
    return np.array(a2var)

def read_pc_coef(isurf):
    #*********************************
    # Read PC coefficient
    #*********************************
    #coefDir = '/work/hk01/utsumi/PMM/stop/data/coef'
    coefDir = stopbaseDir + '/data/coef'
    egvecPath = coefDir + '/egvec.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    egvalPath = coefDir + '/egval.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    varratioPath = coefDir + '/varratio.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    
    a2egvec = np.load(egvecPath)  # (n-th, ncomb)
    a1varratio = np.load(varratioPath)
    a1cumvarratio= np.cumsum(a1varratio)
    return a2egvec, a1varratio, a1cumvarratio

def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)

#def unit(x):
#    return ( x - np.min(x,0) )/( np.max(x,0) - np.min(x,0) )


def load_data(lDTime):
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
    
    trainTc   = read_Tc(lDTime, ldydx, isurf, ch=ch)
    
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


    
    #****************************************************
    # Screen invalid data
    #****************************************************
    a1flag = ma.masked_inside(trainTc, tcmin, tcmax).all(axis=1).mask
    
    a1flag = a1flag * ~ma.masked_invalid(trainTc).all(axis=1).mask
    
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
    
    print lvar
     
    for var in lvar:
        varmin, varmax = dvarminmax[var]
        dtrainx[var] = my_unit(dtrainx[var], varmin, varmax)
        print var, trainTc.shape,dtrainx[var].shape
    
    for var in lvar:
        trainTc = np.hstack([trainTc, dtrainx[var].reshape(-1,1)])
        print var, trainTc.shape, dtrainx[var].reshape(-1,1).shape    
    print trainTc.shape        

    return trainTc, trainStop




print 'Define functions'
#***********************************************************
# Main loop start
#***********************************************************

for act in lact:
    for isurf in lisurf:
       
        ncol  = dncol[act]  
        model = build_model(ncol)

        for Mon in lMon:

            #************************
            if   Mon == 0:
                lDTime = util.ret_lDTime(datetime(Year,1,1), datetime(Year,12,31),timedelta(days=1))
            
            else:
                eDay = calendar.monthrange(Year,Mon)[1]
                lDTime = util.ret_lDTime(datetime(Year,Mon,1), datetime(Year,Mon,eDay), timedelta(days=1))

            np.random.shuffle(lDTime)

            if isurf==1:
                bsize = 10
            else:
                bsize = 30
            llDTime = split2batchs(lDTime,bsize)


            for lDTime in llDTime:
    
                lDTime_train, lDTime_test = train_test_split(lDTime, train_size=0.8, test_size=0.2)
    
                trainX, trainY = load_data(lDTime_train)
     
                #****************************************************
                # Training
                #****************************************************
                # Set test data aside
    
                cpDir = stopbaseDir + '/cp/%s-%s'%(expr,act)
                cpPath= cpDir + '/cp-s%02d-{epoch:04d}.ckpt'%(isurf)
                cp_callback = tf.keras.callbacks.ModelCheckpoint(cpPath, save_weights_only=True, verbose=1, period=5)
    
                early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
    
    
                latest = tf.train.latest_checkpoint(cpDir)
            
                if latest is not None:
                    model.load_weights(latest)
    
    
                EPOCHS = 3
                history = model.fit(
                    x=trainX, y=trainY,
                    batch_size=128,
                    epochs=EPOCHS, validation_split=0.2,
                    callbacks = [cp_callback, early_stop],
                    verbose=1)
                    

                histPath = cpDir + '/hist-%s-%s-s%02d.pickle'%(expr,act, isurf)
                append_history(history, histPath)


                sys.exit() 
   

 
                #*******************************************************
                # Figure
                #*******************************************************
                #def plot_history(history, figPath, key='mse'):
                #    fig = plt.figure(figsize=(16,10))
     
                #    val = plt.plot(history.epoch, history.history['val_'+key],
                #                   '--', label=' Val')
                #    plt.plot(history.epoch, history.history[key], color=val[0].get_color(),
                #             label=' Train')
                #    
                #    plt.xlabel('Epochs')
                #    plt.ylabel(key.replace('_',' ').title())
                #    plt.legend()
                #    
                #    plt.xlim([0,max(history.epoch)])
    
                #    plt.savefig(figPath)
                #    print figPath
              
                #figPath = figDir + '/plot-hist-s%02d.png'%(isurf)
                # 
                #plot_history(history, figPath)
    
    
                testX, testY = load_data(lDTime_test)
                pred = model.predict(testX)
    
                expr= 'a-%s.s-%02d'%(act,isurf)
                corr = np.corrcoef(testY.flatten(), pred.flatten())[0,1]
                print corr
                Figure(testY, pred, 50)
    
                #***************************************************

