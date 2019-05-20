from numpy import *
import numpy as np
import tensorflow as tf
import sys, os
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
#*********************************
# Functions
#*********************************
def FFN(TraX, TraY, TesX, TesY, learning_rate, epochs, batch_size, dim, act): 
    fn1 = tf.nn.sigmoid
    fn2 = tf.nn.relu
    def fn3(x):
        return x/(1+np.abs(x))
    ac  = [fn1,fn3,fn1,fn3,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1] # number of entry = len(dim) - 2
    total_batch = int(len(TraX)/batch_size) + 1
    Xdata = [ TraX[i*batch_size:(i+1)*batch_size] for i in range(total_batch) ]
    Ydata = [ TraY[i*batch_size:(i+1)*batch_size] for i in range(total_batch) ]
    
    tf.reset_default_graph()
    X = tf.placeholder(tf.float32, [None, TraX.shape[1]])
    Y = tf.placeholder(tf.float32, [None, TraY.shape[1]])
        
    W = [ tf.Variable(tf.random_normal([dim[i], dim[i+1]])) for i in range(len(dim) - 1) ]
    b = [ tf.Variable(tf.random_normal([dim[i+1]]))         for i in range(len(dim) - 1) ]
    A = [ X ]
    for i in range(len(dim) - 2):
        A.append(ac[i](tf.matmul(A[-1], W[i]) + b[i]))
    A.append(tf.matmul(A[-1], W[-1]) + b[-1])  
    if act == 0:
        cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1]) ))) 
    if act == 1:
        cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1])*error_function(Y,label,100,12,1,5) ))) 
    gogo = tf.train.AdamOptimizer(learning_rate).minimize(cost)
    real = tf.placeholder(tf.float32, [None, TraY.shape[1]])
    pred = tf.placeholder(tf.float32, [None, TraY.shape[1]])
    rmse = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(real - pred))))
    sess = tf.Session()
    sess.run(tf.global_variables_initializer())
    for epoch in range(epochs):    
        for i in range(total_batch):
            feed1 = {X:Xdata[i], Y:Ydata[i]}
            sess.run(gogo, feed_dict = feed1)
            training_error = sess.run(cost, feed_dict = feed1)
            prediction     = sess.run(A[-1], feed_dict = {X:TesX})
            test_error     = sess.run(rmse, feed_dict = {real:TesY, pred:prediction})
        if epoch % int(epochs/5) == 0:    
            print('Training Error:',training_error,'and','Testing Error:', test_error)
    return prediction

def Figure(Label, Prediction, bins):
    recover_testY = (Max-Min)*Label.flatten()      + Min
    recover_pred  = (Max-Min)*Prediction.flatten() + Min
    pl.figure(figsize=(15,15))
    gs = gridspec.GridSpec(2,2, width_ratios=[1,1], height_ratios=[1,1])
    
    pl.subplot(gs[0,:])
    pl.plot(recover_testY/1000, c='r', label ='Observation')
    pl.plot(recover_pred /1000, c='b', label ='Prediction')
    pl.ylabel('height(km)')
    pl.legend()
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
    pl.plot(revised_interval/1000          , cumulative_number   , color='green', alpha=0.5, label='Prediction')    
    pl.fill_between(revised_interval/1000  , cumulative_number, 0, color='green', alpha=0.5)
    pl.plot(revised_interval1/1000         , cumulative_number1  , color='red'  , alpha=0.5 ,label='Observation')    
    pl.fill_between(revised_interval1/1000 ,cumulative_number1, 0, color='red'  , alpha=0.5)
    pl.ylabel('number of samples')
    pl.xlabel('height(km)')
    pl.legend() 
    pl.title('Distribution')
    pl.legend()
    
    pl.subplot(gs[3])
    pl.scatter(recover_testY/1000, recover_pred/1000,s=3)
    pl.plot(np.arange(18000)/1000,np.arange(18000)/1000,c='black',linestyle = ':')
    pl.axis([0,18,0,18])
    pl.xticks([0,5,10,15])
    pl.yticks([0,5,10,15])
    pl.xlabel('Observation(km)')
    pl.ylabel('Prediction(km)')
    pl.title('Correlation')
    pl.grid()

def error_function(x, data, bins, degree, Min, Max):
    vec       = ((data - np.min(data,0))/(np.max(data,0)-np.min(data,0)))
    interval  = [ i/bins for i in range(bins + 1)]
    frequency = np.array([ ((vec<=interval[i+1]).sum() - (vec<interval[i]).sum())/len(vec) for i in range(bins) ])
    xx        = np.arange(bins)/(bins - 1)
    mat       = np.concatenate([(xx**i).reshape(-1,1) for i in range(degree)], axis=1)
    coef      = np.dot(np.linalg.inv(np.dot(mat.T,mat)), np.dot(mat.T, frequency))
    poly      = 1 - sum([coef[i]*(x**i) for i in range(degree)])
    values    = 1 - sum([coef[i]*(xx**i) for i in range(degree)])
    M, N      = np.max(values), np.min(values)
    return (Max - Min)/(M - N)*(poly - N) + Min 

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
    return np.sqrt( ( ( ((Max-Min)*x+Min).flatten()-((Max-Min)*y+Min).flatten() )**2 ).mean() )

def cc(x,y):
    return np.corrcoef( x.flatten(), y.flatten() )[0,1]

def sort(x):
    return np.sort(x.flatten())

def unit(x):
    return (x-np.min(x,0))/(np.max(x,0)-np.min(x,0))

def f_act(x):
    degree = 12
    y_val = np.sort(unit(label.flatten()))
    X     = (np.arange(len(y_val))/(len(y_val)-1) )
    mat   = np.concatenate([(X**i).reshape(-1,1) for i in range(degree)], axis=1)
    coef  = np.dot(np.linalg.inv(np.dot(mat.T,mat)), np.dot(mat.T, y_val))
    poly  = sum([coef[i]*(x**i) for i in range(degree)])
    return poly


def mk_daylist(rat_train=0.8):
    nall   = 365
    ntrain = int(nall*rat_train)
    np.random.seed(0)
    
    a1idx = range(nall)
    a1idx = np.random.choice(a1idx, len(a1idx), replace=False)
    a1idx_train = a1idx[:ntrain]
    a1idx_valid = a1idx[ntrain:]

    return a1idx_train,a1idx_valid

def read_Tc(lDTime=None, ldydx=None, isurf=None):
    a2tc = deque([])

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        a2tcTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcPath1=srcDir + '/Tc1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            srcPath2=srcDir + '/Tc2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            if not os.path.exists(srcPath1): continue
            atc1 = np.load(srcPath1)
            atc2 = np.load(srcPath2)
            atc  = c_[atc1, atc2]

            try:
                a2tcTmp = c_[a2tcTmp, atc]
            except ValueError:
                a2tcTmp = atc

        if a2tcTmp is None:
            continue
        else:
            a2tcTmp = array(a2tcTmp)


        #**********************
        a2tc.extend(a2tcTmp)

    return array(a2tc)

def read_stop(lDTime=None, isurf=None):
    a1stop = deque([])

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        srcDir = '/work/hk01/utsumi/PMM/stop/data/stop/%04d/%02d/%02d'%(Year,Mon,Day)
        srcPath=srcDir + '/stop.%02dsurf.npy'%(isurf)
        if not os.path.exists(srcPath): continue
        astop = np.load(srcPath)
        a1stop.extend(astop)
    return array(a1stop)

def read_pc_coef(isurf):
    #*********************************
    # Read PC coefficient
    #*********************************
    coefDir = '/work/hk01/utsumi/PMM/stop/data/coef'
    egvecPath = coefDir + '/egvec.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    egvalPath = coefDir + '/egval.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    varratioPath = coefDir + '/varratio.%02dch.%03dpix.%02dsurf.npy'%(ntc1+ntc2, len(ldydx),isurf)
    
    a2egvec = np.load(egvecPath)  # (n-th, ncomb)
    a1varratio = np.load(varratioPath)
    a1cumvarratio= np.cumsum(a1varratio)
    return a2egvec, a1varratio, a1cumvarratio

def read_data(lDTime=None, ldydx=None, isurf=None):
    #*********************************
    # Read PC coefficient
    #*********************************
    a2egvec, a1varratio, a1cumvarratio = read_pc_coef(isurf)

    #*********************************
    # Read mean and std of Tc
    #*********************************
    meanDir = '/work/hk01/utsumi/PMM/stop/data/coef'
    meanPath1= meanDir + '/mean1.0dy.0dx.%02dsurf.npy'%(isurf)
    meanPath2= meanDir + '/mean2.0dy.0dx.%02dsurf.npy'%(isurf)
    stdPath1 = meanDir + '/std1.0dy.0dx.%02dsurf.npy'%(isurf)
    stdPath2 = meanDir + '/std2.0dy.0dx.%02dsurf.npy'%(isurf)
    
    amean1 = np.load(meanPath1)
    amean2 = np.load(meanPath2)
    
    astd1 = np.load(stdPath1)
    astd2 = np.load(stdPath2)
    
    amean = deque([])
    astd  = deque([])
    for (dy,dx) in ldydx:
        amean.extend(amean1)
        amean.extend(amean2)
        astd.extend(astd1)
        astd.extend(astd2)
    amean = array(amean)
    astd  = array(astd)
    
    #*********************************
    # Read Tc data
    #*********************************
    a2tc = read_Tc(lDTime, ldydx=ldydx, isurf=isurf)
    
    #*********************************
    # Read storm top height data
    #*********************************
    a1stop = read_stop(lDTime, isurf=isurf)
    
    #**********************
    # Screen invalid data
    #**********************
    try:
        a1flag = ma.masked_inside(a2tc,50,350).mask.all(axis=1)
        a2tc   = a2tc[a1flag,:]
        a1stop = a1stop[a1flag]
    except AxisError:
        pass
    
    #*********************************
    # Normalize Tc
    #*********************************
    a2tc = a2tc - amean
    
    #*********************************
    # Tc --> PC
    #*********************************
    a2pc = np.dot(a2tc, a2egvec[:npc_use,:].T)

    return a2pc, a1stop
#***********************************************************************
#lisurf = range(1,15+1)
ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-3,-2,-1,0,1,2,3]]

ldays_train,ldays_valid = mk_daylist(rat_train=0.8)
lDTime_train = [datetime(2017,1,1)+timedelta(days=i) for i in ldays_train]
lDTime_valid = [datetime(2017,1,1)+timedelta(days=i) for i in ldays_valid]
ntc1 = 9
ntc2 = 4
isurf = 3
act  = 0  # 0: Normal. 1: Apply Error function
learning_rate = 0.005
a2egvec, a1varratio, a1cumvarratio = read_pc_coef(isurf)
npc_use = ma.masked_greater(a1cumvarratio,0.99).argmax()

ncomb = npc_use
dim = [ncomb, 30, 30, 30,10, 1]


#*********************************
# FFN Model
#*********************************
fn1 = tf.nn.sigmoid
fn2 = tf.nn.relu
def fn3(x):
    return x/(1+np.abs(x))
ac  = [fn1,fn3,fn1,fn3,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1] # number of entry = len(dim) - 2

tf.reset_default_graph()
X = tf.placeholder(tf.float32, [None, ncomb])
Y = tf.placeholder(tf.float32, [None, 1])
    
W = [ tf.Variable(tf.random_normal([dim[i], dim[i+1]])) for i in range(len(dim) - 1) ]
b = [ tf.Variable(tf.random_normal([dim[i+1]]))         for i in range(len(dim) - 1) ]
A = [ X ]
for i in range(len(dim) - 2):
    A.append(ac[i](tf.matmul(A[-1], W[i]) + b[i]))
A.append(tf.matmul(A[-1], W[-1]) + b[-1])  
if act == 0:
    cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1]) ))) 
if act == 1:
    cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1])*error_function(Y,label,100,12,1,5) ))) 
gogo = tf.train.AdamOptimizer(learning_rate).minimize(cost)
real = tf.placeholder(tf.float32, [None, 1])
pred = tf.placeholder(tf.float32, [None, 1])
rmse = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(real - pred))))

#*********************************
# Read Test data
#*********************************
lDTime_valid = lDTime_valid[:30] # test
TesX, TesY = read_data(lDTime_valid, ldydx=ldydx, isurf=isurf) 
TesY = TesY.reshape(-1,1)
TesX, TesY = unit(TesX), unit(TesY)
#-----------------------------
#*********************************
# Strat session
#*********************************
epochs = 2
nbatch = 40
batchsize= len(lDTime_train)/nbatch
NUM_CPU=2
NUM_THREADS=4
saver = tf.train.Saver(max_to_keep=3)
#with tf.Session() as sess:
with tf.Session(config=tf.ConfigProto(
  inter_op_parallelism_threads=NUM_CPU
 ,intra_op_parallelism_threads=NUM_THREADS)) as sess:

    sess.run(tf.global_variables_initializer())
    for epoch in range(epochs):    
        for ibatch in range(nbatch):
            if ibatch ==nbatch-1:
                lDTime_trainTmp= lDTime_train[ibatch*batchsize:]
            else:
                lDTime_trainTmp= lDTime_train[ibatch*batchsize:(ibatch+1)*batchsize]
        
            TraX, TraY = read_data(lDTime_train, ldydx=ldydx, isurf=isurf) 
            TraY = TraY.reshape(-1,1)
            TraX = unit(TraX)
            TraY = unit(TraY)

            feed1 = {X:TraX, Y:TraY}
            sess.run(gogo, feed_dict = feed1)
            training_error = sess.run(cost, feed_dict = feed1)
            prediction     = sess.run(A[-1], feed_dict = {X:TesX})
            test_error     = sess.run(rmse, feed_dict = {real:TesY, pred:prediction})
            if ibatch % 3 == 0:    
                print 'epoch=%d  ibatch=%d'%(epoch,ibatch)
                print('Training Error:',training_error,'and','Testing Error:', test_error)
                saveDir = '/work/hk01/utsumi/PMM/stop/ml-param'
                savePath= saveDir + '/stop.%02d'%(isurf)
                sv = saver.save(sess, savePath)
                print sv



sys.exit()



label = a1stop_train  # for error function
Max, Min = np.max(label), np.min(label)

dim = [a2pc_train.shape[1], 30, 30, 30,10, a1stop_train.shape[1]]
rmse2, cc2 = [], []
for i in range(20):
    #prediction   = FFN(trainX, trainY, testX, testY, 0.005, 50, 1024*4, dim, 1)
    prediction   = FFN(a2pc_train, a1stop_train, a2pc_valid, a1stop_valid, 0.005, 50, 1024*4, dim, 0)
    rmse2.append(Rmse(prediction, a1stop_valid))
    cc2.append(cc(prediction, a1stop_valid))



#*********************************





train = np.concatenate([ np.loadtxt('train3_12345.txt' ) 
                        ])
label = np.concatenate([ np.loadtxt('label3_12345.txt' )
                         ]).reshape(-1,1)

train.shape, label.shape

cov = 1/len(train)*np.dot((train-train.mean(0)).T, train-train.mean(0))
u,s,v = np.linalg.svd(cov)
restriction = 39
print(s[:restriction].sum()/s.sum())
U = u[:,:restriction]
reduction = np.dot(train, U)
Max, Min = np.max(label), np.min(label)
print(reduction.shape)





'''
# Apply error function
Train, Label = unit(reduction), unit(label)
ntrain       = int(0.7*len(reduction))
trainX, trainY, testX, testY = Train[:ntrain], Label[:ntrain], Train[ntrain:], Label[ntrain:]
print(trainX.shape, trainY.shape, testX.shape, testY.shape)
dim = [trainX.shape[1], 30, 30, 30,10, trainY.shape[1]]
rmse2, cc2 = [], []
for i in range(20):
    prediction   = FFN(trainX, trainY, testX, testY, 0.005, 50, 1024*4, dim, 1)
    rmse2.append(Rmse(prediction, testY))
    cc2.append(cc(prediction, testY))



# Apply error function
Train, Label = unit(reduction), unit(label)
ntrain       = int(0.7*len(reduction))
trainX, trainY, testX, testY = Train[:ntrain], Label[:ntrain], Train[ntrain:], Label[ntrain:]
print(trainX.shape, trainY.shape, testX.shape, testY.shape)
dim = [trainX.shape[1], 30, 30, 30,10, trainY.shape[1]]
rmse2, cc2 = [], []
for i in range(20):
    prediction   = FFN(trainX, trainY, testX, testY, 0.005, 50, 1024*4, dim, 1)
    rmse2.append(Rmse(prediction, testY))
    cc2.append(cc(prediction, testY))

'''
