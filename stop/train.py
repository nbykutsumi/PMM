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
#def rmse(x,y):
#    x = x.flatten()
#    y = y.flatten()
#    return np.sqrt((((x-y))**2).mean())

#def Rmse(x,y):
#    return np.sqrt( ( ( ((Max-Min)*x+Min).flatten()-((Max-Min)*y+Min).flatten() )**2 ).mean() )

def cc(x,y):
    return np.corrcoef( x.flatten(), y.flatten() )[0,1]

def sort(x):
    return np.sort(x.flatten())

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
   
    if a2tc.shape[0] !=0: 
        #*********************************
        # Normalize Tc
        #*********************************
        #a2tc = a2tc - amean
        a2tc = (a2tc - amean)/astd
        #*********************************
        # Tc --> PC
        #*********************************
        a2pc = np.dot(a2tc, a2egvec[:npc_use,:].T)

    else:
        a2pc = array([])
        a1stop=array([])

    return a2pc, a1stop

def unit(x,Min,Max):
    return (x-Min)/(Max-Min)

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

def error_func(a1obs, coef):
    degree = len(coef)-1
    for i in range(degree+1):
        if i==0:
            y = coef[i]*(a1obs**i)
        else:
            y = y + coef[i]*(a1obs**i)
    return y


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
#act  = 1  # 0: Normal. 1: Apply Error function
learning_rate = 0.005
a2egvec, a1varratio, a1cumvarratio = read_pc_coef(isurf)
npc_use = ma.masked_greater(a1cumvarratio,0.99).argmax()

ncomb = npc_use
dim = [ncomb, 30, 30, 30,10, 1]
MinStop, MaxStop = 0, 32000  # [m]
degree = 12  # Polyfit power degree
coef_b = 5

#*********************************
# Read Test data
#*********************************
lDTime_valid = lDTime_valid[:30] # test
TesX, TesY = read_data(lDTime_valid, ldydx=ldydx, isurf=isurf) 
TesY = TesY.reshape(-1,1)
TesY = unit(TesY, MinStop,MaxStop)

print 'TesX.min,max',TesX.min(), TesX.max()
print 'TesY.min,max',TesY.min(), TesY.max()
#sys.exit()
#*********************************
# Coefficient for error weight function
#*********************************
lDTime_pop = lDTime_train[:40] # test
PopX,PopY  = read_data(lDTime_pop, ldydx=ldydx, isurf=isurf) 
PopY = PopY.reshape(-1,1)
PopY = unit(PopY, MinStop, MaxStop)
coef_poly  = mk_coef_polyfit(PopY, degree, coef_b)


#*********************************
# FFN Model
#*********************************
fn1 = tf.nn.sigmoid
fn2 = tf.nn.relu
def fn3(x):
    return x/(1+np.abs(x))
ac  = [fn1,fn3,fn1,fn3,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1,fn1] # number of entry = len(dim) - 2

tf.reset_default_graph()
X = tf.placeholder(tf.float32, [None, ncomb],name='input')
Y = tf.placeholder(tf.float32, [None, 1])
    
W = [ tf.Variable(tf.random_normal([dim[i], dim[i+1]]),name='w%d'%(i)) for i in range(len(dim) - 1) ]
b = [ tf.Variable(tf.random_normal([dim[i+1]]),name='b%d'%(i))         for i in range(len(dim) - 1) ]
A = [ X ]
for i in range(len(dim) - 2):
    A.append(ac[i](tf.matmul(A[-1], W[i]) + b[i]))
A.append(tf.matmul(A[-1], W[-1]) + b[-1])  
if act == 0:
    cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1]) ))) 
if act == 1:
    cost = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(Y - A[-1])*error_func(Y,coef_poly) ))) 
gogo = tf.train.AdamOptimizer(learning_rate).minimize(cost)
real = tf.placeholder(tf.float32, [None, 1])
#pred = tf.placeholder(tf.float32, [None, 1])
pred = tf.add(tf.matmul(A[-2], W[-1]), b[-1], name='pred')
rmse = tf.sqrt(tf.reduce_mean(tf.reduce_mean(tf.square(real - pred))))

#*********************************
# Strat session
#*********************************
epochs = 20
#epochs = 1
batchsize= 1
nbatch = int(len(lDTime_train)/batchsize)
NUM_CPU=2
NUM_THREADS=4
saver = tf.train.Saver(max_to_keep=3)

with tf.Session(config=tf.ConfigProto(
  inter_op_parallelism_threads=NUM_CPU
 ,intra_op_parallelism_threads=NUM_THREADS)) as sess:

    sess.run(tf.global_variables_initializer())
    for epoch in range(epochs):    
        np.random.shuffle(lDTime_train)  # Shuffle

        for ibatch in range(nbatch):
            if ibatch ==nbatch-1:
                lDTime_trainTmp= lDTime_train[ibatch*batchsize:]
            else:
                lDTime_trainTmp= lDTime_train[ibatch*batchsize:(ibatch+1)*batchsize]

            TraX, TraY = read_data(lDTime_trainTmp, ldydx=ldydx, isurf=isurf) 
            if TraX.shape[0]==0: continue

            TraY = TraY.reshape(-1,1)
            #TraX = unit(TraX)
            
            TraY = unit(TraY, MinStop, MaxStop)

            feed1 = {X:TraX, Y:TraY}
            sess.run(gogo, feed_dict = feed1)
            training_error = sess.run(cost, feed_dict = feed1)
            #prediction     = sess.run(A[-1], feed_dict = {X:TesX})
            prediction     = sess.run(pred, feed_dict = {X:TesX})
            test_error     = sess.run(rmse, feed_dict = {real:TesY, pred:prediction})
            if ibatch % 4 == 0:    
                print 'act=%d epoch=%d  ibatch=%d'%(act, epoch,ibatch)
                print('Training Error:',training_error,'and','Testing Error:', test_error)
                saveDir = '/work/hk01/utsumi/PMM/stop/ml-param-%d'%(act)
                util.mk_dir(saveDir)
                savePath= saveDir + '/stop.%02d'%(isurf)
                sv = saver.save(sess, savePath)
                print sv




