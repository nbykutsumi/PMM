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

def unit(x):
    return (x-np.min(x,0))/(np.max(x,0)-np.min(x,0))


#*********************************
# Prediction
#*********************************
ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-3,-2,-1,0,1,2,3]]
ldays_train,ldays_valid = mk_daylist(rat_train=0.8)
lDTime_train = [datetime(2017,1,1)+timedelta(days=i) for i in ldays_train]
lDTime_valid = [datetime(2017,1,1)+timedelta(days=i) for i in ldays_valid]
lDTime_valid = lDTime_valid[:30] # test
ntc1 = 9
ntc2 = 4
isurf = 3

a2egvec, a1varratio, a1cumvarratio = read_pc_coef(isurf)
npc_use = ma.masked_greater(a1cumvarratio,0.99).argmax()
ncomb = npc_use

ckptDir = '/work/hk01/utsumi/PMM/stop/ml-param'
ckptPath= ckptDir + '/stop.%02d'%(isurf)

#*********************************
# Read Test data
#*********************************
TesX, TesY = read_data(lDTime_valid, ldydx=ldydx, isurf=isurf)
TesY = TesY.reshape(-1,1)
TesX, TesY = unit(TesX), unit(TesY)
#-----------------------------
NUM_CPU=2
NUM_THREADS=4
with tf.Session(config=tf.ConfigProto(
  inter_op_parallelism_threads=NUM_CPU
 ,intra_op_parallelism_threads=NUM_THREADS)) as sess:

    saver = tf.train.import_meta_graph(ckptPath + '.meta')
    saver.restore(sess, ckptPath)


