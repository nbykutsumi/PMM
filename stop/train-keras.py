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

argvs = sys.argv
if len(argvs)>1:
    print argvs
    lsurf = argvs[1:]
else:
    #lsurf = ['ocean','seaice','vege','snow','swater','coast','siedge']
    #lsurf = ['vege']
    #lsurf = ['ocean']
    #lsurf = ['seaice','vege','snow','swater','coast','siedge']
    #lsurf = ['seaice','vege']
    #lsurf = ['snow','swater','coast','siedge']
    #lsurf = ['swater']    
    #lsurf = ['siedge']    
    lsurf = ['seaice']    
    #lsurf = range(1,14+1)
    #lsurf = range(2,14+1)
    #lsurf = range(2,5+1)
    #lsurf = range(6,9+1)
    #lsurf = range(10,13+1)
    #lsurf = range(14,14+1)
    #lsurf = [1]

lsurfTmp = []
for x in lsurf:
    if x.isdigit(): x = int(x)
    lsurfTmp.append(x)
lsurf = lsurfTmp

onlypred = 0
savemodel= 1
#restmodel = 1
restmodel = 0
#restmodel = 1
#branch = 1
branch = 0

nloop = 20

#expr  = 'a01'
#expr  = 'mape01'
#expr  = 'wtpdf02'
#expr  = 'wtpdf02'
#expr  = 'wtpdf03'
#expr  = 'wtpdf04'  # b=20
#expr  = 'wtpdf04'   # b=40
#expr  = 'wtmape01'
#expr  = 'wtmape02' # 4-layers
#expr  = 'wtmape03' # 3-layers
#expr  = 'wtmape03-02' # 3-layers
#expr  = 'wtmape03-03' # 3-layers
#expr  = 'wtmape03-04' # 3-layers
#expr  = 'wtmape04' # 3-layers, 3X7
#expr  = 'wtmape05-01' # 3-layers, 3X3, batch-norm
#expr  = 'wtmape05-02' # 3-layers, 3X3, batch-norm
#expr  = 'wtpdf05-01' # 3-layers, 3X3, batch-norm
#expr  = 'wtpdf05-02' # 3-layers, 3X3, batch-norm
#expr  = 'mse01-01' # 3-layers, 3X3, batch-norm, mse
#expr  = 'mse01-02' # 3-layers, 3X3, batch-norm, mse
#expr  = 'mse01-03' # 3-layers, 3X3, batch-norm, mse
#expr  = 'mse01-04' # 3-layers, 3X3, batch-norm, mse
#expr  = 'mse01-05' # 3-layers, 3X3, batch-norm, mse
#expr  = 'mse02-01' # 3-layers, 3X3, batch-norm, mse, (64,64)
#expr  = 'mse02-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  Adam (default)
#expr  = 'mse03-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  # Adam, clipnorm=1.
#expr  = 'mse04-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  # Adam, clipnorm=0.1
#expr  = 'mse05-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  # Adam, clipnorm=0.01
#exprOrg  = 'mse06-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  # Adam, clipnorm=1e-3
#exprOrg  = 'mse07-01' # 3-layers, 1X1, batch-norm, mse, (64,64,64)  # Adam, clipnorm=1e-3
#exprOrg  = 'mse08-01' # 3-layers, 1X1, batch-norm, mse, (64,64)  # Adam, clipnorm=1e-3
#exprOrg  = 'mse09-01' # 3-layers, 1X1, batch-norm, mse, (64,64,64)  # Adam, clipnorm=1e-3
exprOrg  = 'mse10-01' # 3-layers, 1X1, batch-norm, mse, (64,64,64)  # Adam, clipnorm=1e-3, standing water (swater) with EPOCH=5
#expr  = 'best01-01' # 3-layers, 3X3, batch-norm, mse
#expr  = 'ovs01'
#expr  = 'msle01'
#expr  = 'male01'
#lact = ['H','L','LT']
#lact = ['L','LT']
#lact = ['LTQZ']
lact = ['HTQZ']
#lact = ['L']


dparent={
    'ocean':['mse01-02','HTQZ','mse'],
    'seaice':['mse01-01','HTQZ','mse'],
    'vege':['mse01-02','HTQZ','mse'],
#    'snow':,
#    'swater':,
#    'coast':,
#    'siedge':
    }

daddvar = {'T':['t2m'], 'Q':['tqv'], 'Z':'gtopo'}
dvarminmax = {'t2m':[200,350],'tqv':[0,120],'gtopo':[-100,8800]}
tcmin, tcmax = 50, 350
stopmin = 0
stopmax = 32000
#clipnorm = 1.
#clipnorm = 0.01
clipnorm = 1e-3
ldy = [0]
#ldy   = [-1,0,1]
#ldx   = [-3,-2,-1,0,1,2,3]
#ldx   = [-1,0,1]
ldx = [0]
ny = len(ldy)
nx = len(ldx)
cy = (ny-1)/2
cx = (nx-1)/2

dncol ={'H'  : ny*nx*13
       ,'L'  : ny*nx*9
       ,'LT' : ny*nx*9 + 1
       ,'LTQ': ny*nx*9 + 2
       ,'LTZ': ny*nx*9 + 2
       ,'LTQZ':ny*nx*9 + 3
       ,'HTQ': ny*nx*13 + 2
       ,'HTZ': ny*nx*13 + 2
       ,'HTQZ':ny*nx*13 + 3

       } 

dnc   ={'H': 13
       ,'L': 9
       ,'LT':9 + 1
       ,'LTQ':9 + 2
       ,'LTZ':9 + 2
       ,'LTQZ':9 + 3
       } 


lDTimeSkip = util.ret_lDTime(datetime(2017,9,25),datetime(2017,9,29),timedelta(days=1))

#lunits = [96,96,64,64]
#lunits = [96,96,64]
lunits = [64,64,64]
#lunits = [64,64]

#dunits = {
#    1: [64,64,64],
#    2: [32,32,32],  # opt
#    3: [32,32,32],  # opt
#    4: [32,32,32],  # opt
#    5: [96,32,32],  # opt
#    6: [96,32,32],  # opt
#    7: [96,32,32],  # opt
#    8: [32,32],  # opt
#    9: [32,32,32,32],  # opt
#    10:[32,32,64,96],  # opt
#    11:[32,32,32,32],  # opt
#    12:[32,32],  # opt
#    13:[64,64,64],
#    14:[32,32,32,32],  # opt
#    }


#***********************************************************
# Functions
#***********************************************************
def ret_lisurf(surf):
    if surf in ['ocean','seaice','vege','snow','swater','coast','siedge']:
        if   surf=='ocean':  lisurf=[1]
        elif surf=='seaice': lisurf=[2]
        elif surf=='vege':   lisurf=[3,4,5,6,7]
        elif surf=='snow':   lisurf=[8,9,10,11]
        elif surf=='swater': lisurf=[12]
        elif surf=='coast':  lisurf=[13]
        elif surf=='siedge': lisurf=[14]
    elif surf in range(1,14+1): 
        lisurf = [surf]
    else:
        print 'check surf',surf
        sys.exit()

    return lisurf



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
    pl.title('surf=%s act=%s box=%dx%d'%(surf, act,  len(ldy),len(ldx)))
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



def read_Tc_2d(lDTime=None, ldydx=None, surf=None, samplerate=None, ch='LH'):
    print ''
    print 'read_Tc_2d', surf
    lisurf = ret_lisurf(surf)

    a2tc = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        for isurf in lisurf:
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


def read_var_collect(varName=None, lDTime=None, ldydx=None, surf=None, samplerate=None):
    print ''
    print 'read_var_collect', surf

    lisurf = ret_lisurf(surf)

    a2var = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        for isurf in lisurf:
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



def load_data_2d(lDTime, surf=None, shuffleflag=True, samplerate=None):
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

    trainTc   = read_Tc_2d(lDTime, ldydx, surf=surf, ch=ch, samplerate=samplerate)
    trainStop = read_var_collect('stop', lDTime, [[0,0]], surf=surf, samplerate=samplerate)
   
  
    if 'T' in act:
        lvar.append('t2m')
        dtrainx['t2m'] = read_var_collect('t2m', lDTime, [[0,0]], surf=surf, samplerate=samplerate)
   
    if 'Q' in act:
        lvar.append('tqv')
        dtrainx['tqv'] = read_var_collect('tqv', lDTime, [[0,0]], surf=surf, samplerate=samplerate)

    if 'Z' in act:
        lvar.append('gtopo')
        dtrainx['gtopo'] = read_var_collect('tqv', lDTime, [[0,0]], surf=surf, samplerate=samplerate)

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


def correct_imbalance(trainx, trainy, ybnds, a=0.5):
    ahist,_ = np.histogram(trainy, bins=ybnds)
    nmax    = ahist.max()
    anhist  = ahist / float(nmax)  # normalized

    a1idx = np.arange(len(trainy)).astype('int32')

    a1idxsample = np.array([])
    for i in range(len(ybnds)-1):
        if anhist[i] == 0: continue

        ymin,ymax = ybnds[i],ybnds[i+1]
        a1idxTmp = ma.masked_where(ma.masked_outside(trainy, ymin, ymax).mask, a1idx).compressed()

        #-- resampling index array --
        nsample = int((1-anhist[i])*nmax * a)
        a1idxTmp2 = np.random.choice(a1idxTmp, nsample, replace=True)
       
        a1idxsample = np.concatenate([a1idxsample, a1idxTmp2])

    a1idxsample = a1idxsample.astype('int32')
    np.random.shuffle(a1idxsample)
    trainx_add = trainx[a1idxsample]
    trainy_add = trainy[a1idxsample]

    trainx = np.concatenate([trainx, trainx_add], axis=0)
    trainy = np.concatenate([trainy, trainy_add], axis=0)
         
    return trainx, trainy

def split_data(X,Y, trainfrac=0.8):
    nrec = X.shape[0]
    ntrain = int(nrec*trainfrac)
    trainX = X[:ntrain]
    trainY = Y[:ntrain]
    validX = X[ntrain:]
    validY = Y[ntrain:]
    return trainX, trainY, validX, validY

print 'Define functions'
#***********************************************************
# Cost functions
#***********************************************************

def gauss(x, *params):

    num_func = int(len(params)/3)
    y_list = []
    for i in range(num_func):
        #Y = np.zeros_like(x)
        param_range = list(range(3*i,3*(i+1),1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
        y = amp * tf.exp( -((x - ctr)/wid)**2)
        y_list.append(y)

    #y_sum = np.zeros_like(x)
    #y_sum = tf.zeros(x.shape)
    for iy, y in enumerate(y_list):
        if iy==0:
            y_sum = y
        else:
            y_sum = y_sum + y

    y_sum = y_sum + params[-1]

    return y_sum




def loss_pdf(params, pdfmax, a=1, b=10):
    def loss_pdf_inner(y_true, y_pred):
        fy = gauss(y_true, *params)
        gy = (b-a)*(-fy / pdfmax + 1) + a
        #return tf.reduce_mean(tf.abs((y_pred - y_true)/y_true)*gy)
        return tf.reduce_mean(tf.abs(y_pred - y_true)*gy)

    return loss_pdf_inner

def loss_wtmape(params, pdfmax, a=1, b=10):
    def loss_wtmape_inner(y_true, y_pred):
        fy = gauss(y_true, *params)
        gy = (b-a)*(-fy / pdfmax + 1) + a
        return tf.reduce_mean(tf.abs((y_pred - y_true)/y_true)*gy)

    print 'wtmape'
    return loss_wtmape_inner

def loss_mean_abs_log_error(y_true, y_pred):
    return tf.reduce_mean(tf.abs( tf.math.log(y_pred+1) - tf.math.log(y_true+1) ))

def test_func(y_pred, params, pdfmax, a=1, b=5):
    fy = gauss(y_pred, *params)
    gy = (b-a)*(-fy / pdfmax + 1) + a
    return gy





#***********************************************************
# Read parameter
#***********************************************************
paramDir = stopbaseDir + '/pdfparam'
dpdfparam = {}
dpdfmax   = {}
for surf in lsurf:
    paramPath = paramDir + '/gauss.param.surf-%s.npy'%(surf)
    param = np.load(paramPath)
    dpdfparam[surf] = param[:-1]
    dpdfmax[surf]   = param[-1]
    print surf,dpdfparam[surf] 

##-- test ---
#x = np.arange(200,16000+1,200)
#params = dpdfparam[7]
#pdfmax = dpdfmax[7]
#y = gauss(x, *params)
#y2 = test_func(x, params, pdfmax)
#
#print np.array(y)
#fig = plt.figure()
#plt.plot(x, y, color='k')
#plt.plot(x, y2,color='r')
#plt.savefig('/home/utsumi/temp/stop/temp.png')
#sys.exit()
##----------
#***********************************************************
# Main loop start
#***********************************************************


for act in lact:
    for iloop in range(nloop):
        expr = exprOrg + '-%02d'%(iloop)
        for surf in lsurf:
            if branch ==0:
                print ''
                print ''
                print 'brancth=',branch
                print ''
                print ''
                cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)
                util.mk_dir(cpDir)
                cpPath= cpDir + '/cp-s-%s.ckpt'%(surf)
        
                #lunits= dunits[surf]
                ncol  = dncol[act]  
                #nc  = dnc[act]  
                model = build_model_2d(ncol, lunits)
                #model = build_model_cnn([ny,nx,nc])
                if savemodel==1:
                    modelPath = cpDir + '/model-s-%s.json'%(surf)
                    json_string = model.to_json()
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
    
            else: # branch==1
                if onlypred==1:
                    print 'check branch',brach
                    sys.exit()
    
                print ''
                print ''
                print 'brancth=',branch
                print ''
                print ''
    
                pexpr, pact, loss = dparent[surf]
                parentDir= stopbaseDir + '/cp/%s-%s-ssn0'%(pexpr, pact) 
                parmodelPath=parentDir + '/model-s-%s.json'%(surf)
                parcpPath   =parentDir + '/cp-s-%s.ckpt'%(surf)
    
                act = pact
    
                cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)
                cpPath= cpDir + '/cp-s-%s.ckpt'%(surf)
                modelPath = cpDir + '/model-s-%s.json'%(surf)
    
                if restmodel ==1:
                    with open(modelPath,'r') as f:
                        json_string = f.read()
                    model = model_from_json(json_string)
        
                    model.load_weights(parcpPath)
    
                else:
                    print ''
                    print ''
                    print 'Read model and weights'
                    print parmodelPath
                    print parcpPath
                    print ''
                    print ''
    
                    with open(parmodelPath,'r') as f:
                        json_string = f.read()
                    model = model_from_json(json_string)
                    
                    model.load_weights(parcpPath)
    
    
                optimizer = tf.keras.optimizers.Adam(clipnorm=clipnorm)
                model.compile(optimizer=optimizer
                              ,loss = loss
                              ,metrics=['mae']
                             )
    
            #************************
    
            if surf==1:
                bsize = 366
                EPOCHS = 2
                samplerate = 0.2
            elif surf=='ocean':
                bsize = 366
                EPOCHS = 2
                samplerate = 0.2
            elif surf=='vege':
                bsize = 366
                EPOCHS = 2
                samplerate = None
            elif surf=='swater':
                bsize = 366
                EPOCHS = 5
                samplerate = None
            elif surf=='seaice':
                bsize = 366
                EPOCHS = 5
                samplerate = None
            elif surf=='siedge':
                bsize = 366
                EPOCHS = 15
                samplerate = None


            else:
                bsize = 366
                EPOCHS = 2
                samplerate = None

  
    
            lDTime_train = np.load('/home/utsumi/bin/PMM/stop/ldtime-%04d-train.npy'%(Year), allow_pickle=True)
            lDTime_test  = np.load('/home/utsumi/bin/PMM/stop/ldtime-%04d-test.npy'%(Year), allow_pickle=True)
            #lDTime_train, lDTime_test = train_test_split(lDTime_all, train_size=0.8, test_size=0.2)
    
            llDTime = split2batchs(lDTime_train, bsize)
    
            #************************
            if onlypred==1:
                print ''
                print ''
                print '***********************'
                print 'Only prediction'
                print '***********************'
                print ''
                print ''
                if surf in [1, 'ocean', 'vege']: 
                    testsamplerate = 0.01
                else:
                    testsamplerate = None
                testX, testY = load_data_2d(lDTime_test, surf, testsamplerate)
                pred = model.predict(testX)
                print 'surf=',surf    
                #expr= 'a-%s.s-%02d'%(act,surf)
                corr = np.corrcoef(testY.flatten(), pred.flatten())[0,1]
                print corr
                figPath = figDir + '/train.%s-%s-ssn%s-surf-%s.png'%(expr,act,season, surf)
                Figure(testY, pred, 50, figPath)
                sys.exit()    
            #************************
    
    
            for klDTime,lDTime in enumerate(llDTime):
            
                trainX, trainY = load_data_2d(lDTime, surf=surf, samplerate=samplerate, shuffleflag=True)
        
                #trainX, trainY, validX, validY = split_data(trainX, trainY, trainfrac=0.8) 
        
                #histbef,_ = np.histogram(trainY, np.arange(0,26000,500))  # test
                #histaft,xbnds = np.histogram(np.log(trainY))  # test
                #
                print ''
                print 'surf=',surf,'shapes train',trainX.shape, trainY.shape
                print ''
        
                ###-- oversampling -------
                ##ybnds = range(0,6000+1,200)
                ##trainX, trainY = correct_imbalance(trainX, trainY, ybnds, a=0.5)
        
                ##print ''
                ##print 'surf=',surf,'shapes train',trainX.shape, trainY.shape
        
                ##histaft,_= np.histogram(trainY, np.arange(0,26000,500))  # test
        
                #fig = plt.figure(figsize=(10,5))
                #ax  = fig.add_axes([0.2,0.2,0.6,0.6])
                ##x = np.arange(0,26000,500)[:-1] + 250
                #x = 0.5*(xbnds[:-1] + xbnds[1:])
                ##ax.plot(x, histbef, color='k')
                #ax.plot(x, histaft, color='r')
                ##ax.set_xlim([0,15000])
                #plt.savefig('/home/utsumi/temp/stop/temp.hist.png')
                #plt.clf()
                #sys.exit() 
                ##********************
        
                #****************************************************
                # Training
                #****************************************************
                #cpPath= cpDir + '/cp-s%02d-{epoch:04d}.ckpt'%(surf)
                #cpPath= cpDir + '/cp-s%02d-{val_loss:.2f}.ckpt'%(surf)
                #cpPath= cpDir + '/cp-s%02d.ckpt'%(surf)
                cp_callback = tf.keras.callbacks.ModelCheckpoint(cpPath, save_weights_only=True, verbose=1, monitor='val_loss', save_best_only=True, mode='min')
            
                early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)
            
           
            
                history = model.fit(
                    x=trainX, y=trainY,
                    batch_size=128,
                    #epochs=EPOCHS, validation_data=(validX, validY),
                    epochs=EPOCHS, validation_split=0.2,
                    callbacks = [cp_callback, early_stop],
                    verbose=1)
        
                #history = model.fit(
                #    x=trainX, y=trainY,
                #    #batch_size=128,
                #    batch_size=256,
                #    epochs=EPOCHS, validation_split=0.2,
                #    callbacks = [cp_callback, early_stop],
                #    verbose=1)
                    
        
                histPath = cpDir + '/hist-%s-%s-s%s.pickle'%(expr,act, surf)
                if restmodel ==0:
                    if (iloop==0)and(os.path.exists(histPath)):
                        
                        os.remove(histPath)
        
                append_history(history, histPath)
        
        
                #*******************************************************
                # Figure
                #*******************************************************
                model.load_weights(cpPath)
    
                if surf in [1, 'ocean', 'vege']: 
                    testsamplerate = 0.2
                else:
                    testsamplerate = None
                testX, testY = load_data_2d(lDTime_test, surf, testsamplerate)
                pred = model.predict(testX)
                print 'surf=',surf    
                #expr= 'a-%s.s-%02d'%(act,surf)
                corr = np.corrcoef(testY.flatten(), pred.flatten())[0,1]
                print corr
                figPath = figDir + '/train.%s-%s-ssn%s-surf-%s.png'%(expr,act,season, surf)
                Figure(testY, pred, 50, figPath)
            
                    #***************************************************

