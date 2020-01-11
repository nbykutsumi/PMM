from numpy import *
import numpy as np
import pickle
import os, sys, socket, glob
import optuna

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

season = 0
#expr  = 'a01'
expr  = 'opt01'
#lact = ['H','L','LT']
#lact = ['L','LT']
act = 'LTQZ'
lisurf = [2,3,4,5,6,7,8,9,10,11,12,14]
#lisurf = [6]
d = {}
for isurf in lisurf:
    cpDir = stopbaseDir + '/cp/%s-%s-ssn%s'%(expr,act, season)

    #studyPath = cpDir + '/study-surf%d-'%(isurf)+ now.strftime('%Y%m%d_%H%M%S') + '.pickle'
    ssearch   = cpDir + '/study-surf%d-*'%(isurf)+ '.pickle'
    lstudyPath = glob.glob(ssearch)
    if len(lstudyPath) ==0:
        print 'No file'
        print ssearch
        continue
    
    studyPath = lstudyPath[0]
    with open(studyPath,'rb') as f:
        d[isurf] = pickle.load(f)

    #trials = d[isurf].trials
    #for i in range(len(trials)):
    #    nlayers= trials[i].params['nlayers']
    #    value = trials[i].value
    #    params = trials[i].params
    #    number = trials[i].number
    #    #print number,nlayers, value, params

    print ''
    print '--------- isurf=',isurf,'------------'    
    best = d[isurf].best_trial
    nlayers = best.params['nlayers']
    params  = best.params
    print isurf,params
