import numpy as np
from datetime import datetime, timedelta
import myfunc.util as util
import os, sys, socket

#myhost = socket.gethostname()
#if myhost =='shui':
#    tankDir   = '/tank'
#
#elif myhost == 'well':
#    tankDir   = '/home/utsumi

lDTimeSkip = util.ret_lDTime(datetime(2017,9,25),datetime(2017,9,29),timedelta(days=1))

Year = 2017
iDTime = datetime(Year,1,1)
eDTime = datetime(Year,12,31)
lDTime = util.ret_lDTime(iDTime, eDTime, timedelta(days=1))
lDTime = [x for x in lDTime if not x in lDTimeSkip]

trainfrac = 0.8
nrec = len(lDTime)
ntrain= int(nrec * trainfrac)

np.random.seed(0)
np.random.shuffle(lDTime)
lDTime_train = lDTime[:ntrain]
lDTime_test  = lDTime[ntrain:]

print len(lDTime_train),len(lDTime_test)

trainPath = '/home/utsumi/bin/PMM/stop/ldtime-%04d-train.npy'%(Year)
testPath = '/home/utsumi/bin/PMM/stop/ldtime-%04d-test.npy'%(Year)

np.save(trainPath, lDTime_train)
np.save(testPath, lDTime_test)
print testPath

