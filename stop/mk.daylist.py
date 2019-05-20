import numpy as np
rtrain = 0.8
nall   = 365
ntrain = int(nall*rtrain)
np.random.seed(0)

a1idx = range(nall)
a1idx = np.random.choice(a1idx, len(a1idx), replace=False)
a1idx_train = a1idx[:ntrain]
a1idx_valid = a1idx[ntrain:]
print len(a1idx_train)
print len(a1idx_valid)
