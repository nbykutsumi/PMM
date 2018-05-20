from numpy import *
from gv_fsub import *
import numpy as np

dx   =4
a2in = arange(3*5).reshape(3,5)

a2in[0,1] = -99
a1idxTmp = [0,1,2]

a1ave = gv_fsub.mean_slice_negativemask(a2in.T, a1idxTmp, dx)


print a2in
print a1ave
