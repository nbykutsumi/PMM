from numpy import *
from gv_fsub import *
import numpy as np
import myfunc.IO.GPM as GPM
import myfunc.util as util
import os, sys


a2in = arange(4*4).reshape(4,4)

a1groundBin = [0,1,2,3]
nxout = 3

a2out= gv_fsub.extract_slice_clusterprof(a2in.T, a1groundBin, nxout).T

print a2in
print ''
print a1groundBin
print ''
print a2out
