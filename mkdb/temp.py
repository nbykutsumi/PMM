import matplotlib
from numpy import *
from f_match_fov import *

a3in = arange(2*3*5).reshape(2,3,5)
a2x  = c_[ones(3).reshape(3,-1)*2,ones(3).reshape(3,-1)*4]
a2y  = c_[arange(3).reshape(3,-1), arange(3).reshape(3,-1)]


a3out= f_match_fov.extract_3d(a3in.T, a2x.T, a2y.T, -9999, -9999.).T

print '--- in ---'
print a3in
print '--- idx --'
print a2x
print a2y
print '--- out --'
print a3out




#a2in = arange(25).reshape(5,5).astype(int32)
#a2x  = c_[ones(5).reshape(5,-1)*2,ones(5).reshape(5,-1)*2]
#a2y  = c_[arange(5).reshape(5,-1), arange(5).reshape(5,-1)]
#
#a2out= f_match_fov.extract_2d(a2in.T, a2x.T, a2y.T, -9999, -9999.).T
#print '--- in ---'
#print a2in
#print '--- idx --'
#print a2x
#print a2y
#print '--- out --'
#print a2out
