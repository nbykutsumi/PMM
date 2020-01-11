from numpy import *
import numpy as np

axis =2
x = np.arange(2*2*3).reshape(2,2,3)
y = x + 1

xnumprof = np.ones([2,2,3]).astype('int32') 

xnumprof[0,0,0] = 0.5
xnumprof[0,1,:] = 0.5
ynumprof = xnumprof

ny,nx,nz = x.shape
xnummax = xnumprof.max(axis=2).reshape(ny,nx,1)
ynummax = ynumprof.max(axis=2).reshape(ny,nx,1)
a3xnumprof = ma.masked_invalid(xnumprof / xnummax )
a3ynumprof = ma.masked_invalid(ynumprof / ynummax )

x = ma.masked_where(a3xnumprof<0.9, x)
y = ma.masked_where(a3ynumprof<0.9, y)

xm = x.mean(axis=axis).reshape(ny,nx,1)
ym = y.mean(axis=axis).reshape(ny,nx,1)
A  = ((x-xm)*(y-ym)).sum(axis=axis)
B  = ((x-xm)**2).sum(axis=axis)
C  = ((y-ym)**2).sum(axis=axis)


CC= A/( np.sqrt(B*C) )

print 'x',x
print ''
print 'xm',xm
print ''
print 'A',A
print ''
print CC
