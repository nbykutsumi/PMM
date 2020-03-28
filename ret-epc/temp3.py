import matplotlib
matplotlib.use('Agg')
from numpy import *
import matplotlib.pyplot as plt
import numpy as np

path1='/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.relsurf01.minrec1000.maxrec10000/2014/06/01/vfracConvrad.001463.npy'
path2='/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.relsurf02.minrec1000.maxrec10000/2014/06/01/vfracConvrad.001463.npy'

a1 = ma.masked_less(np.load(path1),0)
a2 = ma.masked_less(np.load(path2),0)

plt.scatter(a1,a2)
plt.savefig('/home/utsumi/temp/ret/temp.png')
