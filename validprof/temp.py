import matplotlib
matplotlib.use('Agg')
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import os, sys

#-- EPC ---
radpath2= '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.relsurf01.minrec1000.maxrec10000/2014/06/01/vfracConvrad.001454.npy'
pmwpath2= '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.relsurf01.minrec1000.maxrec10000/2014/06/01/vfracConvpmw.001454.npy'

rad2 = ma.masked_less(np.load(radpath2),0)
pmw2 = ma.masked_less(np.load(pmwpath2),0)

plt.scatter(rad2,pmw2)
plt.savefig('/home/utsumi/temp/ret/temp.epc.png')

sys.exit()


#-- GPROF ---
radpath1= '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift/2014/06/01/vfracConvrad.001454.npy'
radpath2= '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift.test/2014/06/01/vfracConvrad.001454.npy'

pmwpath2= '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift.test/2014/06/01/vfracConvpmw.001454.npy'


rad1 = ma.masked_less(np.load(radpath1),0)
rad2 = ma.masked_less(np.load(radpath2),0)

pmw2 = ma.masked_less(np.load(pmwpath2),0)

plt.scatter(rad2,pmw2)
plt.savefig('/home/utsumi/temp/ret/temp.gprof.png')




