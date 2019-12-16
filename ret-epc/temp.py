import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/org.AMSR2.smp1000/10007'

estPath = srcDir + '/nsurfNScmb.est.10007.npy'
obsPath = srcDir + '/nsurfNScmb.obs.10007.npy'

est = np.load(estPath)
obs = np.load(obsPath)

plt.scatter(obs, est)
plt.savefig('/home/utsumi/temp/ret/temp.png')
plt.clf()
