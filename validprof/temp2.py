# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import JPLDB
from numpy import *
%matplotlib inline

idx_db = 5810

db = JPLDB.JPLDB('AMSR2')
dbdir = '/media/disk2/share/PMM/JPLDB/EPC_DB/AMSR2_EPC_DATABASE'
dbpath= dbdir + '/db_%05d.bin'%(idx_db)
db.set_file(dbpath)

a2profcmb = db.get_var('precip_water_prof_NS')
a2profgpr = db.get_var('prof_GPROF')

a1preccmb = db.get_var('precip_NS_cmb')
a1precgpr = db.get_var('precip_GPROF')

l=[]
for i in range(a1preccmb.shape[0]):
    if (a1preccmb[i]>2)&(a1precgpr[i]>2):
        l.append(i)

i=l[3]
fig = plt.figure()
a1profcmb = ma.masked_less(a2profcmb[i],0)
plt.plot(a1profcmb)
plt.show()

fig = plt.figure()
a1profgpr = ma.masked_less(a2profgpr[i],0)
plt.plot(a1profgpr)

plt.show()



# %%
