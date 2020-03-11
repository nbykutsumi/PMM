# %%
%matplotlib inline
import numpy as np
from numpy import *
import glob
import myfunc.util as util
from datetime import datetime, timedelta
import JPLDB
import pandas as pd
import matplotlib.pyplot as plt

sensor = 'SSMIS'
#db_idx = 5810
db_idx = 5815
db = JPLDB.JPLDB(sensor)
srcDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/AMSR2_EPC_DATABASE'
srcPth = srcDir + '/db_%05d.bin'%(db_idx)

db.set_file(srcPth)
amn  = db.get_var('mn')
atb   = db.get_var('tb')
atb_min= db.get_var('tb_min')
atb_max= db.get_var('tb_max')
#print atb_min[:3]
#print ''
#print atb_max[:3]
asfc  = db.get_var('sfc_class')
aelev = db.get_var('elev')

#print asfc
#print asfc.min(), asfc.max()
#print ''
#print aelev
#print aelev.min(), aelev.max()
avfrac_conv_ns   =db.get_var('vfrac_conv_NS')
avfrac_conv_ms   =db.get_var('vfrac_conv_MS')
#print ''
#print avfrac_conv_ns
#print avfrac_conv_ns.min(), avfrac_conv_ns.max()
#print ''
#print avfrac_conv_ms
#print avfrac_conv_ms.min(), avfrac_conv_ms.max()
#print ''


aprecip_ns_cmb       =db.get_var('precip_NS_cmb')
aprecip_max_ns_cmb   =db.get_var('precip_max_NS_cmb')
aprecip_ms_cmb       =db.get_var('precip_MS_cmb')
aprecip_max_ms_cmb   =db.get_var('precip_max_MS_cmb')
avfrac_conv_ns_cmb   =db.get_var('vfrac_conv_NS_cmb')
avfrac_conv_ms_cmb   =db.get_var('vfrac_conv_MS_cmb')
#print ''
#print avfrac_conv_ns_cmb
#print avfrac_conv_ns_cmb.min(), avfrac_conv_ns_cmb.max()
#print ''
#print avfrac_conv_ms_cmb
#print avfrac_conv_ms_cmb.min(), avfrac_conv_ms_cmb.max()
#print ''



atype_precip_ns  = db.get_var('type_precip_NS')
atype_precip_ms  = db.get_var('type_precip_MS')
ashallow_rain_ns = db.get_var('shallow_rain_NS')
ashallow_rain_ms = db.get_var('shallow_rain_MS')
astopmax1 = db.get_var('storm_height_max_ku')
astopmax2 = db.get_var('storm_height_max_ka')
ats  = db.get_var('ts')
at2m = db.get_var('t2m')
aprof= db.get_var('precip_water_prof_NS')

#print ''
#print at2m
#print at2m.min(), at2m.max()

print ''
print  aprof
print aprof.min(), aprof.max()

dshallow_rain_ns = pd.DataFrame(ashallow_rain_ns, columns=['shallow_ns_%i'%i for i in range(5)])
dshallow_rain_ms = pd.DataFrame(ashallow_rain_ms, columns=['shallow_ms_%i'%i for i in range(5)])
dtype_precip_ns =pd.DataFrame(atype_precip_ns, columns=['type_ns_%i'%i for i in range(3)])
dtype_precip_ms =pd.DataFrame(atype_precip_ms, columns=['type_ms_%i'%i for i in range(3)])

df = pd.DataFrame()
df = df.assign(mn=amn)

df = df.assign(vfarc_conv_ns=avfrac_conv_ns)
df = df.assign(vfarc_conv_ms=avfrac_conv_ms)
df = df.assign(precip_ns_cmb=aprecip_ns_cmb)
df = df.assign(precip_max_ns_cmb=aprecip_max_ns_cmb)
df = df.assign(precip_ms_cmb=aprecip_ms_cmb)
df = df.assign(precip_max_ms_cmb=aprecip_max_ms_cmb)

df = df.assign(vfarc_conv_ns_cmb=avfrac_conv_ns_cmb)
df = df.assign(vfarc_conv_ms_cmb=avfrac_conv_ms_cmb)
df = pd.concat([df, dtype_precip_ns],axis=1)
df = pd.concat([df, dshallow_rain_ns],axis=1)
df = pd.concat([df, dtype_precip_ms],axis=1)
df = pd.concat([df, dshallow_rain_ms],axis=1)
df = df.assign(stopmax1=astopmax1)
df = df.assign(stopmax2=astopmax2)
df = df.assign(ts=ats)
df = df.assign(t2m=at2m)
dprof = pd.DataFrame(aprof, columns=['col-%i'%i for i in range(88)])
df = pd.concat([df,dprof],axis=1)
df.index.name = 'i'

#dprof[dprof<0] = np.nan

# %%
